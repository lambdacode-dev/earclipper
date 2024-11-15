/****************************************************************************************
 * Ear Clipper implementation by Xianming Chen on Nov 6, 2024 for quilter.ai job intervew
 * Use 64 bit fixed point arithmetic - enough for typical PCB board with 0.1 µm precision.
 *
 * Set this value to false to use flaoting point arithmetic.
 *    constexpr bool use_fixed_point_arithmetic = true; 
 *
 *. Tested on Mac OS 14.5, Ubuntu 24.04 and Docker Ubuntu 22.04, with g++ 13, g++ 11, and
 *  Apple clang 16. To build: 
 *    make earclipper
 *  or,
 *    CXXFLAGS="-std=c++17" make earclipper
 **/

#include <iostream>
#include <fstream>
#include <cassert>
#include <list>
#include <iomanip>
#include <unordered_set>

constexpr bool use_fixed_point_arithmetic = true; 
constexpr int scale = use_fixed_point_arithmetic ? 10'000'000 : 1;
template<bool T>
using NumType = std::conditional_t<T, int64_t, double>;
using Num = NumType<use_fixed_point_arithmetic>;

/****************************************************************************************
 * Some simple linear algebra utilities for 2D points/vectors
 */
struct Point {
    Num x = 0, y = 0;
};
using Vector = Point;

inline Vector operator-(Point const& a, Point const& b) {
    return {b.x - a.x, b.y - a.y};
}

inline Num cross_product(Vector const& a, Vector const& b) {
    return a.x*b.y - a.y*b.x;
}

inline Num triangle_area(Point const& a, Point const& b, Point const& c) {
    return cross_product(b-a, c-a);
}

inline bool in_close_interval(Num n, Num a, Num b) {
    return a <= b ?
        (n >= a && n <= b):
        (n >= b && n <= a);
}

// Consider point coincident with edge end points as NOT inside triangle.
// This allows for polygons with holes where holes are connected to outer polygon
// via conincident edges of opposite directions (e.g. square_disk.csv)
inline bool inside_triangle(Point const& v, Point const& a, Point const& b, Point const& c) {
    auto vab = triangle_area(v, a, b);
    if(vab == 0)
        return false;

    auto vbc = triangle_area(v, b, c);
    if(vbc == 0)
        return  false;

    if(vab > 0 != vbc > 0)
        return false;

    auto vca = triangle_area(v, c, a);
    if(vca == 0)
        return false;
    
    return (vbc > 0 == vca > 0);
}
/****************************************************************************************/


namespace std {
template<>
class hash<list<Point>::iterator> {
public:
    size_t operator()(list<Point>::iterator i) const {
        return hash<Point*>{}(&*i);
    }
};
}

inline std::ostream& operator <<(std::ostream& os, Point const& p) {
    return os << double(p.x)/scale << "," << double(p.y)/scale;
}

inline void read_from_file(char const* filename, std::list<Point>& points) {
    std::ifstream file(filename);
    while(file) {
        char comma;
        double x, y;
        file >> x >> comma >> y;
        if(file)
            points.push_back({Num(x*scale),Num(y*scale)});
    }
}

int main (int argc, char** argv) {
    using namespace std;
    if(argc != 2) {
        cerr << "Usage: " << argv[0] << " polygon_csv_filename\n";
        return 1;
    }

    list<Point> points;
    read_from_file(argv[1], points);
    assert(points.size() >= 3);

    // if given last point = first: remove to cirular iterate with next()
    if(points.begin()->x == points.rbegin()->x && points.begin()->y == points.rbegin()->y)
        points.pop_back();  

    auto next = [&points](std::list<Point>::iterator itr) {
        if(++itr == points.end()) 
            itr = points.begin();
        return itr;
    };
    auto prev = [&points](std::list<Point>::iterator itr) {
        if(itr == points.begin()) 
            itr = points.end();
        return --itr;
    };

    // total signed area from integrating the piecewise linear function defined by points
    auto integrate_polygon = [&next, &points]() { 
        Num total = 0;
        if(points.size() >= 3) {
            auto p1 = points.begin(), p0 = p1++;
            do {
                total += (p0->y + p1->y) * (p1->x - p0->x);
                p0 = next(p0);
                p1 = next(p1);
            } while(p0 != points.begin());
        }
        return -total;  // negate so positive area for ccw polygon
    };

    Num area_from_integral = integrate_polygon(), 
        area_from_triangulation = 0;

    unordered_set<list<Point>::iterator> eartip_points;
    unordered_set<list<Point>::iterator> concav_points;

    auto check_convex = [&points, &prev, &next, area_from_integral] (list<Point>::iterator p1) {
        auto p0 = prev(p1);
        auto p2 = next(p1);
        auto area = triangle_area(*p0, *p1, *p2); 
        return area != 0 && area > 0 == area_from_integral > 0;;
    };

    auto check_ear = [&concav_points, &points, &prev, &next] (list<Point>::iterator p1) {
        auto p0 = prev(p1);
        auto p2 = next(p1);
        for(auto v : concav_points) {
            if(inside_triangle(*v, *p0, *p1, *p2))
                return false;
        }
        return true;
    };

    // classify into concave and convex points.
    for(auto p1 = points.begin(); p1 != points.end(); ++p1) {
        if(check_convex(p1))
            eartip_points.insert(p1); // not ear yet
        else 
            concav_points.insert(p1); // count middle point of degenerate triangle
    }
    // filter out non eartip from convex points
    for(auto itr = eartip_points.begin(); itr != eartip_points.end(); ) {
        if(check_ear(*itr))
            ++itr;
        else
            itr = eartip_points.erase(itr);
    }

    // clip ear
    while(!eartip_points.empty() && points.size() >= 3) {
        auto itr = eartip_points.begin();
        auto p1 = *itr;
        auto p0 = prev(p1);
        auto p2 = next(p1);
        auto area = triangle_area(*p0, *p1, *p2);
        assert(area == 0 || check_convex(p1));
        if(area) {
            area_from_triangulation += area;
            cout << *p0 << endl << *p1 << endl << *p2 << endl << endl;
        }
        eartip_points.erase(itr);
        points.erase(p1);
        for(auto p : {p0, p2}) {
            if(check_convex(p)) {
                concav_points.erase(p);
                if(check_ear(p))
                    eartip_points.insert(p);
                else
                    eartip_points.erase(p);
            }
        }
    }
    assert( abs(area_from_triangulation - area_from_integral) <= (use_fixed_point_arithmetic ? 0 : 1e-5) );
    cout << "Using " << (use_fixed_point_arithmetic ? "fixed" : "floating") << " point arithmetic\n";
    cout << "area_from_integral      = " << std::fixed << std::setprecision(20) << abs(area_from_integral/double(scale)/scale/2.0) << endl; 
    cout << "area_from_triangulation = " << std::fixed << std::setprecision(20) << abs(area_from_triangulation/double(scale)/scale/2.0) << endl; 
    return 0;
}
