/****************************************************************************************
 * Ear Clipper implementation by Xianming Chen 
 * Use 64 bit fixed point arithmetic - enough for typical PCB board with 0.1 µm precision.
 *
 * Set this value to false to use floating point arithmetic.
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
constexpr Num epsilon = use_fixed_point_arithmetic ? 0 : Num(1e-8);

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
    auto area = cross_product(b-a, c-a);
    if(abs(area) <= epsilon)
        area = 0;
    return area;
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


/************** IO **********************************************************************/
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


class EarClipper {
    std::list<Point> points;
    using PointPtr = std::list<Point>::iterator;
    std::unordered_set<PointPtr> eartip_points;
    std::unordered_set<PointPtr> concav_points;

    Num area_from_integral = 0, area_from_triangulation = 0;

    // total signed area from integrating the piecewise linear function defined by points
    Num integrate_polygon () { 
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
    }
    bool check_convex(PointPtr p1) {
        auto p0 = prev(p1);
        auto p2 = next(p1);
        auto area = triangle_area(*p0, *p1, *p2); 
        return area != 0 && area > 0 == area_from_integral > 0;;
    }

    bool check_ear (PointPtr p1) {
        auto p0 = prev(p1);
        auto p2 = next(p1);
        for(auto v : concav_points) {
            if(inside_triangle(*v, *p0, *p1, *p2))
                return false;
        }
        return true;
    }

    void find_concave_and_eartips() {
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
    }

    // circular iterator
    PointPtr next(PointPtr itr) {
        if(++itr == points.end()) 
            itr = points.begin();
        return itr;
    }
    PointPtr prev(PointPtr itr) {
        if(itr == points.begin()) 
            itr = points.end();
        return --itr;
    }
public:
    EarClipper(std::list<Point>&& _points) : points(_points) {
        // if given last point = first: remove to cirular iterate with next()
        if(points.begin()->x == points.rbegin()->x && points.begin()->y == points.rbegin()->y)
            points.pop_back();  

        assert(points.size() >= 3);
        area_from_integral = integrate_polygon();
        find_concave_and_eartips();
    }
    void operator()() {
        while(!eartip_points.empty() && points.size() >= 3) {
            auto itr = eartip_points.begin();
            auto p1 = *itr;
            auto p0 = prev(p1);
            auto p2 = next(p1);
            auto area = triangle_area(*p0, *p1, *p2);
            assert(area == 0 || check_convex(p1));
            if(area) {
                area_from_triangulation += area;
                std::cout << *p0 << std::endl << *p1 << std::endl << *p2 << std::endl << std::endl;
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
        assert( abs(area_from_triangulation - area_from_integral) <= (use_fixed_point_arithmetic ? 0 : epsilon) );
        std::cout << "Using " << (use_fixed_point_arithmetic ? "fixed" : "floating") << " point arithmetic\n";
        std::cout << "area_from_integral      = " << std::fixed << std::setprecision(20) << abs(area_from_integral/double(scale)/scale/2.0) << std::endl; 
        std::cout << "area_from_triangulation = " << std::fixed << std::setprecision(20) << abs(area_from_triangulation/double(scale)/scale/2.0) << std::endl; 
    }
};

int main (int argc, char** argv) {
    using namespace std;
    if(argc != 2) {
        cerr << "Usage: " << argv[0] << " polygon_csv_filename\n";
        return 1;
    }

    list<Point> points;
    read_from_file(argv[1], points);
    EarClipper clipper(std::move(points));
    clipper();

    return 0;
}
