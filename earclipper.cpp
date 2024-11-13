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

    auto next = [&points](std::list<Point>::iterator& itr) {
        if(++itr == points.end()) 
            itr = points.begin();
        return itr;
    };

    // total signed area from integrating the piecewise linear function defined by points
    auto integrate_polygon = [&next, &points]() { 
        Num total = 0;
        if(points.size() >= 3) {
            auto p1 = points.begin(), p0 = p1++;
            do {
                total += (p0->y + p1->y) * (p1->x - p0->x);
                next(p0);
                next(p1);
            } while(p0 != points.begin());
        }
        return -total;  // negate so positive area for ccw polygon
    };

    Num area_from_integral = integrate_polygon(), 
        area_from_triangulation = 0;

    while(points.size() >= 3) {
        auto p0 = points.begin();
        auto p1 = ++points.begin();
        auto p2 = ++++points.begin();
        do {
            auto area = triangle_area(*p0, *p1, *p2);
            bool degenerate = (area == 0);
            bool is_ear = !degenerate && (area > 0 == area_from_integral > 0); // p1 could be ear tip only if convex corner
            for(auto v = p2; is_ear && next(v) != p0; ) {
                is_ear = !inside_triangle(*v, *p0, *p1, *p2);
            }
            if(degenerate || is_ear) {
                if(is_ear) {
                    area_from_triangulation += area;
                    cout << *p0 << endl << *p1 << endl << *p2 << endl << endl;
                }
                points.erase(p1);
                break;
            }
            else {
                next(p0);
                next(p1);
                next(p2);
            }
        }
        while(p0 != points.begin());
    }
    assert( abs(area_from_triangulation - area_from_integral) <= (use_fixed_point_arithmetic ? 0 : 1e-5) );
    cout << "Using " << (use_fixed_point_arithmetic ? "fixed" : "floating") << " point arithmetic\n";
    cout << "area_from_integral      = " << std::fixed << std::setprecision(20) << abs(area_from_integral/double(scale)/scale/2.0) << endl; 
    cout << "area_from_triangulation = " << std::fixed << std::setprecision(20) << abs(area_from_triangulation/double(scale)/scale/2.0) << endl; 
    return 0;
}