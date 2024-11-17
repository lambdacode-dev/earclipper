/****************************************************************************************
 * Use 64 bit fixed point arithmetic - enough for typical PCB board with 0.1 µm precision.
 *
 * Set this value to false to use floating point arithmetic.
 *    constexpr bool use_fixed_point_arithmetic = true; 
 **/
#include <iostream>
#include <fstream>
#include <iomanip>

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
template<typename List>
void read_from_file(char const* filename, List& points) {
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
