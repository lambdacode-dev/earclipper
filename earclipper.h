#ifndef EARCLIPPER_H
#define EARCLIPPER_H
/****************************************************************************************
 * Ear Clipper implementation by Xianming Chen 
 **/

#include "la2d.h"
#include "integrate_polygon.h"
#include <list>
#include <cassert>
#include <unordered_set>

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
        area_from_integral = integrate_polygon(points);
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
        std::cout << "Using " << (use_fixed_point_arithmetic ? "fixed" : "floating") << " point arithmetic\n";
        std::cout << "area_from_integral      = " << std::fixed << std::setprecision(20) << abs(area_from_integral/double(scale)/scale/2.0) << std::endl; 
        std::cout << "area_from_triangulation = " << std::fixed << std::setprecision(20) << abs(area_from_triangulation/double(scale)/scale/2.0) << std::endl; 
        assert( abs(area_from_triangulation - area_from_integral) <= (use_fixed_point_arithmetic ? 0 : epsilon) );
    }
};
#endif
