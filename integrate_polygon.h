#ifndef INTEGRATE_POLYGON_H 
#define INTEGRATE_POLYGON_H  
#include "la2d.h"

// total signed area from integrating the piecewise linear function defined by points
template<typename PointList>
Num integrate_polygon (PointList const& points) { 
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
#endif
