/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_Polygon2D_H_
#define model_Polygon2D_H_

#include <cfloat>
#include <tuple>
#include <map>
#include <Eigen/Dense>
#include <assert.h>
//#include <earcut.hpp>

namespace model
{
    
    struct Polygon2D
    {
        
        
        template<typename PointType>
        static double area(const std::vector<PointType>& points)
        {
            double A(0.0);
            for(size_t k=0;k<points.size();++k)
            {
                const size_t k1(k<points.size()-1? k+1 : 0);
                A+=(points[k](0)-points[0](0))*(points[k1](1)-points[k](1))-(points[k](1)-points[0](1))*(points[k1](0)-points[k](0));
            }
            return 0.5*A;
        }
        
//        template<typename PointType>
//        static double inside(const std::vector<PointType>& points, const PointType& test) 
//        {
//            double A(area(points));
//            for(size_t k=0;k<points.size();++k)
//            {
//                const size_t k1(k<points.size()-1? k+1 : 0);
//                A+=(points[k](0)-points[0](0))*(points[k1](1)-points[k](1))-(points[k](1)-points[0](1))*(points[k1](0)-points[k](0));
//            }
//            return 0.5*A;
//        }
        
        template<typename PointType>
        static double isLeft(const PointType& P0, const PointType& P1, const PointType& P2 )
        {// return type changed from int to double http://geomalgorithms.com/a03-_inclusion.html
            return ( (P1(0) - P0(0)) * (P2(1) - P0(1))
                    - (P2(0) -  P0(0)) * (P1(1) - P0(1)) );
        }
        
        template<typename PointType>
        static int windingNumber(const PointType& P, const std::vector<PointType>& V)
        {// http://geomalgorithms.com/a03-_inclusion.html
            int    wn = 0;    // the  winding number counter
            
            // loop through all edges of the polygon
            for (size_t i=0; i<V.size(); i++) {   // edge from V[i] to  V[i+1]
                const size_t j(i<V.size()-1? i+1 : 0);
                if (V[i](1) <= P(1)) {          // start y <= P(1)
                    if (V[j](1)  > P(1))      // an upward crossing
                        if (isLeft( V[i], V[j], P) > 0)  // P left of  edge
                            ++wn;            // have  a valid up intersect
                }
                else {                        // start y > P(1) (no test needed)
                    if (V[j](1)  <= P(1))     // a downward crossing
                        if (isLeft( V[i], V[j], P) < 0)  // P right of  edge
                            --wn;            // have  a valid down intersect
                }
            }
            return wn;
        }
        
    };
    
} /* namespace model */
#endif
