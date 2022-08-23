/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SutherlandHodgman_H_
#define model_SutherlandHodgman_H_

#include <cfloat>
#include <tuple>
#include <map>
#include <Eigen/Dense>
#include <assert.h>
//#include <earcut.hpp>

namespace model
{
    
    class SutherlandHodgman
    {
        
        
        template<typename PointType>
        static bool inside(const PointType& p, const PointType& p1, const PointType& p2)
        {// check if a point is on the LEFT side of an edge
            return (p2(1) - p1(1)) * p(0) + (p1(0) - p2(0)) * p(1) + (p2(0) * p1(1) - p1(0) * p2(1)) < 0;
        }
        
        
        template<typename PointType>
        static PointType intersection(const PointType& cp1, const PointType& cp2, const PointType& s, const PointType& e)
        {// calculate intersection point
            const PointType dc(cp1-cp2);
            const PointType dp(s-e);
            
            const double n1 = cp1(0) * cp2(1) - cp1(1) * cp2(0);
            const double n2 = s(0) * e(1) - s(1) * e(0);
            const double n3 = 1.0 / (dc(0) * dp(1) - dc(1) * dp(0));
            
            return { (n1 * dp(0) - n2 * dc(0)) * n3, (n1 * dp(1) - n2 * dc(1)) * n3 };
        }
        
        
    public:
        
        template<typename PointType>
        static std::vector<PointType> clip(const std::vector<PointType>& subjectPolygon, const std::vector<PointType>& clipPolygon)
        {// Sutherland-Hodgman clipping
            
            //   const int   N = 99; // clipped (new) polygon size
            const int   N = subjectPolygon.size()+clipPolygon.size(); // clipped (new) polygon size
            PointType cp1, cp2, s, e, inputPolygon[N], newPolygon[N];
            
            // copy subject polygon to new polygon and set its size
            for(size_t i = 0; i < subjectPolygon.size(); i++)
            {
                newPolygon[i] = subjectPolygon[i];
            }
            
            size_t newPolygonSize = subjectPolygon.size();
            
            for(size_t j = 0; j < clipPolygon.size(); j++)
            {
                // copy new polygon to input polygon & set counter to 0
                for(size_t k = 0; k < newPolygonSize; k++)
                {
                    inputPolygon[k] = newPolygon[k];
                }
                int counter = 0;
                
                // get clipping polygon edge
                cp1 = clipPolygon[j];
                cp2 = clipPolygon[(j + 1) % clipPolygon.size()];
                
                for(size_t i = 0; i < newPolygonSize; i++)
                {
                    // get subject polygon edge
                    s = inputPolygon[i];
                    e = inputPolygon[(i + 1) % newPolygonSize];
                    
                    if(inside(s, cp1, cp2) && inside(e, cp1, cp2))
                    {// Case 1: Both vertices are inside:
                        // Only the second vertex is added to the output list
                        newPolygon[counter++] = e;
                    }
                    else if(!inside(s, cp1, cp2) && inside(e, cp1, cp2))
                    {// Case 2: First vertex is outside while second one is inside:
                        // Both the point of intersection of the edge with the clip boundary
                        // and the second vertex are added to the output list
                        newPolygon[counter++] = intersection(cp1, cp2, s, e);
                        newPolygon[counter++] = e;
                    }
                    else if(inside(s, cp1, cp2) && !inside(e, cp1, cp2))
                    {// Case 3: First vertex is inside while second one is outside:
                        // Only the point of intersection of the edge with the clip boundary
                        // is added to the output list
                        newPolygon[counter++] = intersection(cp1, cp2, s, e);
                    }
                    else if(!inside(s, cp1, cp2) && !inside(e, cp1, cp2))
                    {// Case 4: Both vertices are outside
                        // No vertices are added to the output list
                    }
                }
                // set new polygon size
                newPolygonSize = counter;
            }
            return std::vector<PointType>(&newPolygon[0],&newPolygon[0]+newPolygonSize);
        }
        
    };
    
} /* namespace model */
#endif

//// check if a point is on the LEFT side of an edge
//template<typename PointType>
//bool inside(const PointType& p, const PointType& p1, const PointType& p2)
//{
//    return (p2.y - p1.y) * p.x + (p1.x - p2.x) * p.y + (p2.x * p1.y - p1.x * p2.y) < 0;
//}
//
//// calculate intersection point
//template<typename PointType>
//PointType intersection(const PointType& cp1, const PointType& cp2, const PointType& s, const PointType& e)
//{
//    point2D dc = { cp1.x - cp2.x, cp1.y - cp2.y };
//    point2D dp = { s.x - e.x, s.y - e.y };
//
//    float n1 = cp1.x * cp2.y - cp1.y * cp2.x;
//    float n2 = s.x * e.y - s.y * e.x;
//    float n3 = 1.0 / (dc.x * dp.y - dc.y * dp.x);
//
//    return { (n1 * dp.x - n2 * dc.x) * n3, (n1 * dp.y - n2 * dc.y) * n3 };
//}
//
//// Sutherland-Hodgman clipping
//template<typename PointType>
//std::vector<PointType> SutherlandHodgman(const std::vector<PointType>& subjectPolygon, const std::vector<PointType>& clipPolygon)
//{
//
//    //   const int   N = 99; // clipped (new) polygon size
//    const int   N = subjectPolygon.size()+clipPolygon.size(); // clipped (new) polygon size
//    point2D cp1, cp2, s, e, inputPolygon[N], newPolygon[N];
//
//    // copy subject polygon to new polygon and set its size
//    for(int i = 0; i < subjectPolygon.size(); i++)
//    {
//        newPolygon[i] = subjectPolygon[i];
//    }
//
//    int newPolygonSize = subjectPolygon.size();
//
//    for(int j = 0; j < clipPolygon.size(); j++)
//    {
//        // copy new polygon to input polygon & set counter to 0
//        for(int k = 0; k < newPolygonSize; k++)
//        {
//            inputPolygon[k] = newPolygon[k];
//        }
//        int counter = 0;
//
//        // get clipping polygon edge
//        cp1 = clipPolygon[j];
//        cp2 = clipPolygon[(j + 1) % clipPolygon.size()];
//
//        for(int i = 0; i < newPolygonSize; i++)
//        {
//            // get subject polygon edge
//            s = inputPolygon[i];
//            e = inputPolygon[(i + 1) % newPolygonSize];
//
//            if(inside(s, cp1, cp2) && inside(e, cp1, cp2))
//            {// Case 1: Both vertices are inside:
//                // Only the second vertex is added to the output list
//                newPolygon[counter++] = e;
//            }
//            else if(!inside(s, cp1, cp2) && inside(e, cp1, cp2))
//            {// Case 2: First vertex is outside while second one is inside:
//                // Both the point of intersection of the edge with the clip boundary
//                // and the second vertex are added to the output list
//                newPolygon[counter++] = intersection(cp1, cp2, s, e);
//                newPolygon[counter++] = e;
//            }
//            else if(inside(s, cp1, cp2) && !inside(e, cp1, cp2))
//            {// Case 3: First vertex is inside while second one is outside:
//                // Only the point of intersection of the edge with the clip boundary
//                // is added to the output list
//                newPolygon[counter++] = intersection(cp1, cp2, s, e);
//            }
//            else if(!inside(s, cp1, cp2) && !inside(e, cp1, cp2))
//            {// Case 4: Both vertices are outside
//                // No vertices are added to the output list
//            }
//        }
//        // set new polygon size
//        newPolygonSize = counter;
//    }
//    return std::vector<point2D>(&newPolygon[0],&newPolygon[0]+newPolygonSize);
//}
//
//void writePolygon(const std::vector<point2D>& poly, const std::string& filename)
//{
//    std::ofstream outputFile(filename);
//    for(const auto& pt : poly)
//    {
//        outputFile<<std::setprecision(15)<<std::scientific<<pt.x<<" "<<pt.y<<"\n";
//    }
//}
//
//int main(int argc, char** argv)
//{
//    // subject polygon
//    std::vector<point2D> subjectPolygon{
//        {50,150}, {200,50}, {350,150},
//        {350,300},{250,300},{200,250},
//        {150,350},{100,250},{100,200}
//    }; // must be counter-clockwise
//
//    // clipping polygon
//    std::vector<point2D> clipPolygon{ {100,100}, {300,100}, {300,300}, {100,300} }; // must be counter-clockwise
//
//    // The new clipped polygon
//    std::vector<point2D> newPolygon(SutherlandHodgman(subjectPolygon,clipPolygon));
//
//    // print clipped polygon points
//    cout << "Clipped polygon points:" << endl;
//    for(int i = 0; i < newPolygon.size(); i++)
//    {
//        cout << "(" << newPolygon[i].x << ", " << newPolygon[i].y << ")" << endl;
//    }
//
//    writePolygon(subjectPolygon,"subject.txt");
//    writePolygon(clipPolygon,"clipper.txt");
//    writePolygon(newPolygon,"result.txt");
//
//
//    return 0;
//}
