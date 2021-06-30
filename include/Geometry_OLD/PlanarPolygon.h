/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PlanarPolygon_H_
#define model_PlanarPolygon_H_

#include <cfloat>
#include <tuple>
#include <map>
#include <Eigen/Dense>
#include <assert.h>
#include <earcut.hpp>

namespace model
{
    
    class PlanarPolygon : public std::deque<std::deque<std::array<double, 2>>>
    {
        
        using Point2d = std::array<double, 2>;
        
        const Eigen::Matrix3d R;
        const Eigen::Vector3d n;
        
        /**********************************************************************/
        static Eigen::Matrix3d getR(const Eigen::Vector3d& x,
                                    const Eigen::Vector3d& z);
        
        /**********************************************************************/
        std::array<double, 2> projectRotate(const Eigen::Vector3d& v);
        
    public:
        
        
        PlanarPolygon(const Eigen::Vector3d& x,        // global coordinate of the x axis on the plane
                      const Eigen::Vector3d& z);      // global z coordinate
       
        
        /**********************************************************************/
        void assignPoints(const std::deque<Eigen::Vector3d>& points);
        
        /**********************************************************************/
        void addHole(const std::deque<Eigen::Vector3d>& hole);
        
        /**********************************************************************/
        std::deque<std::array<size_t, 3>> triangulate() const;
        
//        bool isInside(const Eigen::Vector3d& P)
//        {
//            return isInside(projectRotate(P))
//        }
//
//
//
//
//
//
//        bool isInside(const std::array<double, 2>& P)
//        {
//
//            std::array<double, 2> ca;
//            Eigen::Map<Eigen::Vector2d> c(ca.data());
//            for(const auto& point)
//
//
//            return isInside(projectRotate(P))
//        }
        
        
        
        
    };
    
} /* namespace model */
#endif
