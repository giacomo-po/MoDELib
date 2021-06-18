/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_TriangularMesh_H_
#define model_TriangularMesh_H_

#include <deque>
#include <Eigen/Dense>
#include <iostream>
#include <triangle.h>

namespace model
{
	
    class TriangularMesh : public std::deque<Eigen::Vector2d>
    /*                  */,public std::deque<Eigen::Vector3i>
    
    {
        

        
        
    public:
        
        const Eigen::Vector2d& vertex(const size_t& vID);
        std::deque<Eigen::Vector2d>& vertices();
        const std::deque<Eigen::Vector2d>& vertices() const;
        std::deque<Eigen::Vector3i>& triangles();
        const std::deque<Eigen::Vector3i>& triangles() const;
        void reMesh(const std::deque<Eigen::Matrix<double,2,1>>& boundaryPts,
                    const std::deque<Eigen::Matrix<double,2,1>>& internalPts,
                    const double& meshSize);
        
        
    };
	
} // namespace model
#endif
