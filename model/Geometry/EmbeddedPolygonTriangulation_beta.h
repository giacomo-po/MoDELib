/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_EmbeddedPolygonTriangulation_H_
#define model_EmbeddedPolygonTriangulation_H_

#include <Eigen/Dense>
#include <triangle.hpp>

namespace model
{
    
    
    template<int dim,typename...VertexData>
    struct EmbeddedPolygonTriangulationVertex : public std::tuple<VertexData...>
    {
        
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef Eigen::Matrix<double,dim-1,1> VectorLowerDim;

        const VectorDim P;
        const VectorLowerDim PL;
        
        EmbeddedPolygonTriangulationVertex(const VectorDim& P_in,const VectorLowerDim& PL_in) :
        /* init */ P(P_in)
        /* init */,PL(PL_in)
        {
            
        }
        
    };
    
    template<int dim,typename...VertexData>
    struct EmbeddedPolygonTriangulation : public std::vector<EmbeddedPolygonTriangulationVertex<dim,VertexData...>>
    {
      
        typedef Eigen::Matrix<double,dim,1> VectorLowerDim;
        typedef Eigen::Matrix<double,dim-1,1> VectorLowerDim;
        
//        const MeshPlane& meshPlane;
//        
//        MeshPlaneEmbeddedTriangulation(const MeshPlane& meshPlane_in) :
//        /* init */ meshPlane(meshPlane_in)
//        {
//            for(const auto& x : meshPlane.meshIntersections)
//            {
////                std::shared_ptr<MeshPlaneTriangulationVertex<dim,Args...>> v(new MeshPlaneTriangulationVertex<dim,Args...>());
//                this->emplace_back(x->P0,meshPlane.localPosition(x->P0));
//            }
//        }
        
        void triangulate(const double& triangleSize)
        {
            
        }
        
    };
    
    
}
#endif
