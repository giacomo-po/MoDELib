/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_MidPlane_H_
#define model_MidPlane_H_

#include <vector>
#include <array>
#include <cfloat>
#include <Eigen/Dense>
#include <Simplex.h>
#include <IntegrationDomain.h>

namespace model
{
    template <int i>
    class MidPlane
    {
        
        const double mid;
        const double tol;

    public:
        
        /**************************************/
        template <typename FiniteElementType>
        MidPlane(const FiniteElementType& fe,
               const double& tol_in=FLT_EPSILON) :
        /* init list */ mid(0.5*(fe.xMin()(i)+fe.xMax()(i))),
        /* init list */ tol(tol_in)
        {

        }
        
        /**************************************/
        template <typename NodeType>
        bool operator()(const NodeType& node) const
        {
            return std::abs(node.P0(i)-mid)<tol;
            
        }
        
//        /**************************************/
//        template <int dim>
//        bool operator()(const Eigen::Matrix<double,dim,1>& P0) const
//        {
//            return std::abs(P0(i)-min)<tol;
//            
//        }
//        
//        /**********************************************************************/
//        template <typename FiniteElementType, size_t qOrder, template <short unsigned int, size_t> class QuadratureRule>
//        static IntegrationDomain<FiniteElementType,1,qOrder,QuadratureRule> boundary(const FiniteElementType& fe)
//        {
//            
//            AtXmin<i> atx(fe);
//            
//            IntegrationDomain<FiniteElementType,1,qOrder,QuadratureRule> temp;
//            
//            
//            for (typename FiniteElementType::ElementContainerType::const_iterator eIter =fe.elementBegin();
//                 /*                                                            */ eIter!=fe.elementEnd();
//                 /*                                                            */ eIter++)
//            {
//                if(eIter->second.isBoundaryElement())
//                {
//                    const std::vector<int> boundaryFaces=eIter->second.boundaryFaces();
//                    for (int f=0;f<boundaryFaces.size();++f)
//                    {
//                        bool isExternalBoundaryFace(true);
//                        
//                        std::array<const Simplex<FiniteElementType::dim,0>*, SimplexTraits<FiniteElementType::dim,FiniteElementType::dim-1>::nVertices> vertices=eIter->second.simplex.child(boundaryFaces[f]).vertices();
//                        for(int v=0;v<vertices.size();++v)
//                        {
//                            isExternalBoundaryFace *= atx(vertices[v]->P0);
//                        }
//                        
//                        if(isExternalBoundaryFace)
//                        {
//                            temp.emplace_back(&(eIter->second),boundaryFaces[f]);
//                        }
//                    }
//                }
//            }
//            
//            return temp;
//        }
//        
    };
    
    
}	// close namespace
#endif

