/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_AtXmax_H_
#define model_AtXmax_H_

#include <vector>
#include <Eigen/Dense>
#include <model/FEM/Domains/IntegrationDomain.h>

namespace model
{
    template <int i>
    class   AtXmax
    {
        
        const double max;
        
    public:
        
        /**************************************/
        template <typename FiniteElementType>
        AtXmax(const FiniteElementType& fe) :
        max(fe.xMax()(i))
        {
            
        }
        
        /**************************************/
        template <typename NodeType>
        bool operator()(const NodeType& node) const
        {
            return std::abs(node.p0(i)-max)<FLT_EPSILON;
        }
        
        /**************************************/
        template <int dim>
        bool operator()(const Eigen::Matrix<double,dim,1>& p0) const
        {
            return std::abs(p0(i)-max)<FLT_EPSILON;
        }
        
        /**********************************************************************/
        template <typename FiniteElementType, int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule>
        static IntegrationDomain<FiniteElementType,1,qOrder,QuadratureRule> boundary(const FiniteElementType& fe)
        {

            AtXmax<i> atx(fe);
            
            IntegrationDomain<FiniteElementType,1,qOrder,QuadratureRule> temp;
            
            
            for (typename FiniteElementType::ElementContainerType::const_iterator eIter =fe.elementBegin();
                 /*                                                            */ eIter!=fe.elementEnd();
                 /*                                                            */ eIter++)
            {
                if(eIter->second.isBoundaryElement())
                {
                    const std::vector<int> boundaryFaces=eIter->second.boundaryFaces();
                    for (int f=0;f<boundaryFaces.size();++f)
                    {
                        bool isExternalBoundaryFace(true);
                        
                        std::array<const Simplex<FiniteElementType::dim,0>*, SimplexTraits<FiniteElementType::dim,FiniteElementType::dim-1>::nVertices> vertices=eIter->second.simplex.child(boundaryFaces[f]).vertices();
                        for(int v=0;v<vertices.size();++v)
                        {
                            isExternalBoundaryFace *= atx(vertices[v]->P0);
                        }
                        
                        if(isExternalBoundaryFace)
                        {
                            temp.emplace_back(&(eIter->second),boundaryFaces[f]);
                        }
                    }
                }
            }
            
            return temp;
        }
        
        
    };
    
    
}	// close namespace
#endif

