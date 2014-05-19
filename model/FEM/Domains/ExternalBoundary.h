/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_ExternalBoundary_H_
#define model_ExternalBoundary_H_

#include <deque>
#include <utility>      // std::pair, std::make_pair
#include <Eigen/Dense>
#include <model/Mesh/Simplex.h>
#include <model/FEM/Domains/IntegrationDomain.h>

namespace model
{
    
 
    /**************************************************************************/
	/**************************************************************************/
//    template <typename BndType, int dim, int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule>
	class ExternalBoundary //: public IntegrationDomain<dim,1,qOrder,QuadratureRule>
    {

//        const BndType bnd;
        
    public:
//        /**********************************************************************/
//        template <typename FiniteElementType>
//        ExternalBoundary(const FiniteElementType& fe) : bnd(fe)
//        {
//            //const T bnd(fe);
//            
//            for (int n=0;n<fe.elementSize();++n)
//            {
//                if(fe.element(n).isBoundaryElement())
//                {
//                    const std::vector<int> boundaryFaces=fe.element(n).boundaryFaces();
//                    
//
//                    for (int f=0;f<boundaryFaces.size();++f)
//                    {
//
//                        bool isExternalBoundaryFace(true);
//                        
//                        std::array<const Simplex<FiniteElementType::dim,0>*, SimplexTraits<FiniteElementType::dim,FiniteElementType::dim-1>::nVertices> vertices=fe.element(n).simplex.child(boundaryFaces[f]).vertices();
//                        for(int v=0;v<vertices.size();++v)
//                        {
//                            isExternalBoundaryFace *= bnd.onExternalBoundary(vertices[v]->P0);
//                        }
//                        
//                        if(isExternalBoundaryFace)
//                        {
//                            this->push_back(std::make_pair(n,boundaryFaces[f]));
//                        }
//                    }
//                }
//            }
//        }
        
        
        
        template <typename FiniteElementType, int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule>
        static IntegrationDomain<FiniteElementType::dim,1,qOrder,QuadratureRule> domain(const FiniteElementType& fe)
        {
            IntegrationDomain<FiniteElementType::dim,1,qOrder,QuadratureRule> temp;
            for (int n=0;n<fe.elementSize();++n)
            {
                if(fe.element(n).isBoundaryElement())
                {
                    const std::vector<int> boundaryFaces=fe.element(n).boundaryFaces();
                    for (int f=0;f<boundaryFaces.size();++f)
                    {
                        //temp.push_back(std::make_pair(n,boundaryFaces[f]));
                        temp.emplace_back(n,boundaryFaces[f]);
                    }
                }
            }
            return temp;
        }
        
    };

}	// close namespace
#endif