/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_BoundaryAtXmax_H_
#define model_BoundaryAtXmax_H_

#include <Eigen/Dense>
#include <model/FEM/Boundary/AtXmax>

namespace model
{
    template <int i>
    struct BoundaryAtXmax
    {
        
        
        /**********************************************************************/
        template <typename FiniteElementType, int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule>
        static IntegrationDomain<FiniteElementType::dim,1,qOrder,QuadratureRule> domain(const FiniteElementType& fe)
        {
            
            AtXmax<i> atXmax(fe);
            

            

            IntegrationDomain<FiniteElementType::dim,1,qOrder,QuadratureRule> temp;

            
            for (typename FiniteElementType::ElementContainerType::const_iterator eIter =fe.elementBegin();
                 /*                                                            */ eIter!=fe.elementEnd();
                 /*                                                            */ eIter++)
            {
                                    const std::vector<int> boundaryFaces=eIter->second.boundaryFaces();

                
                if(eIter->second.isBoundaryElement())
                {
                    const std::vector<int> boundaryFaces=eIter->second.boundaryFaces();
                    for (int f=0;f<boundaryFaces.size();++f)
                    {
                        
                        bool isExternalBoundaryFace(true);
                        
                        std::array<const Simplex<FiniteElementType::dim,0>*, SimplexTraits<FiniteElementType::dim,FiniteElementType::dim-1>::nVertices> vertices=eIter->second.simplex.child(boundaryFaces[f]).vertices();
                        for(int v=0;v<vertices.size();++v)
                        {
                            isExternalBoundaryFace *= atXmax(vertices[v]->P0);
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


//    struct BoundaryAtXmax
//    {
//
//        double top;
//
//        /**************************************/
//        template <typename FiniteElementType>
//        BoundaryAtXmax(const FiniteElementType& fe) :
//        top(-DBL_MAX)
//        {
//            for (int n=0;n<fe.nodeSize();++n)
//            {
//                if (fe.node(n).p0(fe.dim-1)>top)
//                {
//                    top=fe.node(n).p0(fe.dim-1);
//                }
//            }
//
//            std::cout<<"top of domain is "<<top<<std::endl;
//        }
//
//        template <int dim>
//        bool onExternalBoundary(const Eigen::Matrix<double,dim,1>& P0) const
//        {
//            return P0(dim-1)==top;
//        }
//
//
//    };

//                    for (int f=0;f<boundaryFaces.size();++f)
//                    {
//                        //temp.push_back(std::make_pair(n,boundaryFaces[f]));
//                        temp.emplace_back(n,boundaryFaces[f]);
//                    }
