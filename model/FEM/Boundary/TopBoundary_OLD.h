/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_TopBoundary_H_
#define model_TopBoundary_H_

#include <Eigen/Dense>

namespace model
{
    
    struct TopBoundary
    {
        

        
        template <int dim>
        static bool onTopBoundary(const Eigen::Matrix<double,dim,1>& P0, const double& top)
        {
            return P0(dim-1)==top;
        }
        
        
        template <typename FiniteElementType, int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule>
        static IntegrationDomain<FiniteElementType::dim,1,qOrder,QuadratureRule> domain(const FiniteElementType& fe)
        {
            
            const double t0(clock());
            std::cout<<"TopBoundary, making IntegrationDomain..."<<std::flush;
            double top(-DBL_MAX);
            
            
            for (int n=0;n<fe.nodeSize();++n)
            {
                if (fe.node(n).p0(fe.FiniteElementType::dim-1)>top)
                {
                    top=fe.node(n).p0(FiniteElementType::dim-1);
                }
            }
            
            
//            for (int n=0;n<fe.elementSize();++n)
//            {
////                if(fe.element(n).isBoundaryElement())
//                const int b(fe.element(n).isBoundaryElement());
//                //if(b){
//                std::vector<int> v(1,b);
//                //}
////                                    }
//            }

            for (typename FiniteElementType::ElementContainerType::const_iterator eIter =fe.elementBegin();
                 /*                                                            */ eIter!=fe.elementEnd();
                 /*                                                            */ eIter++)
            {
                //                if(fe.element(n).isBoundaryElement())
                const int b(eIter->second.isBoundaryElement());
                //if(b){
                std::vector<int> v(1,b);
                //}
                //                                    }
            }
            

            IntegrationDomain<FiniteElementType::dim,1,qOrder,QuadratureRule> temp;
//            for (int n=0;n<fe.elementSize();++n)
//            {
//                if(fe.element(n).isBoundaryElement())
//                {
//                    
//                    
////                    const std::vector<int> boundaryFaces=fe.element(n).boundaryFaces();
//
////
////                    for (int f=0;f<boundaryFaces.size();++f)
////                    {
////                        
////                        bool isExternalBoundaryFace(true);
////                        
////                        std::array<const Simplex<FiniteElementType::dim,0>*, SimplexTraits<FiniteElementType::dim,FiniteElementType::dim-1>::nVertices> vertices=fe.element(n).simplex.child(boundaryFaces[f]).vertices();
////                        for(int v=0;v<vertices.size();++v)
////                        {
////                            isExternalBoundaryFace *= onTopBoundary(vertices[v]->P0,top);
////                        }
////                        
////                        if(isExternalBoundaryFace)
////                        {
////                            temp.emplace_back(n,boundaryFaces[f]);
////                        }
////                    }
//                    
//                }
//            }
            
            std::cout<<" done.["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]"<<std::endl;
            
            
            return temp;
        }
        
    };
    

    
    
}	// close namespace
#endif


//    struct TopBoundary
//    {
//
//        double top;
//
//        /**************************************/
//        template <typename FiniteElementType>
//        TopBoundary(const FiniteElementType& fe) :
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
