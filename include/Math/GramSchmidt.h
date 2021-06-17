/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_GRAMSCHMIDT_H_
#define model_GRAMSCHMIDT_H_

#include <iostream>
#include <iomanip>
#include <vector>
#include <Eigen/Dense>

namespace model
{
    
    struct GramSchmidt
    {
        /*! \brief A class template that performs a modified Gram-Schmidt (MGS)
         * orthonormalization (with reorthogonalization) of a given set of vectors
         * in dimension dim.
         */
        
        
        
        /********************************************************************/
        template <int dim>
        static void orthoNormalize(std::vector<Eigen::Matrix<double,dim,1>>& NV,
                              const double& tol=Eigen::NumTraits<double>::dummy_precision())
        {/*
          * Ref:
          * Giraud, L., Langou, J., & Rozloznik, M. (2005). The loss of 
          * orthogonality in the Gram-Schmidt orthogonalization process.
          * Computers & Mathematics with Applications, 50(7), 1069–1075.
          */
            typedef Eigen::Matrix<double,dim,1> VectorDim;
            typedef std::vector<VectorDim> VectorOfNormalsType;
            VectorOfNormalsType NV_new;
            
            for (size_t i=0;i<NV.size();++i)
            {
                VectorDim temp=NV[i];
                
                for (int r=0;r<2;++r) // 2 is the number or reorthogonalizations
                {
                    for (size_t j=0;j<NV_new.size();++j)
                    {
                        //temp-= NV[i].dot(NV_new[j]))*NV_new[j]; // classical GS (CGM)
                        temp -= (temp.dot(NV_new[j])*NV_new[j]).eval(); //modified GS (MGS)
                        
                    }
                }
                
                const double tempNorm(temp.norm());
                if (tempNorm>tol)
                {
                    NV_new.push_back(temp/tempNorm);
                }
            }
            
            if(NV_new.size()>dim)
            {
                std::cout<<"GramSchmidt FAILED. Input vectors:"<<std::endl;
                for(const auto& v : NV)
                {
                    std::cout<<std::setprecision(15)<<std::scientific<<v.transpose()<<std::endl;
                }
                std::cout<<"orthonormal vectors:"<<std::endl;
                for(const auto& v : NV_new)
                {
                    std::cout<<std::setprecision(15)<<std::scientific<<v.transpose()<<std::endl;
                }
                std::cout<<"tol="<<tol<<std::endl;
                assert(0 && "GramSchmidt failed.");
            }
            
            NV=NV_new; // overwrite
        }
        
//        /**********************************************************************/
//        static void makeUnique(std::deque<const LatticePlane*>& NV)
//        {
//            
//            std::deque<const LatticePlane*> temp;
//            
//            for (size_t i=0;i<NV.size();++i)
//            {
//                bool unique=true;
//                for (size_t j=0;j<temp.size();++j)
//                {
//                    unique*=(temp[j]->n.cross(NV[i]->n).squaredNorm());
//                }
//                
//                if (unique)
//                {
//                    temp.push_back(NV[i]);
//                }
//                
//                if(temp.size()>=3)
//                {
//                    break;
//                }
//            }
//            
////            if(temp.size()>3)
////            {
//////                std::cout<<"GramSchmidt FAILED. Input vectors:"<<std::endl;
//////                for(const auto& v : NV)
//////                {
//////                    std::cout<<std::setprecision(15)<<std::scientific<<v.transpose()<<std::endl;
//////                }
//////                std::cout<<"orthonormal vectors:"<<std::endl;
//////                for(const auto& v : *this)
//////                {
//////                    std::cout<<std::setprecision(15)<<std::scientific<<v.transpose()<<std::endl;
//////                }
//////                std::cout<<"tol="<<tol<<std::endl;
//////                std::cout<<"Node="<<ID<<std::endl;
////                
////                assert(0 && "GramSchmidt failed.");
////            }
//            
//            NV=temp; // overwrite
//        }
        
    };
    
}
#endif
