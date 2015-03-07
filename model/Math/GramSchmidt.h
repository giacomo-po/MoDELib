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
#include <Eigen/StdVector>

namespace model
{
    
    template <short unsigned int dim>
    struct GramSchmidt : public std::vector<Eigen::Matrix<double,dim,1>,Eigen::aligned_allocator<Eigen::Matrix<double,dim,1> > >
    {
        /*! \brief A class template that performs a modified Gram-Schmidt (MGS)
         * orthonormalization (with reorthogonalization) of a given set of vectors
         * in dimension dim.
         */
        
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef std::vector<VectorDim,Eigen::aligned_allocator<VectorDim> > VectorOfNormalsType;
        
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        
        /********************************************************************/
        GramSchmidt(const VectorOfNormalsType& NV, const size_t ID, const double& tol=Eigen::NumTraits<double>::dummy_precision())
        {/*
          * Ref:
          * Giraud, L., Langou, J., & Rozloznik, M. (2005). The loss of 
          * orthogonality in the Gram-Schmidt orthogonalization process.
          * Computers & Mathematics with Applications, 50(7), 1069–1075.
          */
            for (size_t i=0;i<NV.size();++i)
            {
                VectorDim temp=NV[i];
                
                for (int r=0;r<2;++r) // 2 is the number or reorthogonalizations
                {
                    for (size_t j=0;j<this->size();++j)
                    {
                        //temp-= NV[i].dot(this->operator[](j))*this->operator[](j); // classical GS (CGM)
                        temp -= (temp.dot(this->operator[](j))*this->operator[](j)).eval(); //modified GS (MGS)
                        
                    }
                }
                
                const double tempNorm(temp.norm());
                if (tempNorm>tol)
                {
                    this->push_back(temp/tempNorm);
                }
            }
            
            if(this->size()>dim)
            {
                std::cout<<"GramSchmidt FAILED. Input vectors:"<<std::endl;
                for(const auto& v : NV)
                {
                    std::cout<<std::setprecision(15)<<std::scientific<<v.transpose()<<std::endl;
                }
                std::cout<<"orthonormal vectors:"<<std::endl;
                for(const auto& v : *this)
                {
                    std::cout<<std::setprecision(15)<<std::scientific<<v.transpose()<<std::endl;
                }
                std::cout<<"tol="<<tol<<std::endl;
                std::cout<<"Node="<<ID<<std::endl;

                assert(0 && "GramSchmidt failed.");
            }
        }
        
    };
    
}
#endif
