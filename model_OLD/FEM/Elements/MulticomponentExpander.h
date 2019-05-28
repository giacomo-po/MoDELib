/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_MulticomponentExpander_H_
#define model_MulticomponentExpander_H_

#include <Eigen/Dense>

namespace model
{
    
	
    template<int c>
    struct MulticomponentExpander
    {
        
        template<typename T, int N>
        static Eigen::Matrix<T,c,c*N> expandSF(const Eigen::Matrix<T,1,N>& m)
        {
            Eigen::Matrix<T,c,c*N> temp;
            for (int j=0;j<N;++j)
            {
                temp.template block<c,c>(0,c*j).setIdentity()*=m(j);
                //temp.template block<c,c>(0,c*j)*=m(j);
                
            }
            return temp;
        }
        
        template<typename T,int dim, int N>
        static Eigen::Matrix<T,c*dim,c*N> expandSFgrad(const Eigen::Matrix<T,dim,N>& m)
        {
            Eigen::Matrix<T,c*dim,c*N> temp(Eigen::Matrix<T,c*dim,c*N>::Zero());
            for (int j=0;j<N;++j)
            {
                for (int i=0;i<c;++i)
                {
                    temp.template block<dim,1>(dim*i,c*j+i).setIdentity()=m.col(j);
                }
            }
            return temp;
        }

//        template<typename T, int N>
//        static Eigen::Matrix<T,(c*(c+1))/2,c*N> expandSFdef(const Eigen::Matrix<T,c,N>& m)
//        {
//            Eigen::Matrix<T,(c*(c+1))/2,c*N> temp(Eigen::Matrix<T,(c*(c+1))/2,c*N>::Zero());
//            
//            for (int n=0;n<N;++n) // loop over shape functions
//            {
//                int row=0;
//                for (int k=0;k<c;++k) // k is the index of the diagonal
//                {
//                    for (int j=0;j<c-k;++j)
//                    {
//                        temp(row,n*c+j+k)+=0.5*m(j  ,n);
//                        temp(row,n*c+j  )+=0.5*m(j+k,n);
//                        row++;
//                    }
//                }
//            }
//            return temp;
//        }
        
        template<typename T, int N>
        static Eigen::Matrix<T,(c*(c+1))/2,c*N> expandSFdef(const Eigen::Matrix<T,c,N>& m)
        {/*!@\param[in] m matrix of shape function gradients
          *\returns the matrix of shape function symmetric gradients 
          * (in engineering sense, that is e_{ij}=u_{i,j} if i==j, but 
          * e_{ij}=u_{i,j}+u_{j,i} if i!=j)
          */
            Eigen::Matrix<T,(c*(c+1))/2,c*N> temp(Eigen::Matrix<T,(c*(c+1))/2,c*N>::Zero());
            
            for (int n=0;n<N;++n) // loop over shape functions
            {
                int row=0;
                for (int k=0;k<c;++k) // k is the index of the diagonal
                {
                    for (int j=0;j<c-k;++j)
                    {
                        if(k==0)
                        {
                            temp(row,n*c+j)=m(j,n);
                        }
                        else
                        {
                            temp(row,n*c+j+k)=m(j  ,n);
                            temp(row,n*c+j  )=m(j+k,n);
                        }
                        row++;
                    }
                }
            }
            return temp;
        }
        
    };
    
    template<>
    struct MulticomponentExpander<1>
    {
        
        template<typename T, int N>
        static const Eigen::Matrix<T,1,N>& expandSF(const Eigen::Matrix<T,1,N>& m)
        {
            return m;
        }
        
        template<typename T,int dim, int N>
        static const Eigen::Matrix<T,dim,N>& expandSFgrad(const Eigen::Matrix<T,dim,N>& m)
        {
            return m;
        }
        
        template<typename T, int N>
        static const Eigen::Matrix<T,1,N>& expandSFdef(const Eigen::Matrix<T,1,N>& m)
        {
            return m;
        }
        
        
    };
    
    
}	// close namespace
#endif
