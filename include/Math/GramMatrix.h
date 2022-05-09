/* This file is part of model.
 *
 * model is distributed without any warranty under the MIT License.
 */


#ifndef model_GramMatrix_h_
#define model_GramMatrix_h_

#include <Eigen/Core>
#include <array>

namespace model
{
    template<typename T,int N>
    class GramMatrix : public Eigen::Matrix<T,N,N>
    {
        
        template<int dim>
        static Eigen::Matrix<T,N,N> getMatrix(const std::array<Eigen::Matrix<T,dim,1>,N>& a)
        {
            Eigen::Matrix<T,dim,N> X(Eigen::Matrix<T,dim,N>::Zero());
            for(int k=0;k<N;++k)
            {
                X.col(k)=a[k];
            }
            return X.transpose()*X;
        }


        
    public:
        
        template<int dim>
        GramMatrix(const std::array<Eigen::Matrix<T,dim,1>,N>& a) :
        /* init */ Eigen::Matrix<T,N,N>(getMatrix(a))
        {

            
        }
        
        
    };
    

    
//    Rational rat(const long int& n,const long int& d)
//    {
//        return Rational(n,d);
//    }
    
} // end namespace



#endif
