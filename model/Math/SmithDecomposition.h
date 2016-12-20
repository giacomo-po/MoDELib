/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SmithDecomposition_h_
#define model_SmithDecomposition_h_

#include <utility>      // std::pair, std::make_pair
//#include <tuple>
#include <assert.h>
#include <Eigen/Core>



namespace model
{
    template <int N>
    class SmithDecomposition
    {
        typedef long long int IntValueType;
        typedef Eigen::Matrix<IntValueType,N,N> MatrixNi;
        
        /**********************************************************************/
        static IntValueType gcd(const IntValueType& a,const IntValueType& b)
        {
            return b>0? gcd(b, a % b) : a;
        }
        
        /**********************************************************************/
        template <typename Derived>
        static IntValueType gcd(const Eigen::MatrixBase<Derived>& v)
        {
            IntValueType g=v(0);
            for(int n=1;n<v.size();++n)
            {
                g=gcd(g,v(n));
            }
            return g;
        }
        
        /**********************************************************************/
        static std::pair<int,int> findNonZero(const MatrixNi& A)
        {
            for(int i=0;i<N;++i)
            {
                for(int j=0;j<N;++j)
                {
                    if(A(i,j)!=0)
                    {
                        return std::make_pair(i,j);
                    }
                }
            }
            return std::make_pair(-1,-1);
        }
        
        /**********************************************************************/
        void  row_signChange(const int& ID)
        {/* unimodular matrix changing the sign of row/col ID
          */
            MatrixNi T(MatrixNi::Identity());
            T(ID,ID)=-1;
            U=T*U;
            D=T*D;
        }
        
        /**********************************************************************/
        void  col_signChange(const int& ID)
        {/* unimodular matrix changing the sign of row/col ID
          */
            MatrixNi T(MatrixNi::Identity());
            T(ID,ID)=-1;
            V=V*T;
            D=D*T;
        }
        
        /**********************************************************************/
        void row_swap(const int& i, const int& j)
        {/* unimodular matrix swapping rows/cols i and j
          */
            MatrixNi T(MatrixNi::Identity());
            T(i,i)=0;
            T(j,j)=0;
            T(i,j)=1;
            T(j,i)=1;
            U=T*U;
            D=T*D;
        }
        
        /**********************************************************************/
        void col_swap(const int& i, const int& j)
        {/* unimodular matrix swapping rows/cols i and j
          */
            MatrixNi T(MatrixNi::Identity());
            T(i,i)=0;
            T(j,j)=0;
            T(i,j)=1;
            T(j,i)=1;
            V=V*T;
            D=D*T;
        }
        
        /**********************************************************************/
        void row_add_multiply(const int& i, const int& j,const IntValueType& n)
        {/* unimodular matrix adding n times row (col) j to row (col) i
          */
            // note that add applied to colums needs transpose!
            MatrixNi T(MatrixNi::Identity());
            T(i,j)=n;
            U=T*U;
            D=T*D;
        }
        
        /**********************************************************************/
        void col_add_multiply(const int& i, const int& j,const IntValueType& n)
        {/* unimodular matrix adding n times row (col) j to row (col) i
          */
            // note that add applied to colums needs transpose!
            MatrixNi T(MatrixNi::Identity());
            T(j,i)=n;
            V=V*T;
            D=D*T;
        }
        
        /**********************************************************************/
        bool reduceFirstRow()
        {
            // Operate on first row so that D(0,j) is divisible by D(0,0)
            bool didOperate=false;
            int j=1;
            while (j<N)
            {
                const IntValueType r=D(0,j)%D(0,0);
                const IntValueType q=D(0,j)/D(0,0);
                if(r)
                {
                    col_add_multiply(j,0,-q);
                    col_swap(0,j);
                    j=1;
                    didOperate=true;
                }
                else
                {
                    j++;
                }
            }
            return didOperate;
        }
        
        /**********************************************************************/
        bool reduceFirstCol()
        {
            // Operate on first row so that D(0,j) is divisible by D(0,0)
            bool didOperate=false;
            int j=1;
            while (j<N)
            {
                const IntValueType r=D(j,0)%D(0,0);
                const IntValueType q=D(j,0)/D(0,0);
                if(r)
                {
                    row_add_multiply(j,0,-q);
                    row_swap(0,j);
                    j=1;
                    didOperate=true;
                }
                else
                {
                    j++;
                }
            }
            return didOperate;
        }
        
        MatrixNi U;
        MatrixNi D;
        MatrixNi V;
        
    public:

//        /**********************************************************************/
//        SmithDecomposition(const Eigen::Matrix<int,N,N>& A) :
//        /* delegate constructor */ SmithDecomposition(A.template cast<IntValueType>())
//        {
//        }
        
        /**********************************************************************/
        SmithDecomposition(const MatrixNi& A) :
        /* init */ U(MatrixNi::Identity()),
        /* init */ D(A),
        /* init */ V(MatrixNi::Identity())
        {/*! The algorithm follows section 2 in GEORGE HAVAS AND BOHDAN S. MAJEWSKI,
          * Integer Matrix Diagonalization,
          * J. Symbolic Computation (1997) 24, 399–408.
          */
            
            // Find a non-zero value of A and make it D(0,0)
            std::pair<int,int> ij=findNonZero(D);
            if(ij.first>=0)
            {
                // Swap row(i) and col (j) to amke D(0,0) non-zero.
                row_swap(0,ij.first);
                col_swap(0,ij.second);
                
                // Operate on first row and col so that all entried are divisible by D(0,0)
                bool didOperate=true;
                while (didOperate)
                {
                    didOperate=reduceFirstRow()+reduceFirstCol();
                }
                
                // If D(0,0) is negative make it positive
                if(D(0,0)<0)
                {
                    row_signChange(0);
                }
                
                // Subtract suitable multiples of the first row from the other
                // rows to replace each entry in the first column, except a(0,0), by zero
                for (int i=1;i<N;++i)
                {
                    const int q=D(i,0)/D(0,0);
                    row_add_multiply(i,0,-q);
                }
                
                // Subtract suitable multiples of the first column from the other
                // columns to replace each entry in the first row, except a(0,0), by zero
                for (int j=1;j<N;++j)
                {
                    const int q=D(0,j)/D(0,0);
                    col_add_multiply(j,0,-q);
                }
                
                
                // Recursive iteration
                const SmithDecomposition<N-1> sd(D.template block<N-1,N-1>(1,1));
                
//                std::cout<<"U old="<<std::endl<<U<<std::endl;
//                std::cout<<"sd.U="<<std::endl<<sd.matrixU()<<std::endl;

//                std::cout<<"V old="<<std::endl<<V<<std::endl;
//                std::cout<<"sd.V="<<std::endl<<sd.matrixV()<<std::endl;

                
                D.template block<N-1,N-1>(1,1)=sd.matrixD();
                U.template block<N-1,N>(1,0)=sd.matrixU()*U.template block<N-1,N>(1,0);
                V.template block<1,N-1>(0,1)=V.template block<1,N-1>(0,1)*sd.matrixV();
                V.template block<N-1,N-1>(1,1)=V.template block<N-1,N-1>(1,1)*sd.matrixV();
                
//                std::cout<<"U new="<<std::endl<<U<<std::endl;
//                std::cout<<"V new="<<std::endl<<V<<std::endl;

                
            }
            
            
//            assert((U*A*V-D).squaredNorm()==0);
        }
        
        /**********************************************************************/
        const MatrixNi& matrixU() const
        {
            return U;
        }
        
        /**********************************************************************/
        const MatrixNi& matrixD() const
        {
            return D;
        }
        
        /**********************************************************************/
        const MatrixNi& matrixV() const
        {
            return V;
        }
        
    };
    
    /**************************************************************************/
    /**************************************************************************/
    template <>
    class SmithDecomposition<1>
    {
        typedef long long int IntValueType;
        typedef Eigen::Matrix<IntValueType,1,1> MatrixNi;
        
        const MatrixNi U;
        const MatrixNi D;
        const MatrixNi V;
        
    public:
        
        /**********************************************************************/
        SmithDecomposition(const MatrixNi& A) :
        /* init */ U(MatrixNi::Identity()),
        /* init */ D(A),
        /* init */ V(MatrixNi::Identity())
        {
        }
        /**********************************************************************/
        const MatrixNi& matrixU() const
        {
            return U;
        }
        
        /**********************************************************************/
        const MatrixNi& matrixD() const
        {
            return D;
        }
        
        /**********************************************************************/
        const MatrixNi& matrixV() const
        {
            return V;
        }
    };
    
}
#endif



//        /**********************************************************************/
//        void iterate()
//        {
//
//
//
//
//
//            int rowStart=0;
//            row_iteration(rowStart);
//            std::cout<<D<<std::endl<<std::endl;
//
//            col_iteration(rowStart);
//            std::cout<<D<<std::endl<<std::endl;
//
//        }
//
//        /**********************************************************************/
//        void row_iteration(const int& rowStart)
//        {
//            const int g=gcd(D.col(0));
//            std::cout<<"g_row="<<g<<std::endl;
//
//            int pivotI=-1;
//            int pivotJ=-1;
//            int pivotN=-1;
//            for(int i=0;i<N;++i)
//            {
//                for(int j=0;j<N;++j)
//                {
//                    if(j!=i)
//                    {
//                        if(((g-D(i,0))%D(j,0))==0)
//                        {
//                            pivotI=i;
//                            pivotJ=j;
//                            pivotN=(g-D(i,0))/D(j,0);
//                        }
//                        else if(((g+D(i,0))%D(j,0))==0)
//                        {
//                            assert(0);
//                            row_signChange(i);
//                            iterate();
//                        }
//                    }
//                    if(pivotI>-1)
//                    {
//                        break;
//                    }
//                }
//                if(pivotI>-1)
//                {
//                    break;
//                }
//            }
//            assert(pivotI>-1);
//
//            row_add_multiply(pivotI,pivotJ,pivotN);
//
//
//            for(int i=0;i<N;++i)
//            {
//                if(i!=pivotI)
//                {
//                    int m=D(i,0)/g;
//                    row_add_multiply(i,pivotI,-m);
//                }
//            }
//
//        }

//        /**********************************************************************/
//        void col_iteration(const int& rowStart)
//        {
//            const int g=gcd(D.row(0));
//            std::cout<<"g_col="<<g<<std::endl;
//            int pivotI=-1;
//            int pivotJ=-1;
//            int pivotN=-1;
//            for(int i=0;i<N;++i)
//            {
//                for(int j=0;j<N;++j)
//                {
//                    if(j!=i)
//                    {
//                        if(((g-D(0,i))%D(0,j))==0)
//                        {
//                            pivotI=i;
//                            pivotJ=j;
//                            pivotN=(g-D(0,i))/D(0,j);
//                        }
//                        else if(((g+D(0,i))%D(0,j))==0)
//                        {
//                            assert(0);
//                            col_signChange(i);
//                            iterate();
//                        }
//                    }
//                    if(pivotI>-1)
//                    {
//                        break;
//                    }
//
//                }
//                if(pivotI>-1)
//                {
//                    break;
//                }
//            }
//            assert(pivotI>-1);
//
//            col_add_multiply(pivotI,pivotJ,pivotN);
//
//
//            for(int i=0;i<N;++i)
//            {
//                if(i!=pivotI)
//                {
//                    int m=D(0,i)/g;
//                    col_add_multiply(i,pivotI,-m);
//                }
//            }
//
//        }
