/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2016 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SmithDecomposition_h_
#define model_SmithDecomposition_h_

#include <utility>      // std::pair, std::make_pair
#include <assert.h>
#include <Eigen/Core>

namespace model
{
    /*!Class template which computes the Smith normal (or canonical) form of a
     * square matrix of integers A. 
     *
     * The decomposition is:
     * A=X*D*Y
     * where
     * X and Y are unimodular integers matrices
     * D is a diagonal ingeter matrix such that D(k,k) divides D(k+1,k+1)
     * The class also computes the matrices U and V such that
     * D=U*A*V
     * where U=inv(X) and V=inv(Y) are also unimodular integer matrices.
     */
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
        static std::pair<int,int> findNonDivisible(const MatrixNi& A)
        {
            for(int i=1;i<N;++i)
            {
                for(int j=1;j<N;++j)
                {
                    if(A(i,j)%A(0,0)) // D(i,j) is no divisible by D(0,0)
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
            X=X*T;
            D=T*D;
        }
        
        /**********************************************************************/
        void  col_signChange(const int& ID)
        {/* unimodular matrix changing the sign of row/col ID
          */
            MatrixNi T(MatrixNi::Identity());
            T(ID,ID)=-1;
            V=V*T;
            Y=T*Y;
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
            X=X*T;
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
            Y=T*Y;
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
            T(i,j)=-n;
            X=X*T;
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
            T(j,i)=-n;
            Y=T*Y;
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

        MatrixNi D;
        MatrixNi U;
        MatrixNi V;
        MatrixNi X;
        MatrixNi Y;

        
    public:
        
        //        /**********************************************************************/
        //        SmithDecomposition(const Eigen::Matrix<int,N,N>& A) :
        //        /* delegate constructor */ SmithDecomposition(A.template cast<IntValueType>())
        //        {
        //        }
        
        /**********************************************************************/
        SmithDecomposition(const MatrixNi& A) :
        /* init */ D(A),
        /* init */ U(MatrixNi::Identity()),
        /* init */ V(MatrixNi::Identity()),
        /* init */ X(MatrixNi::Identity()),
        /* init */ Y(MatrixNi::Identity())
        {/*! The algorithm follows section 2 in GEORGE HAVAS AND BOHDAN S. MAJEWSKI,
          * Integer Matrix Diagonalization,
          * J. Symbolic Computation (1997) 24, 399–408.
          */
            
            std::pair<int,int> rc=std::make_pair(0,0);
            while (rc.first>=0)
            {
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
                        const IntValueType q=D(i,0)/D(0,0);
                        row_add_multiply(i,0,-q);
                    }
                    
                    // Subtract suitable multiples of the first column from the other
                    // columns to replace each entry in the first row, except a(0,0), by zero
                    for (int j=1;j<N;++j)
                    {
                        const IntValueType q=D(0,j)/D(0,0);
                        col_add_multiply(j,0,-q);
                    }
                    
                    // Now D is in the form [D(0,0) 0s;0s D'] where D' is a (N-1)x(N-1) matrix
                    rc=findNonDivisible(D);
                    if(rc.first>0)
                    {
                        row_add_multiply(0,rc.first,1);
                    }
                    
                    
                    
                    
                }
                else
                {
                    rc=std::make_pair(-1,-1); // end main while loop immediately
                }
                

            }

            // Recursive iteration
            const SmithDecomposition<N-1> sd(D.template block<N-1,N-1>(1,1));
            D.template block<N-1,N-1>(1,1)=sd.matrixD();
            U.template block<N-1,N>(1,0)=sd.matrixU()*U.template block<N-1,N>(1,0);
            V.template block<1,N-1>(0,1)=V.template block<1,N-1>(0,1)*sd.matrixV();
            V.template block<N-1,N-1>(1,1)=V.template block<N-1,N-1>(1,1)*sd.matrixV();

            Y.template block<N-1,N>(1,0)=sd.matrixY()*Y.template block<N-1,N>(1,0);
            X.template block<1,N-1>(0,1)=X.template block<1,N-1>(0,1)*sd.matrixX();
            X.template block<N-1,N-1>(1,1)=X.template block<N-1,N-1>(1,1)*sd.matrixX();

            
            // Find a non-zero value of A and make it D(0,0)
            
            
            assert((U*A*V-D).squaredNorm()==0);
            assert((X*D*Y-A).squaredNorm()==0);

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
        
        /**********************************************************************/
        const MatrixNi& matrixX() const
        {
            return X;
        }
        
        /**********************************************************************/
        const MatrixNi& matrixY() const
        {
            return Y;
        }
    };
    
    /**************************************************************************/
    /**************************************************************************/
    template <>
    class SmithDecomposition<1>
    {
        typedef long long int IntValueType;
        typedef Eigen::Matrix<IntValueType,1,1> MatrixNi;

        const MatrixNi D;
        const MatrixNi U;
        const MatrixNi V;
        const MatrixNi X;
        const MatrixNi Y;
        
    public:
        
        /**********************************************************************/
        SmithDecomposition(const MatrixNi& A) :
        /* init */ D(A),
        /* init */ U(MatrixNi::Identity()),
        /* init */ V(MatrixNi::Identity()),
        /* init */ X(MatrixNi::Identity()),
        /* init */ Y(MatrixNi::Identity())

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
        
        /**********************************************************************/
        const MatrixNi& matrixX() const
        {
            return X;
        }
        
        /**********************************************************************/
        const MatrixNi& matrixY() const
        {
            return Y;
        }
    };
    
}
#endif
