/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>.
*
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _model_SparseNullSpace_h_
#define _model_SparseNullSpace_h_

#include <Eigen/Sparse>
#include <Eigen/SparseQR>
#include <Eigen/OrderingMethods>

// http://stackoverflow.com/questions/2181418/computing-the-null-space-of-a-matrix-as-fast-as-possible

namespace model
{

    template<typename SparseMatrixType>
    class SparseNullSpace
    {
        
        typedef typename SparseMatrixType::Scalar Scalar;
        
        SparseMatrixType Z;
        SparseMatrixType Y;
        
    public:
        
        SparseNullSpace(const SparseMatrixType& C, const Scalar& tol=Eigen::NumTraits<Scalar>::dummy_precision())
        {
            Eigen::SparseQR<SparseMatrixType,Eigen::COLAMDOrdering<int> > qr(C.transpose());

            
            int nnz=0; // number of non-zero diagonal elements
            for(int c=0;c<qr.matrixR().cols();++c)
            {
                if(std::fabs(qr.matrixR().coeff(c,c))<= tol)
                {
                     break;
                }
                else
                {
                    nnz++;
                }
            }
            
            SparseMatrixType Q;
            Q=qr.matrixQ();
            Z=Q.rightCols(Q.cols()-nnz);
            Y=Q.leftCols(nnz);

        }
        
        const SparseMatrixType& matrixZ() const
        {
            return Z;
        }

        const SparseMatrixType& matrixY() const
        {
            return Y;
        }

    };
    
}

#endif
