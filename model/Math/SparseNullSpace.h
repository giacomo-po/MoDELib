/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>.
*
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _model_SparseNullSpace_h_
#define _model_SparseNullSpace_h_

#include <iostream>
#include <chrono>

#include <Eigen/Sparse>
#include <Eigen/SparseQR>
//#include <Eigen/SparseLU> // only for square matrics
#include <Eigen/OrderingMethods>
#include <model/MPI/MPIcout.h>

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

//            const auto t1= std::chrono::system_clock::now();
//            model::cout<<"Computing null-space... "<<std::flush;
//            Eigen::SparseLU<SparseMatrixType,Eigen::COLAMDOrdering<int> > lu(C.transpose()); 
//            model::cout<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t1)).count()<<" sec]"<<std::endl;

            
            const auto t0= std::chrono::system_clock::now();
            model::cout<<"Computing QR... "<<std::flush;


            
//            Eigen::SparseQR<SparseMatrixType,Eigen::COLAMDOrdering<int> > qr(C.transpose());
            Eigen::SparseQR<SparseMatrixType,Eigen::COLAMDOrdering<int> > qr(C.transpose());

            model::cout<<"done qr... "<<std::flush;

            
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
            
            model::cout<<"done loop... "<<std::flush;

            
            SparseMatrixType Q;
            Q=qr.matrixQ();
            
            model::cout<<"done Q... "<<std::flush;

            Z=Q.rightCols(Q.cols()-nnz);
            
            model::cout<<"done Z... "<<std::flush;

            Y=Q.leftCols(nnz);
            
            model::cout<<"done Y... "<<std::flush;

            
            model::cout<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<std::endl;


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
