/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_MatrixCompanion_h_
#define model_MatrixCompanion_h_

#include <Eigen/Dense>
#include <assert.h>


namespace model{
    
    template <int pOrder>
    class MatrixCompanion
    {
        enum{Ncoeff=pOrder+1};
        typedef Eigen::Matrix<double,1,Ncoeff> VectorNcoeff;
        typedef Eigen::Matrix<double,pOrder,pOrder> MatrixPorder;
        typedef Eigen::EigenSolver<MatrixPorder> EigenSolverType;
        
        
        

        MatrixPorder getCM(const VectorNcoeff& coeffsIN) const
        {
            
            assert(std::fabs(coeffsIN(0))>FLT_EPSILON && "COEFFICIENT OF HIGHEST DEGREE CANNO BE 0");
            
            MatrixPorder temp(MatrixPorder::Zero());
            temp.row(0)=coeffsIN.template segment<pOrder>(1)/(-coeffsIN(0));
            temp.template block<pOrder-1,pOrder-1>(1,0).setIdentity();
            
            
//            temp<< coeffsIN.template segment<pOrder>(1)/(-coeffsIN(0)),
//                   Eigen::Matrix<double,pOrder-1,pOrder-1>::Identity(),
//            Eigen::Matrix<double,pOrder-1,1>::Zero();
            return temp;
        }
        
        
    public:

        const VectorNcoeff coeffs;
        const MatrixPorder cM; // the matrix companion
        const EigenSolverType es;
        
        MatrixCompanion(const VectorNcoeff& coeffsIN) :
        /* init list */ coeffs(coeffsIN),
        /* init list */ cM(getCM(coeffsIN)),
        /* init list */ es(cM)
        {
        
        }
        
        
        typename EigenSolverType::ComplexScalar root(const int& i) const
        {
            return es.eigenvalues()(i);
        }
        
        
    };
    
    
    
} // close namespace model
#endif