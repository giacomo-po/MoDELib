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
#include <cfloat>


namespace model
{
    
    
    template <int pOrder>
    class MatrixCompanion
    {/*!Class template that computes the roots of a polynomial of degree pOrder
      * as the eigenvalues of the matrix-companion of the polynomial.
      */
        enum{Ncoeff=pOrder+1};
        typedef Eigen::Matrix<double,1,Ncoeff> VectorNcoeff;
        typedef Eigen::Matrix<double,pOrder,pOrder> MatrixPorder;
        typedef Eigen::EigenSolver<Eigen::MatrixXd> EigenSolverType;
        
    public:
        
        const VectorNcoeff coeffs;
        const Eigen::MatrixXd cM; // the matrix companion
        const EigenSolverType es;
        const size_t rootSize;
        
        /**********************************************************************/
        MatrixCompanion(const VectorNcoeff& coeffsIN) :
        /* init list */ coeffs(coeffsIN),
        /* init list */ cM(getCM(coeffsIN)),
        /* init list */ es(cM),
        /* init list */ rootSize(cM.rows())
        {/*!@param[in] coeffsIN the coefficients of the polynomial, from highest 
          * to lowest.
          *
          * The constructor initializes the companion matrix cM, and the 
          * eigensolver es.
          */
        }
        
        /**********************************************************************/
        typename EigenSolverType::ComplexScalar root(const int& i) const
        {/*!@param[in] i an integer in 0<i<=pOrder
          *\returns the i-th root of the polynomial.
          */
            return es.eigenvalues()(i);
        }
        
        
        /**********************************************************************/
        static Eigen::MatrixXd getCM(const VectorNcoeff& coeffsIN)
        {
            
            //assert(std::fabs(coeffsIN(0))>FLT_EPSILON && "COEFFICIENT OF HIGHEST DEGREE CANNOT BE 0");
            
            Eigen::MatrixXd temp;
            
            if ( std::fabs(coeffsIN(0))>FLT_EPSILON)
            {
                temp=MatrixPorder::Zero();
                temp.row(0)=coeffsIN.template segment<pOrder>(1)/(-coeffsIN(0));
                temp.template block<pOrder-1,pOrder-1>(1,0).setIdentity();
            }
            else
            {
                temp=MatrixCompanion<pOrder-1>::getCM(coeffsIN.template segment<pOrder>(1));
            }
            
            return temp;
        }
        
    };
    
    /**********************************************************************/
    /**********************************************************************/
    template <>
    class MatrixCompanion<1>
    {
        enum{pOrder=1};
        enum{Ncoeff=pOrder+1};
        typedef Eigen::Matrix<double,1,Ncoeff> VectorNcoeff;
        typedef Eigen::Matrix<double,pOrder,pOrder> MatrixPorder;
        typedef Eigen::EigenSolver<Eigen::MatrixXd> EigenSolverType;
        
        
        
    public:
        
        const VectorNcoeff coeffs;
        const Eigen::MatrixXd cM; // the matrix companion
        const EigenSolverType es;
        const size_t rootSize;

        
        MatrixCompanion(const VectorNcoeff& coeffsIN) :
        /* init list */ coeffs(coeffsIN),
        /* init list */ cM(getCM(coeffsIN)),
        /* init list */ es(cM),
        /* init list */ rootSize(cM.rows())
        {
            
        }
        
        /**********************************************************************/
        typename EigenSolverType::ComplexScalar root(const int& i) const
        {
            return es.eigenvalues()(i);
        }
        
        /**********************************************************************/
        static MatrixPorder getCM(const VectorNcoeff& coeffsIN)
        {
            assert(std::fabs(coeffsIN(0))>FLT_EPSILON && "COEFFICIENT OF HIGHEST DEGREE CANNOT BE 0");
            MatrixPorder temp(MatrixPorder::Zero());
            temp.row(0)=coeffsIN.segment<pOrder>(1)/(-coeffsIN(0));
            return temp;
        }
        
    };
    
    
} // close namespace model
#endif
