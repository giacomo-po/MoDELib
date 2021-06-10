/* This file is part of PIL, the Particle Interaction Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LatticeTransitionMatrix_h_
#define model_LatticeTransitionMatrix_h_

#include <cfloat> // FLT_EPSILON
#include <Eigen/Dense>
//#include <RoundEigen.h>
#include <SmithDecomposition.h>
#include <Lattice.h>
#include <RationalMatrix.h>


namespace model
{
    /*!Class template that computes the Coincident Site Lattice (CSL) of two 
     * parent lattices using the Smith decomposition method [1].
     *
     * [1] Coincidence Lattices and Associated Shear Transformations
     */
    template <int dim>
    struct LatticeTransitionMatrix : public RationalMatrix<dim>
    {
        static_assert(dim>0,"dim must be > 0.");
        //        static constexpr double roundTol=FLT_EPSILON;
        typedef Eigen::Matrix<  double,dim,1> VectorDimD;
        typedef Eigen::Matrix<long int,dim,1> VectorDimI;
        typedef Eigen::Matrix<double,dim,dim> MatrixDimD;
        typedef Lattice<dim> LatticeType;
        typedef long long int IntValueType;
        typedef Eigen::Matrix<IntValueType,dim,dim> MatrixInt;

        
        
        /**********************************************************************/
        static IntValueType gcd(const IntValueType& a,const IntValueType& b)
        {
            return b>0? gcd(b, a % b) : a;
        }

        
        /**********************************************************************/
        static RationalMatrix<dim> getTransitionMatrix(const LatticeType& A,
                                                       const LatticeType& B)
        {
            const MatrixDimD R(B.latticeBasis*A.reciprocalBasis.transpose());
            
            // Check that R is a proper rotation
            const MatrixDimD RRT=R*R.transpose();
            const double RRTmInorm=(RRT-Eigen::Matrix<double,dim,dim>::Identity()).norm()/Eigen::Matrix<double,dim,dim>::Identity().norm();
            if(RRTmInorm>FLT_EPSILON)
            {
                std::cout<<"R="<<std::endl<<R<<std::endl;
                std::cout<<"R*R^T="<<std::endl<<RRT<<std::endl;
                std::cout<<"norm(R*R^T-I)/norm(I)="<<RRTmInorm<<", tol="<<FLT_EPSILON<<std::endl;
                assert(0 && "R IS NOT ORTHOGONAL.");
            }
            // make sure that C2G is proper
            assert(std::fabs(R.determinant()-1.0) < FLT_EPSILON && "R IS NOT PROPER.");
            
            
            // Compute the transition matrix T=inv(A)*B
            const MatrixDimD T=A.reciprocalBasis.transpose()*B.latticeBasis;
            
            // For the two lattices to have coincident sites, R must be a rational matrix
            // Compute the integer matrix P and the integer sigma such that T=P/sigma
            return RationalMatrix<dim>(T);
        }
        
        
        /**********************************************************************/
        LatticeTransitionMatrix(const LatticeType& A_in,
            const LatticeType& B_in) :
        /* init */ RationalMatrix<dim>(getTransitionMatrix(A_in,B_in))
        {
            
        }

    };
    
    
} // end namespace
#endif


