/* This file is part of PIL, the Particle Interaction Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LatticeTransitionMatrix_cpp_
#define model_LatticeTransitionMatrix_cpp_
#include <LatticeModule.h>
#include <LatticeTransitionMatrix.h> // FLT_EPSILON

namespace model
{
    /*!Class template that computes the Coincident Site Lattice (CSL) of two 
     * parent lattices using the Smith decomposition method [1].
     *
     * [1] Coincidence Lattices and Associated Shear Transformations
     */
    /**********************************************************************/
    template <int dim>
    typename LatticeTransitionMatrix<dim>::IntValueType LatticeTransitionMatrix<dim>::gcd(const IntValueType &a, const IntValueType &b)
    {
        return b > 0 ? gcd(b, a % b) : a;
    }

    /**********************************************************************/
    template <int dim>
    RationalMatrix<dim> LatticeTransitionMatrix<dim>::getTransitionMatrix(const LatticeType &A,
                                                   const LatticeType &B)
    {
        const MatrixDimD R(B.latticeBasis * A.reciprocalBasis.transpose());

        // Check that R is a proper rotation
        const MatrixDimD RRT = R * R.transpose();
        const double RRTmInorm = (RRT - Eigen::Matrix<double, dim, dim>::Identity()).norm() / Eigen::Matrix<double, dim, dim>::Identity().norm();
        if (RRTmInorm > FLT_EPSILON)
        {
            std::cout << "R=" << std::endl
                      << R << std::endl;
            std::cout << "R*R^T=" << std::endl
                      << RRT << std::endl;
            std::cout << "norm(R*R^T-I)/norm(I)=" << RRTmInorm << ", tol=" << FLT_EPSILON << std::endl;
            assert(0 && "R IS NOT ORTHOGONAL.");
        }
        // make sure that C2G is proper
        assert(std::fabs(R.determinant() - 1.0) < FLT_EPSILON && "R IS NOT PROPER.");

        // Compute the transition matrix T=inv(A)*B
        const MatrixDimD T = A.reciprocalBasis.transpose() * B.latticeBasis;

        // For the two lattices to have coincident sites, R must be a rational matrix
        // Compute the integer matrix P and the integer sigma such that T=P/sigma
        return RationalMatrix<dim>(T);
    }

    /**********************************************************************/
    template <int dim>
    LatticeTransitionMatrix<dim>::LatticeTransitionMatrix(const LatticeType &A_in,
                            const LatticeType &B_in) : /* init */ RationalMatrix<dim>(getTransitionMatrix(A_in, B_in))
    {
    }
    template struct LatticeTransitionMatrix<3>;

} // end namespace
#endif


