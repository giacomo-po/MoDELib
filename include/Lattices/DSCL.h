/* This file is part of PIL, the Particle Interaction Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DSCL_h_
#define model_DSCL_h_

#include <cfloat> // FLT_EPSILON
#include <Eigen/Dense>
#include <SmithDecomposition.h>
#include <Lattice.h>
#include <LatticeTransitionMatrix.h>
#include <LLL.h>
#include <RLLL.h>


namespace model
{
    /*!Class template that computes the Displacement Shift Complete Lattice (DSCL)
     * of two parent lattices using the Smith decomposition method [1].
     *
     * [1] Coincidence Lattices and Associated Shear Transformations
     */
    template <int dim>
    class DSCL : public LatticeTransitionMatrix<dim>
    /*       */ ,public Lattice<dim>
    {
        static_assert(dim>0,"dim must be > 0.");
        typedef Eigen::Matrix<  double,dim,1> VectorDimD;
        typedef Eigen::Matrix<long int,dim,1> VectorDimI;
        typedef Eigen::Matrix<double,dim,dim> MatrixDimD;
        typedef Lattice<dim> LatticeType;
        typedef long long int IntValueType;
        typedef Eigen::Matrix<IntValueType,dim,dim> MatrixInt;
        
        /**********************************************************************/
        static MatrixDimD getLatticeBasis(const LatticeTransitionMatrix<dim>& rm,
                                          const LatticeType& A,
                                          const LatticeType& B,
                                          const bool& useRLLL);
        
    public:
        

        /**********************************************************************/
        DSCL(const LatticeType& A_in,
             const LatticeType& B_in,
             const bool& useRLLL=true) ;
        
    };
}
#endif