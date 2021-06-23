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
        static IntValueType gcd(const IntValueType& a,const IntValueType& b);

        
        /**********************************************************************/
        static RationalMatrix<dim> getTransitionMatrix(const LatticeType& A,
                                                       const LatticeType& B);
        
        
        /**********************************************************************/
        LatticeTransitionMatrix(const LatticeType& A_in,
            const LatticeType& B_in) ;
    };
    
    
} // end namespace
#endif


