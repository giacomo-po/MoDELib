/* This file is part of model.
 *
 * model is distributed without any warranty under the MIT License.
 */


#ifndef model_BiCrystal_h_
#define model_BiCrystal_h_

#include <LatticeModule.h>
#include <SmithDecomposition.h>
#include <RationalMatrix.h>
#include <LLL.h>
#include <RLLL.h>


namespace model
{
    /*!Class template that computes the Coincident Site Lattice (CSL) of two
     * parent lattices using the Smith decomposition method [1].
     *
     * [1] Coincidence Lattices and Associated Shear Transformations
     */
    template <int dim>
    class BiCrystal : public RationalMatrix<dim>
    /*             */,public SmithDecomposition<dim>
    {
        typedef typename LatticeCore<dim>::IntScalarType IntScalarType;
        typedef typename LatticeCore<dim>::VectorDimD VectorDimD;
        typedef typename LatticeCore<dim>::MatrixDimD MatrixDimD;
        typedef typename LatticeCore<dim>::VectorDimI VectorDimI;
        typedef typename LatticeCore<dim>::MatrixDimI MatrixDimI;

        
        static MatrixDimI getM(const RationalMatrix<dim>& rm, const SmithDecomposition<dim>& sd);
        static MatrixDimI getN(const RationalMatrix<dim>& rm, const SmithDecomposition<dim>& sd);
        static MatrixDimD getCSLBasis(const Lattice<dim>& A,
                                      const Lattice<dim>& B,
                                      const SmithDecomposition<dim>& sd,
                                      const MatrixDimI& M,
                                      const MatrixDimI& N,
                                      const bool& useRLLL);
        static MatrixDimD getDSCLBasis(const Lattice<dim>& A,
                                       const Lattice<dim>& B,
                                       const SmithDecomposition<dim>& sd,
                                       const MatrixDimI& M,
                                       const MatrixDimI& N,
                                       const bool& useRLLL);
        
    public:

        const Lattice<dim>& A;
        const Lattice<dim>& B;
        const MatrixDimI M;
        const MatrixDimI N;
        const int sigmaA;
        const int sigmaB;
        const int sigma;
        const Lattice<dim> csl;
        const Lattice<dim> dscl;
        const Lattice<dim> Ap;
        const Lattice<dim> Bp;
        
        /**********************************************************************/
        BiCrystal(const Lattice<dim>& A,
                  const Lattice<dim>& B,
                  const bool& useRLLL=true);
        
//        LatticeDirection<dim> AtoCSLvector(const LatticeVector<dim>& v) const;
        
    };
    
    
} // end namespace
#endif

