/* This file is part of PIL, the Particle Interaction Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_CSL_h_
#define model_CSL_h_

#include <cfloat> // FLT_EPSILON
#include <Eigen/Dense>
#include <model/Math/RoundEigen.h>
#include <model/Math/SmithDecomposition.h>
//#include <model/Math/BestRationalApproximation.h>
#include <model/LatticeMath/Lattice.h>
#include <model/LatticeMath/RationalRotation.h>


namespace model
{
    
    template <int dim>
    class CSL : public Lattice<dim>
    {
        static_assert(dim>0,"dim must be > 0.");
        //        static constexpr double roundTol=FLT_EPSILON;
        typedef Eigen::Matrix<  double,dim,1> VectorDimD;
        typedef Eigen::Matrix<long int,dim,1> VectorDimI;
        typedef Eigen::Matrix<double,dim,dim> MatrixDimD;
        typedef Lattice<dim> LatticeType;
        
        
        int _sigma;

    public:
        
        const LatticeType& A;
        const LatticeType& B;
        
        /**********************************************************************/
        CSL(const LatticeType& A_in,
            const LatticeType& B_in) :
        /* init */ _sigma(1),
        /* init */ A(A_in),
        /* init */ B(B_in)
        {
            update();
        }
        
        /**********************************************************************/
        void update()
        {
            RationalRotation<dim> rr(B.covBasis()*A.contraBasis().transpose());
            _sigma=rr.sigma();
            SmithDecomposition<dim> sd(rr.integerMatrix());
        }

        /**********************************************************************/
        const long int& sigma() const
        {
            return _sigma;
        }
    };
    
    
} // end namespace
#endif


