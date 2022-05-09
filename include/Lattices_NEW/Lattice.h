/* This file is part of model.
 *
 * model is distributed without any warranty under the MIT License.
 */


#ifndef model_Lattice_h_
#define model_Lattice_h_

#include <LatticeModule.h>
#include <StaticID.h>
#include <BestRationalApproximation.h>

namespace model
{
    template <int dim>
    class Lattice : public StaticID<Lattice<dim>>
    {
        static constexpr double roundTol=FLT_EPSILON;
        typedef typename LatticeCore<dim>::VectorDimD VectorDimD;
        typedef typename LatticeCore<dim>::MatrixDimD MatrixDimD;
        typedef typename LatticeCore<dim>::VectorDimI VectorDimI;
        typedef typename LatticeCore<dim>::MatrixDimI MatrixDimI;

        static MatrixDimD getLatticeBasis(const MatrixDimD& A,const MatrixDimD& Q);
        
    public:
        
        const MatrixDimD    latticeBasis;
        const MatrixDimD reciprocalBasis;
        const MatrixDimD C2G;

        Lattice(const MatrixDimD& A,const MatrixDimD& Q=MatrixDimD::Identity()) ;
        LatticeDirection<dim> latticeDirection(const VectorDimD& d) const;
        
        template<typename OtherDerived>
        LatticeDirection<dim> latticeDirection(const Eigen::MatrixBase<OtherDerived>& other) const
        {
            return LatticeDirection<dim>(LatticeVector<dim>(other,*this));
        }
        
        ReciprocalLatticeDirection<dim> reciprocalLatticeDirection(const VectorDimD& d) const;
        RationalLatticeDirection<dim> rationalLatticeDirection(const VectorDimD& d,
                                                              const typename BestRationalApproximation::LongIntType& maxDen=1000) const;
        
        template<typename OtherDerived>
        LatticeVector<dim> latticeVector(const Eigen::MatrixBase<OtherDerived>& other) const
        {
            return LatticeVector<dim>(other,*this);
        }
        
        LatticeVector<dim> latticeVector(const VectorDimD& p) const;
        ReciprocalLatticeVector<dim> reciprocalLatticeVector(const VectorDimD& p) const;

    };
}
#endif
