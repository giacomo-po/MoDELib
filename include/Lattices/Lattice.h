/* This file is part of PIL, the Particle Interaction Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_Lattice_h_
#define model_Lattice_h_

#include <iomanip> // FLT_EPSILON
#include <cfloat> // FLT_EPSILON
#include <Eigen/Dense>
#include <StaticID.h>
//#include <RoundEigen.h>
//#include <LatticeBase.h>
//#include <ReciprocalLatticeVector.h>
#include <LatticeVector.h>
#include <ReciprocalLatticeVector.h>
#include <LatticeDirection.h>
#include <ReciprocalLatticeDirection.h>
#include <BestRationalApproximation.h>

#include <RationalLatticeDirection.h>

namespace model
{

    /**************************************************************************/
    /**************************************************************************/
    template <int dim>
    class Lattice : public StaticID<Lattice<dim>>
    {
        static constexpr double roundTol=FLT_EPSILON;
        typedef Eigen::Matrix<  double,dim,1> VectorDimD;
        typedef Eigen::Matrix<long int,dim,1> VectorDimI;
        typedef Eigen::Matrix<double,dim,dim> MatrixDimD;
        typedef LatticeVector<dim> LatticeVectorType;
        typedef LatticeDirection<dim> LatticeDirectionType;
        typedef RationalLatticeDirection<dim> RationalLatticeDirectionType;
        typedef ReciprocalLatticeVector<dim> ReciprocalLatticeVectorType;
        typedef ReciprocalLatticeDirection<dim> ReciprocalLatticeDirectionType;

        static MatrixDimD getLatticeBasis(const MatrixDimD& A,const MatrixDimD& Q);

        
    public:
        
        const MatrixDimD    latticeBasis;
        const MatrixDimD reciprocalBasis;
        const MatrixDimD C2G;
        
        Lattice(const MatrixDimD& A,const MatrixDimD& Q);

        static Eigen::Matrix<long int,dim,1> rationalApproximation(VectorDimD nd);
        
        LatticeVectorType snapToLattice(const VectorDimD& d) const;
        
        LatticeDirectionType latticeDirection(const VectorDimD& d) const;
      
        ReciprocalLatticeDirectionType reciprocalLatticeDirection(const VectorDimD& d) const;
        
        RationalLatticeDirectionType rationalLatticeDirection(const VectorDimD& d,
                                                              const typename BestRationalApproximation::LongIntType& maxDen=1000) const;
        
        LatticeVectorType latticeVector(const VectorDimD& p) const;
        
        ReciprocalLatticeVectorType reciprocalLatticeVector(const VectorDimD& p) const;
        
    };

    
}
#endif

//        /**********************************************************************/
//        Lattice<dim,dim> reciprocal() const
//        {
//            return Lattice<dim,dim>(reciprocalBasis(),latticeBasis());
//        }

//        /**********************************************************************/
//        const MatrixDimD& latticeBasis() const __attribute__ ((deprecated)) // make latticeBasis public since it is const
//        {
//            return latticeBasis;
//        }
//
//        /**********************************************************************/
//        const MatrixDimD& reciprocalBasis() const __attribute__ ((deprecated)) // make reciprocalBasis public since it is const
//        {
//            return reciprocalBasis;
//        }
//
//        const VectorDimD basisVector(const size_t& c) const
//        {
//            assert(c<dim);
//            return latticeBasis().col(c);
//        }

//        /**********************************************************************/
//        void setLatticeBasis(const MatrixDimD& A)
//        {
//            latticeBases=getLatticeBases(A);
////            latticeBasis=A;
////            LatticeBaseType::setLatticeBasis(A);
////            reciprocalBasis=latticeBasis().inverse().transpose();
////            reciprocalBasis=latticeBasis.inverse().transpose();
////            std::cout<<"Lattice basis (in columns) =\n"<<latticeBasis<<std::endl;
////            std::cout<<"Lattice reciprocal basis (in columns) =\n"<<reciprocalBasis<<std::endl;
//        }

