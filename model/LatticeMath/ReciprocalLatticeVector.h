/* This file is part of PIL, the Particle Interaction Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_ReciprocalLatticeVector_h_
#define model_ReciprocalLatticeVector_h_

#include <iostream>
#include <cfloat> // FLT_EPSILON
#include <Eigen/Dense>
//#include <RoundEigen.h>
#include <Lattice.h>

//#include <LatticeBase.h>



//template <int dim>
//class ReciprocalLatticeVector;

namespace model
{
    
    template <int dim>
    class LatticeVector;
    
    template <int dim>
    class ReciprocalLatticeVector : public Eigen::Matrix<long int,dim,1>
    {
        static_assert(dim>0,"dim must be > 0.");

        typedef Eigen::Matrix<long int,dim,1> BaseType;
        typedef Lattice<dim> LatticeType;
        
        /**********************************************************************/
        BaseType& base()
        {
            return *this;
        }
        
        /**********************************************************************/
        const BaseType& base() const
        {
            return *this;
        }
        
    public:
        
        static constexpr double roundTol=FLT_EPSILON;
        typedef Eigen::Matrix<  double,dim,1> VectorDimD;
        typedef Eigen::Matrix<long int,dim,1> VectorDimI;
        typedef Eigen::Matrix<double,dim,dim> MatrixDimD;
        typedef ReciprocalLatticeVector<dim> ReciprocalLatticeVectorType;
        typedef LatticeVector<dim> LatticeVectorType;
        typedef LatticeDirection<dim> LatticeDirectionType;

        const LatticeType& lattice;
        
        /**********************************************************************/
        ReciprocalLatticeVector(const LatticeType& lat) :
        /* init base */ BaseType(VectorDimI::Zero()),
        /* init      */ lattice(lat)
        ///* base init */ BaseType(LatticeBaseType::d2contra(d))
        {/*!@param[in] d vector in real space
          * Constructs *this by mapping d to the lattice
          */
            
        }
        
        /**********************************************************************/
        ReciprocalLatticeVector(const VectorDimD& d,
        /*                   */ const LatticeType& lat) :
        /* init base */ BaseType(d2cov(d,lat)),
        /* init base */ lattice(lat)
        {/*!@param[in] d vector in real space
          * Constructs *this by mapping d to the lattice
          */
            
        }
        
        /**********************************************************************/
        template<typename OtherDerived>
        ReciprocalLatticeVector(const Eigen::MatrixBase<OtherDerived>& other,
                                const LatticeType& lat) :
        /* init base */ BaseType(other),
        /* init base */ lattice(lat)
        {/*!@param[in] d vector in real space
          * Constructs *this by mapping d to the lattice
          */
        }
        
        /**********************************************************************/
        ReciprocalLatticeVector(const ReciprocalLatticeVectorType& other) = default;
        ReciprocalLatticeVector(ReciprocalLatticeVectorType&& other) =default;

        
        /**********************************************************************/
        ReciprocalLatticeVectorType& operator=(const ReciprocalLatticeVectorType& other)
        {
            assert(&lattice==&other.lattice && "ReciprocalLatticeVectorType belong to different Lattices.");
            base()=other.base();
            return *this;
        }
        
        /**********************************************************************/
        ReciprocalLatticeVectorType& operator=(ReciprocalLatticeVectorType&& other)
        {
            assert(&lattice==&other.lattice && "LatticeVectors belong to different Lattices.");
            base()=other.base();
            return *this;
        }
        
        /**********************************************************************/
        ReciprocalLatticeVectorType operator+(const ReciprocalLatticeVectorType& other) const
        {
            assert(&lattice==&other.lattice && "LatticeVectors belong to different Lattices.");
            return ReciprocalLatticeVectorType(static_cast<VectorDimI>(*this)+static_cast<VectorDimI>(other),lattice);

        }
        
        /**********************************************************************/
        ReciprocalLatticeVectorType& operator+=(const ReciprocalLatticeVectorType& other)
        {
            assert(&lattice==&other.lattice && "LatticeVectors belong to different Lattices.");
            base()+=other.base();
            return *this;
        }
        
        /**********************************************************************/
        ReciprocalLatticeVectorType operator-(const ReciprocalLatticeVectorType& other) const
        {
            assert(&lattice==&other.lattice && "LatticeVectors belong to different Lattices.");
            return ReciprocalLatticeVectorType(static_cast<VectorDimI>(*this)-static_cast<VectorDimI>(other),lattice);
        }
        
        /**********************************************************************/
        ReciprocalLatticeVectorType& operator-=(const ReciprocalLatticeVectorType& other)
        {
            assert(&lattice==&other.lattice && "LatticeVectors belong to different Lattices.");
            base()-=other.base();
            return *this;
        }
        
        /**********************************************************************/
        ReciprocalLatticeVectorType operator*(const long int& scalar) const
        {
            return ReciprocalLatticeVectorType(static_cast<VectorDimI>(*this)*scalar,lattice);
            
        }
        
        /**********************************************************************/
        long int dot(const LatticeVectorType& other) const
        {
            assert(&lattice==&other.lattice && "LatticeVectors belong to different Lattices.");
            return static_cast<VectorDimI>(*this).dot(static_cast<VectorDimI>(other));
        }
        
//        /**********************************************************************/
//        LatticeVectorType cross(const ReciprocalLatticeVectorType& other) const
//        {
//            assert(&lattice==&other.lattice && "LatticeVectors belong to different Lattices.");
//            return LatticeVectorType(static_cast<VectorDimI>(*this).cross(static_cast<VectorDimI>(other)),lattice);
//        }

        
        /**********************************************************************/
        LatticeDirectionType cross(const ReciprocalLatticeVectorType& other) const
        {
            assert(&lattice==&other.lattice && "LatticeVectors belong to different Lattices.");
            return LatticeDirectionType(LatticeVectorType(static_cast<VectorDimI>(*this).cross(static_cast<VectorDimI>(other)),lattice));
        }
        
        /**********************************************************************/
        VectorDimD cartesian() const
        {
            return lattice.reciprocalBasis*this->template cast<double>();
        }
        
        /**********************************************************************/
        double planeSpacing() const
        {
            return 1.0/cartesian().norm();
        }
        
        /**********************************************************************/
        VectorDimD interplaneVector() const
        {
            const VectorDimD c(cartesian());
            return c/c.squaredNorm();
        }
        
        /**********************************************************************/
        static VectorDimI d2cov(const VectorDimD& d,
                                const LatticeType& lat)
        {
            const VectorDimD nd(lat.latticeBasis.transpose()*d);
//            const VectorDimD rd(RoundEigen<double,dim>::round(nd));
            const VectorDimD rd(nd.array().round());
            if((nd-rd).norm()>roundTol)
            {
                std::cout<<"d2cov, nd="<<nd.transpose()<<std::endl;
                std::cout<<"d2cov, rd="<<rd.transpose()<<std::endl;
                assert(0 && "Input vector is not a reciprocal lattice vector");
            }
            //            assert((nd-rd).norm()<roundTol && "Input vector is not a lattice vector");
            return rd.template cast<long int>();
        }
        
    };
    
    template<int dim>
    ReciprocalLatticeVector<dim> operator*(const long int& scalar, const ReciprocalLatticeVector<dim>& L)
    {
        return L*scalar;
    }
    
} // end namespace
#endif
