/* This file is part of PIL, the Particle Interaction Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LatticeVector_h_
#define model_LatticeVector_h_

#include <cfloat> // FLT_EPSILON
#include <Eigen/Dense>
#include <model/Math/RoundEigen.h>
//#include <model/LatticeMath/LatticeBase.h>
//#include <model/LatticeMath/ReciprocalLatticeVector.h>
#include <model/LatticeMath/Lattice.h>
#include <model/LatticeMath/ReciprocalLatticeVector.h>


namespace model
{
    
    template <int dim>
    struct LatticeVector : public Eigen::Matrix<long int,dim,1>
    {
        static_assert(dim>0,"dim must be > 0.");

        typedef Lattice<dim> LatticeType;
        typedef Eigen::Matrix<long int,dim,1> BaseType;
//        typedef LatticeBase<dim> LatticeBaseType;
        typedef LatticeVector<dim> LatticeVectorType;
        typedef ReciprocalLatticeVector<dim> ReciprocalLatticeVectorType;
//        typedef Eigen::Matrix<long int,dim,1> VectorDimI;

        Eigen::Matrix<long int,dim,1>& base()
        {
            return *this;
        }
        
        const Eigen::Matrix<long int,dim,1>& base() const
        {
            return *this;
        }

        
    public:
        
        static constexpr double roundTol=FLT_EPSILON;
        typedef Eigen::Matrix<  double,dim,1> VectorDimD;
        typedef Eigen::Matrix<long int,dim,1> VectorDimI;
        typedef Eigen::Matrix<double,dim,dim> MatrixDimD;
        
        const LatticeType& lattice;
//        const MatrixDimD& covBasis;
//        const MatrixDimD& contraBasis;

        /**********************************************************************/
        LatticeVector(const LatticeType& lat) :
        /* init base */ BaseType(VectorDimI::Zero()),
        /* init      */ lattice(lat)
        ///* base init */ BaseType(LatticeBaseType::d2contra(d))
        {/*!@param[in] d vector in real space
          * Constructs *this by mapping d to the lattice
          */
            
        }
        
        /**********************************************************************/
        LatticeVector(const VectorDimD& d,
                      const LatticeType& lat) :
//                      const MatrixDimD& invA) :
        /* init base */ BaseType(d2contra(d,lat)),
        /* init      */ lattice(lat)
//        /* init      */ covBasis(covBasis_in),
//        /* init      */ contraBasis(contraBasis_in)
//        /* init base */ covBasisInv(Ainv)
        ///* base init */ BaseType(LatticeBaseType::d2contra(d))
        {/*!@param[in] d vector in real space
          * Constructs *this by mapping d to the lattice
          */
            
        }
        
        /**********************************************************************/
        template<typename OtherDerived>
        LatticeVector(const Eigen::MatrixBase<OtherDerived>& other,
                      const LatticeType& lat
//                      const MatrixDimD& covBasis_in,
//                      const MatrixDimD& contraBasis_in
                      ) :
        /* init base */ BaseType(other),
        /* init      */ lattice(lat)
//
//        /* init base */ covBasis(covBasis_in),
//        /* init base */ contraBasis(contraBasis_in)
        //        /* init base */ covBasisInv(Ainv)
        ///* base init */ BaseType(LatticeBaseType::d2contra(d))
        {/*!@param[in] d vector in real space
          * Constructs *this by mapping d to the lattice
          */
//            std::cout<<"LatticeVector"<<std::endl;
//            model::cout<<"&covBasis="<<&covBasis<<std::endl;
//            model::cout<<"&contraBasis="<<&contraBasis<<std::endl;

        }
        
        /**********************************************************************/
        LatticeVector(const LatticeVectorType& other) = default;
        LatticeVector(LatticeVectorType&& other) =default;
        

        
        /**********************************************************************/
        LatticeVectorType& operator=(const LatticeVectorType& other)
        {
            assert(&lattice==&other.lattice && "LatticeVectors belong to different Lattices.");
//            assert(&covBasis==&other.covBasis && "LatticeVectors have different bases.");
//            assert(&contraBasis==&other.contraBasis && "LatticeVectors have different bases.");
            base()=other.base();
            return *this;
        }
        
        /**********************************************************************/
        LatticeVectorType& operator=(LatticeVectorType&& other)
        {
            assert(&lattice==&other.lattice && "LatticeVectors belong to different Lattices.");
//            assert(&covBasis==&other.covBasis && "LatticeVectors have different bases.");
//            assert(&contraBasis==&other.contraBasis && "LatticeVectors have different bases.");
            base()=other.base();
            return *this;
        }
        
        /**********************************************************************/
        LatticeVectorType operator+(const LatticeVectorType& other) const
        {
            assert(&lattice==&other.lattice && "LatticeVectors belong to different Lattices.");
//            assert(&covBasis==&other.covBasis && "LatticeVectors have different bases.");
//            assert(&contraBasis==&other.contraBasis && "LatticeVectors have different bases.");
//            return LatticeVectorType(static_cast<VectorDimI>(*this)+static_cast<VectorDimI>(other),covBasis,contraBasis);
                        return LatticeVectorType(static_cast<VectorDimI>(*this)+static_cast<VectorDimI>(other),lattice);
        }
        
        /**********************************************************************/
        LatticeVectorType& operator+=(const LatticeVectorType& other)
        {
            assert(&lattice==&other.lattice && "LatticeVectors belong to different Lattices.");
//            assert(&covBasis==&other.covBasis && "LatticeVectors have different bases.");
//            assert(&contraBasis==&other.contraBasis && "LatticeVectors have different bases.");
            base()+=other.base();
            return *this;
        }
        
        /**********************************************************************/
        LatticeVectorType operator-(const LatticeVectorType& other) const
        {
            assert(&lattice==&other.lattice && "LatticeVectors belong to different Lattices.");
//            assert(&covBasis==&other.covBasis && "LatticeVectors have different bases.");
//            assert(&contraBasis==&other.contraBasis && "LatticeVectors have different bases.");
//            return LatticeVectorType(static_cast<VectorDimI>(*this)-static_cast<VectorDimI>(other),covBasis,contraBasis);
            return LatticeVectorType(static_cast<VectorDimI>(*this)-static_cast<VectorDimI>(other),lattice);
        }
        
        /**********************************************************************/
        LatticeVectorType& operator-=(const LatticeVectorType& other)
        {
            assert(&lattice==&other.lattice && "LatticeVectors belong to different Lattices.");
//            assert(&covBasis==&other.covBasis && "LatticeVectors have different bases.");
//            assert(&contraBasis==&other.contraBasis && "LatticeVectors have different bases.");
            base()-=other.base();
            return *this;
        }
        
        /**********************************************************************/
        LatticeVectorType operator*(const long int& scalar) const
        {
//            return LatticeVectorType(static_cast<VectorDimI>(*this)*scalar,covBasis,contraBasis);
            return LatticeVectorType(static_cast<VectorDimI>(*this)*scalar,lattice);

        }
        
        /**********************************************************************/
        long int dot(const ReciprocalLatticeVectorType& other) const
        {
            assert(&lattice==&other.lattice && "LatticeVectors belong to different Lattices.");
//            assert(&covBasis==&other.covBasis && "LatticeVectors have different bases.");
//            assert(&contraBasis==&other.contraBasis && "LatticeVectors have different bases.");
            return static_cast<VectorDimI>(*this).dot(static_cast<VectorDimI>(other));
        }
        
        /**********************************************************************/
        ReciprocalLatticeVectorType cross(const LatticeVectorType& other) const
        {
            assert(&lattice==&other.lattice && "LatticeVectors belong to different Lattices.");
//            assert(&covBasis==&other.covBasis && "LatticeVectors have different bases.");
//            assert(&contraBasis==&other.contraBasis && "LatticeVectors have different bases.");
//            return ReciprocalLatticeVectorType(static_cast<VectorDimI>(*this).cross(static_cast<VectorDimI>(other)),covBasis,contraBasis);
            return ReciprocalLatticeVectorType(static_cast<VectorDimI>(*this).cross(static_cast<VectorDimI>(other)),lattice);
        }

        /**********************************************************************/
        VectorDimD cartesian() const
        {
//            return covBasis*this->template cast<double>();
            return lattice.covBasis()*this->template cast<double>();

        }
        
        /**********************************************************************/
        static VectorDimI d2contra(const VectorDimD& d,
                                   const LatticeType& lat)
        {
            const VectorDimD nd(lat.contraBasis().transpose()*d);
            const VectorDimD rd(RoundEigen<double,dim>::round(nd));
            if((nd-rd).norm()>roundTol)
            {
                std::cout<<"d2contra, nd="<<nd.transpose()<<std::endl;
                std::cout<<"d2contra, rd="<<rd.transpose()<<std::endl;
                assert(0 && "Input vector is not a lattice vector");
            }
            return rd.template cast<long int>();
        }
        
    };
    
    template<int dim>
    LatticeVector<dim> operator*(const long int& scalar, const LatticeVector<dim>& L)
    {
        return L*scalar;
    }
    
} // end namespace
#endif


//        /**********************************************************************/
//        LatticeVector(const LatticeVectorType& other) :
//        /* base copy */ BaseType(static_cast<VectorDimI>(other))
//        {
//            assert(&covBasis==&other.covBasis && "LatticeVectors have different bases.");
//            assert(&contraBasis==&other.contraBasis && "LatticeVectors have different bases.");
//        }
//
//        /**********************************************************************/
//        LatticeVector(LatticeVectorType&& other) :
//        /* base copy */ BaseType(std::move(static_cast<VectorDimI>(other)))
//        {
//            assert(&covBasis==&other.covBasis && "LatticeVectors have different bases.");
//            assert(&contraBasis==&other.contraBasis && "LatticeVectors have different bases.");
//        }


//        template<typename OtherDerived>
//        LatticeVector(const Eigen::MatrixBase<OtherDerived>& other) :
//        /* base init */ BaseType(other)
//        {// This constructor allows  to construct LatticeVector from Eigen expressions
//
//
//        }

//        template<typename OtherDerived>
//        LatticeVector& operator=(const Eigen::MatrixBase <OtherDerived>& other)
//        {// This method allows to assign Eigen expressions to LatticeVector
//            BaseType::operator=(other);
//            return *this;
//        }

//        VectorDimD cartesian() const
//        {
//            return LatticeBaseType::covBasis()*this->template cast<double>();
//        }
