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
#include <model/Math/RoundEigen.h>
#include <model/LatticeMath/Lattice.h>

//#include <model/LatticeMath/LatticeBase.h>



//template <int dim>
//class ReciprocalLatticeVector;

namespace model
{
    
    template <int dim>
    class LatticeVector;
    
    template <int dim>
    struct ReciprocalLatticeVector : public Eigen::Matrix<long int,dim,1>
    {
        static_assert(dim>0,"dim must be > 0.");

        typedef Eigen::Matrix<long int,dim,1> BaseType;
        typedef Lattice<dim> LatticeType;
        

//        typedef LatticeBase<dim> LatticeBaseType;
//        typedef Eigen::Matrix<  double,dim,1> VectorDimD;
        
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
        typedef LatticeVector<dim> LatticeVectorType;
        typedef ReciprocalLatticeVector<dim> ReciprocalLatticeVectorType;

        const LatticeType& lattice;

//        const MatrixDimD& covBasis;
//        const MatrixDimD& contraBasis;

        
//        // Empty constructor
//        ReciprocalLatticeVector() :
//        /* init base */ BaseType(BaseType::Zero())
//        {
//            
//        }
        
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
//                      const MatrixDimD& covBasis_in,
//                      const MatrixDimD& contraBasis_in
        const LatticeType& lat) :
        //                      const MatrixDimD& invA) :
//        /* init base */ BaseType(d2cov(d,covBasis_in)),
        /* init base */ BaseType(d2cov(d,lat)),
        /* init base */ lattice(lat)
//        /* init base */ covBasis(covBasis_in),
//        /* init base */ contraBasis(contraBasis_in)
        //        /* init base */ covBasisInv(Ainv)
        ///* base init */ BaseType(LatticeBaseType::d2contra(d))
        {/*!@param[in] d vector in real space
          * Constructs *this by mapping d to the lattice
          */
            
        }
        
        /**********************************************************************/
        template<typename OtherDerived>
        ReciprocalLatticeVector(const Eigen::MatrixBase<OtherDerived>& other,
                                const LatticeType& lat) :
//                      const MatrixDimD& covBasis_in,
//                      const MatrixDimD& contraBasis_in) :
        /* init base */ BaseType(other),
        /* init base */ lattice(lat)
//        /* init base */ covBasis(covBasis_in),
//        /* init base */ contraBasis(contraBasis_in)
        //        /* init base */ covBasisInv(Ainv)
        ///* base init */ BaseType(LatticeBaseType::d2contra(d))
        {/*!@param[in] d vector in real space
          * Constructs *this by mapping d to the lattice
          */
//            std::cout<<"Creating ReciprocalLatticeVector "<<this<<std::endl;
//            std::cout<<"&covBasis="<<&covBasis_in<<std::endl;
//            std::cout<<"&contraBasis="<<&contraBasis_in<<std::endl;
        }
        
        /**********************************************************************/
        ReciprocalLatticeVector(const ReciprocalLatticeVectorType& other) = default;
        ReciprocalLatticeVector(ReciprocalLatticeVectorType&& other) =default;

        
        /**********************************************************************/
        ReciprocalLatticeVectorType& operator=(const ReciprocalLatticeVectorType& other)
        {
            assert(&lattice==&other.lattice && "ReciprocalLatticeVectorType belong to different Lattices.");
            //            assert(&covBasis==&other.covBasis && "LatticeVectors have different bases.");
            //            assert(&contraBasis==&other.contraBasis && "LatticeVectors have different bases.");
            base()=other.base();
            return *this;
        }
        
        /**********************************************************************/
        ReciprocalLatticeVectorType& operator=(ReciprocalLatticeVectorType&& other)
        {
            assert(&lattice==&other.lattice && "LatticeVectors belong to different Lattices.");
            //            assert(&covBasis==&other.covBasis && "LatticeVectors have different bases.");
            //            assert(&contraBasis==&other.contraBasis && "LatticeVectors have different bases.");
            base()=other.base();
            return *this;
        }
        
        /**********************************************************************/
        ReciprocalLatticeVectorType operator+(const ReciprocalLatticeVectorType& other) const
        {
            assert(&lattice==&other.lattice && "LatticeVectors belong to different Lattices.");
//            assert(&covBasis==&other.covBasis && "LatticeVectors have different bases.");
//            assert(&contraBasis==&other.contraBasis && "LatticeVectors have different bases.");
//            return ReciprocalLatticeVectorType(static_cast<VectorDimI>(*this)+static_cast<VectorDimI>(other),covBasis,contraBasis);
            return ReciprocalLatticeVectorType(static_cast<VectorDimI>(*this)+static_cast<VectorDimI>(other),lattice);

        }
        
        /**********************************************************************/
        ReciprocalLatticeVectorType& operator+=(const ReciprocalLatticeVectorType& other)
        {
            assert(&lattice==&other.lattice && "LatticeVectors belong to different Lattices.");
            //            assert(&covBasis==&other.covBasis && "LatticeVectors have different bases.");
            //            assert(&contraBasis==&other.contraBasis && "LatticeVectors have different bases.");
            base()+=other.base();
            return *this;
        }
        
        /**********************************************************************/
        ReciprocalLatticeVectorType operator-(const ReciprocalLatticeVectorType& other) const
        {
            assert(&lattice==&other.lattice && "LatticeVectors belong to different Lattices.");
//            assert(&covBasis==&other.covBasis && "LatticeVectors have different bases.");
//            assert(&contraBasis==&other.contraBasis && "LatticeVectors have different bases.");
//            return ReciprocalLatticeVectorType(static_cast<VectorDimI>(*this)-static_cast<VectorDimI>(other),covBasis,contraBasis);
            return ReciprocalLatticeVectorType(static_cast<VectorDimI>(*this)-static_cast<VectorDimI>(other),lattice);

        }
        
        /**********************************************************************/
        ReciprocalLatticeVectorType& operator-=(const ReciprocalLatticeVectorType& other)
        {
            assert(&lattice==&other.lattice && "LatticeVectors belong to different Lattices.");
            //            assert(&covBasis==&other.covBasis && "LatticeVectors have different bases.");
            //            assert(&contraBasis==&other.contraBasis && "LatticeVectors have different bases.");
            base()-=other.base();
            return *this;
        }
        
        /**********************************************************************/
        ReciprocalLatticeVectorType operator*(const long int& scalar) const
        {
            //            return LatticeVectorType(static_cast<VectorDimI>(*this)*scalar,covBasis,contraBasis);
            return ReciprocalLatticeVectorType(static_cast<VectorDimI>(*this)*scalar,lattice);
            
        }

        

        
        
        /**********************************************************************/
        long int dot(const LatticeVectorType& other) const
        {
//            std::cout<<"ReciprocalLatticeVector "<<this<<std::endl;
//            std::cout<<"&covBasis="<<&covBasis<<std::endl;
//            std::cout<<"&contraBasis="<<&contraBasis<<std::endl;
//            std::cout<<"&other.covBasis="<<&other.covBasis<<std::endl;
//            std::cout<<"&other.contraBasis="<<&other.contraBasis<<std::endl;
            assert(&lattice==&other.lattice && "LatticeVectors belong to different Lattices.");
//            assert(&covBasis==&other.covBasis && "LatticeVectors have different bases.");
//            assert(&contraBasis==&other.contraBasis && "LatticeVectors have different bases.");
            return static_cast<VectorDimI>(*this).dot(static_cast<VectorDimI>(other));
        }
        
        /**********************************************************************/
        LatticeVectorType cross(const ReciprocalLatticeVectorType& other) const
        {
            assert(&lattice==&other.lattice && "LatticeVectors belong to different Lattices.");
//            assert(&covBasis==&other.covBasis && "LatticeVectors have different bases.");
//            assert(&contraBasis==&other.contraBasis && "LatticeVectors have different bases.");
//            return LatticeVectorType(static_cast<VectorDimI>(*this).cross(static_cast<VectorDimI>(other)),covBasis,contraBasis);
            return LatticeVectorType(static_cast<VectorDimI>(*this).cross(static_cast<VectorDimI>(other)),lattice);

        }
        
        /**********************************************************************/
        VectorDimD cartesian() const
        {
            return lattice.contraBasis()*this->template cast<double>();
        }
        
        /**********************************************************************/
        static VectorDimI d2cov(const VectorDimD& d,
                                const LatticeType& lat)
        {
            const VectorDimD nd(lat.covBasis().transpose()*d);
            const VectorDimD rd(RoundEigen<double,dim>::round(nd));
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
