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
#include <model/LatticeMath/ReciprocalLatticeVector.h>

namespace model
{
    


    
    
    template <int dim>
    struct LatticeVector : public Eigen::Matrix<long int,dim,1>
    {
        static_assert(dim>0,"dim must be > 0.");
        
        typedef Eigen::Matrix<long int,dim,1> BaseType;
        typedef LatticeBase<dim> LatticeBaseType;
        typedef LatticeVector<dim> LatticeVectorType;
        typedef ReciprocalLatticeVector<dim> ReciprocalLatticeVectorType;
//        typedef Eigen::Matrix<long int,dim,1> VectorDimI;



        
    public:
        
        static constexpr double roundTol=FLT_EPSILON;
        
        typedef Eigen::Matrix<  double,dim,1> VectorDimD;
        typedef Eigen::Matrix<long int,dim,1> VectorDimI;

        typedef Eigen::Matrix<double,dim,dim> MatrixDimD;
        
        const MatrixDimD& covBasis;
//        const MatrixDimD& covBasisInv;
        
//        LatticeVector() : //BaseType(VectorDimI::Zero())
////        /* init base */ BaseType(BaseType::Zero())
//        /* init base */ covBasis(A)
//        {
//            
//        }
        
        /**********************************************************************/
        LatticeVector(const VectorDimD& d,
                      const MatrixDimD& A,
                      const MatrixDimD& invA) :
        /* init base */ BaseType(d2contra(d,invA)),
        /* init base */ covBasis(A)
//        /* init base */ covBasisInv(Ainv)
        ///* base init */ BaseType(LatticeBaseType::d2contra(d))
        {/*!@param[in] d vector in real space
          * Constructs *this by mapping d to the lattice
          */
            
        }
        
        /**********************************************************************/
        template<typename OtherDerived>
        LatticeVector(const Eigen::MatrixBase<OtherDerived>& other,
                      const MatrixDimD& A) :
        /* init base */ BaseType(other),
        /* init base */ covBasis(A)
        //        /* init base */ covBasisInv(Ainv)
        ///* base init */ BaseType(LatticeBaseType::d2contra(d))
        {/*!@param[in] d vector in real space
          * Constructs *this by mapping d to the lattice
          */
            
        }
        
        /**********************************************************************/
        LatticeVectorType operator+(const LatticeVectorType& other)
        {
            assert(&covBasis==&other.covBasis && "LatticeVectors have different bases.");
            return LatticeVectorType(static_cast<VectorDimI>(*this)+static_cast<VectorDimI>(other),covBasis);
        }
        
        /**********************************************************************/
        LatticeVectorType operator-(const LatticeVectorType& other)
        {
            assert(&covBasis==&other.covBasis && "LatticeVectors have different bases.");
            return LatticeVectorType(static_cast<VectorDimI>(*this)-static_cast<VectorDimI>(other),covBasis);
        }
        
        long int dot(const ReciprocalLatticeVectorType& other) const
        {
            return BaseType::dot(other);
        }

        /**********************************************************************/
        VectorDimD cartesian() const
        {
            return covBasis*this->template cast<double>();
        }
        
        /**********************************************************************/
        static VectorDimI d2contra(const VectorDimD& d,
                                   const MatrixDimD& invA)
        {
            const VectorDimD nd(invA*d);
            const VectorDimD rd(RoundEigen<double,dim>::round(nd));
            if((nd-rd).norm()>roundTol)
            {
                std::cout<<"d2contra, nd="<<nd.transpose()<<std::endl;
                std::cout<<"d2contra, rd="<<rd.transpose()<<std::endl;
                assert(0 && "Input vector is not a lattice vector");
            }
            return rd.template cast<long int>();
        }
        
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
        
    };
    
    
} // end namespace
#endif
