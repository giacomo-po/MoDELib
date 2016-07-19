/* This file is part of PIL, the Particle Interaction Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_ReciprocalLatticeVector_h_
#define model_ReciprocalLatticeVector_h_

#include <cfloat> // FLT_EPSILON
#include <Eigen/Dense>
#include <model/Math/RoundEigen.h>

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
//        typedef LatticeBase<dim> LatticeBaseType;
//        typedef Eigen::Matrix<  double,dim,1> VectorDimD;
        
    public:
        
        static constexpr double roundTol=FLT_EPSILON;
        typedef Eigen::Matrix<  double,dim,1> VectorDimD;
        typedef Eigen::Matrix<long int,dim,1> VectorDimI;
        typedef Eigen::Matrix<double,dim,dim> MatrixDimD;
        typedef LatticeVector<dim> LatticeVectorType;
        typedef ReciprocalLatticeVector<dim> ReciprocalLatticeVectorType;

        const MatrixDimD& covBasis;
        const MatrixDimD& contraBasis;

        
//        // Empty constructor
//        ReciprocalLatticeVector() :
//        /* init base */ BaseType(BaseType::Zero())
//        {
//            
//        }
        
        /**********************************************************************/
        ReciprocalLatticeVector(const VectorDimD& d,
                      const MatrixDimD& covBasis_in,
                      const MatrixDimD& contraBasis_in) :
        //                      const MatrixDimD& invA) :
        /* init base */ BaseType(d2cov(d,covBasis_in)),
        /* init base */ covBasis(covBasis_in),
        /* init base */ contraBasis(contraBasis_in)
        //        /* init base */ covBasisInv(Ainv)
        ///* base init */ BaseType(LatticeBaseType::d2contra(d))
        {/*!@param[in] d vector in real space
          * Constructs *this by mapping d to the lattice
          */
            
        }
        
        /**********************************************************************/
        template<typename OtherDerived>
        ReciprocalLatticeVector(const Eigen::MatrixBase<OtherDerived>& other,
                      const MatrixDimD& covBasis_in,
                      const MatrixDimD& contraBasis_in) :
        /* init base */ BaseType(other),
        /* init base */ covBasis(covBasis_in),
        /* init base */ contraBasis(contraBasis_in)
        //        /* init base */ covBasisInv(Ainv)
        ///* base init */ BaseType(LatticeBaseType::d2contra(d))
        {/*!@param[in] d vector in real space
          * Constructs *this by mapping d to the lattice
          */
            
        }
        
        /**********************************************************************/
        ReciprocalLatticeVectorType operator+(const ReciprocalLatticeVectorType& other) const
        {
            assert(&covBasis==&other.covBasis && "LatticeVectors have different bases.");
            assert(&contraBasis==&other.contraBasis && "LatticeVectors have different bases.");
            return ReciprocalLatticeVectorType(static_cast<VectorDimI>(*this)+static_cast<VectorDimI>(other),covBasis,contraBasis);
        }
        
        /**********************************************************************/
        ReciprocalLatticeVectorType operator-(const ReciprocalLatticeVectorType& other) const
        {
            assert(&covBasis==&other.covBasis && "LatticeVectors have different bases.");
            assert(&contraBasis==&other.contraBasis && "LatticeVectors have different bases.");
            return ReciprocalLatticeVectorType(static_cast<VectorDimI>(*this)-static_cast<VectorDimI>(other),covBasis,contraBasis);
        }
        
        /**********************************************************************/
        VectorDimD cartesian() const
        {
            return contraBasis*this->template cast<double>();
        }
        
        
        /**********************************************************************/
        long int dot(const LatticeVectorType& other) const
        {
            assert(&covBasis==&other.covBasis && "LatticeVectors have different bases.");
            assert(&contraBasis==&other.contraBasis && "LatticeVectors have different bases.");
            return static_cast<VectorDimI>(*this).dot(static_cast<VectorDimI>(other));
        }
        
        /**********************************************************************/
        LatticeVectorType cross(const ReciprocalLatticeVectorType& other) const
        {
            assert(&covBasis==&other.covBasis && "LatticeVectors have different bases.");
            assert(&contraBasis==&other.contraBasis && "LatticeVectors have different bases.");
            return LatticeVectorType(static_cast<VectorDimI>(*this).cross(static_cast<VectorDimI>(other)),covBasis,covBasis);
        }
        
        
        /**********************************************************************/
        static VectorDimI d2cov(const VectorDimD& d,
                                const MatrixDimD& covBasis)
        {
            const VectorDimD nd(covBasis.transpose()*d);
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
    
} // end namespace
#endif
