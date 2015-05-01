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
#include <model/LatticeMath/LatticeBase.h>


namespace model
{
    template <int dim>
    struct ReciprocalLatticeVector : public Eigen::Matrix<long int,dim,1>
    {
        static_assert(dim>0,"dim must be > 0.");

        typedef Eigen::Matrix<long int,dim,1> BaseType;
        typedef LatticeBase<dim> LatticeBaseType;
        typedef LatticeVector<dim> LatticeVectorType;
        typedef Eigen::Matrix<  double,dim,1> VectorDimD;
        
    public:
        

        
        // Empty constructor
        ReciprocalLatticeVector() :
        /* init base */ BaseType(BaseType::Zero())
        {
            
        }
        
        ReciprocalLatticeVector(const VectorDimD& d) :
        /* init */ BaseType(LatticeBaseType::d2cov(d))
        {/*!@param[in] d vector in real space
          * Constructs *this by mapping d to the lattice
          */
            
        }
        
        // This constructor allows  to construct LatticeVector from Eigen expressions
        template<typename OtherDerived>
        ReciprocalLatticeVector(const Eigen::MatrixBase<OtherDerived>& other) :
        /* base init */ BaseType(other)
        { }
        
        // This method allows to assign Eigen expressions to LatticeVector
        template<typename OtherDerived>
        ReciprocalLatticeVector& operator=(const Eigen::MatrixBase <OtherDerived>& other)
        {
            BaseType::operator=(other);
            return *this;
        }
        
//        ReciprocalLatticeVector& operator=(const ReciprocalLatticeVector& other)
//        {
//            BaseType::operator=(other);
//            return *this;
//        }
        
//        const VectorDimI& cov() const
//        {
//            return coVariant;
//        }
        
//        long int dot(const ReciprocalLatticeVectorType& other) const
//        {
//            return contra().dot(other.cov());
//        }
        
        long int dot(const LatticeVectorType& other) const
        {
            return BaseType::dot(other);
        }
        
        VectorDimD cartesian() const
        {
            return LatticeBaseType::contraBasis()*this->template cast<double>();
        }
        
        
    };
    
} // end namespace
#endif
