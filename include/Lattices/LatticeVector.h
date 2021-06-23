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
//#include <RoundEigen.h>
//#include <LatticeBase.h>
//#include <ReciprocalLatticeVector.h>
#include <Lattice.h>
#include <ReciprocalLatticeVector.h>


namespace model
{
    
    template <int dim>
    class LatticeVector : public Eigen::Matrix<long int,dim,1>
    {
        static_assert(dim>0,"dim must be > 0.");
        
        typedef Lattice<dim> LatticeType;
        typedef Eigen::Matrix<long int,dim,1> BaseType;
        
        /**********************************************************************/
        BaseType& base();
        
        /**********************************************************************/
        const BaseType& base() const;
        
        
    public:
        
        static constexpr double roundTol=FLT_EPSILON;
        typedef Eigen::Matrix<  double,dim,1> VectorDimD;
        typedef Eigen::Matrix<long int,dim,1> VectorDimI;
        typedef Eigen::Matrix<double,dim,dim> MatrixDimD;
        typedef LatticeVector<dim> LatticeVectorType;
        typedef ReciprocalLatticeVector<dim> ReciprocalLatticeVectorType;
        typedef ReciprocalLatticeDirection<dim> ReciprocalLatticeDirectionType;
        
        const LatticeType& lattice;
        
        /**********************************************************************/
        LatticeVector(const LatticeType& lat) ;
        
        /**********************************************************************/
        LatticeVector(const VectorDimD& d,
                      const LatticeType& lat) ;
        
        /**********************************************************************/
        template<typename OtherDerived>
        LatticeVector(const Eigen::MatrixBase<OtherDerived>& other,
                      const LatticeType& lat) :
        /* init base */ BaseType(other),
        /* init      */ lattice(lat)
        { /*!@param[in] d vector in real space
           * Constructs *this by mapping d to the lattice
           */
        }
        
        /**********************************************************************/
        LatticeVector(const LatticeVectorType& other) = default;
        LatticeVector(LatticeVectorType&& other) =default;
        
        
        
        /**********************************************************************/
        LatticeVectorType& operator=(const LatticeVectorType& other);
        
        /**********************************************************************/
        LatticeVectorType& operator=(LatticeVectorType&& other);
        
        /**********************************************************************/
        LatticeVectorType operator+(const LatticeVectorType& other) const;
        
        /**********************************************************************/
        LatticeVectorType& operator+=(const LatticeVectorType& other);
        
        /**********************************************************************/
        LatticeVectorType operator-(const LatticeVectorType& other) const;
        
        /**********************************************************************/
        LatticeVectorType& operator-=(const LatticeVectorType& other);
        
        /**********************************************************************/
        LatticeVectorType operator*(const long int& scalar) const
        {
            return LatticeVectorType(static_cast<VectorDimI>(*this) * scalar, lattice);
        }
        
        /**********************************************************************/
        long int dot(const ReciprocalLatticeVectorType& other) const;
        
      
        /**********************************************************************/
        ReciprocalLatticeDirectionType cross(const LatticeVectorType& other) const;
        
        
        /**********************************************************************/
        VectorDimD cartesian() const;
        
        /**********************************************************************/
        static VectorDimI d2contra(const VectorDimD& d,
                                   const LatticeType& lat);
        
    };
        
    template<int dim>
    LatticeVector<dim> operator*(const long int& scalar, const LatticeVector<dim>& L);
    
    
} // end namespace
#endif
