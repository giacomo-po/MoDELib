/* This file is part of model.
 *
 * model is distributed without any warranty under the MIT License.
 */

#ifndef model_LatticeVector_h_
#define model_LatticeVector_h_

#include <LatticeModule.h>

namespace model
{
    template <int dim>
    class LatticeVectorBase : public Eigen::Matrix<typename LatticeCore<dim>::IntScalarType,dim,1>
    {
        typedef Eigen::Matrix<typename LatticeCore<dim>::IntScalarType,dim,1> BaseType;
        BaseType& base();
        const BaseType& base() const;

    public:
        
//        static constexpr double roundTol=FLT_EPSILON;
        typedef typename LatticeCore<dim>::IntScalarType IntScalarType;
        typedef typename LatticeCore<dim>::VectorDimD VectorDimD;
        typedef typename LatticeCore<dim>::MatrixDimD MatrixDimD;
        typedef typename LatticeCore<dim>::VectorDimI VectorDimI;
        typedef typename LatticeCore<dim>::MatrixDimI MatrixDimI;
        
        const Lattice<dim>& lattice;
        
        LatticeVectorBase(const Lattice<dim>& lat) ;
        LatticeVectorBase(const VectorDimD& d, const Lattice<dim>& lat) ;
        template<typename OtherDerived>
        LatticeVectorBase(const Eigen::MatrixBase<OtherDerived>& other, const Lattice<dim>& lat) :
        /* init base */ BaseType(other),
        /* init      */ lattice(lat)
        { /*!@param[in] d vector in real space
           * Constructs *this by mapping d to the lattice
           */
        }
        LatticeVectorBase(const LatticeVectorBase<dim>& other) = default;
        LatticeVectorBase(LatticeVectorBase<dim>&& other) =default;
        LatticeVectorBase<dim>& operator=(const LatticeVectorBase<dim>& other);
        LatticeVectorBase<dim>& operator=(LatticeVectorBase<dim>&& other);
        LatticeVectorBase<dim> operator+(const LatticeVectorBase<dim>& other) const;
        LatticeVectorBase<dim>& operator+=(const LatticeVectorBase<dim>& other);
        LatticeVectorBase<dim> operator-(const LatticeVectorBase<dim>& other) const;
        LatticeVectorBase<dim>& operator-=(const LatticeVectorBase<dim>& other);
        LatticeVectorBase<dim> operator*(const IntScalarType& scalar) const
        {
            return LatticeVectorBase<dim>(static_cast<VectorDimI>(*this) * scalar, lattice);
        }
        IntScalarType dot(const ReciprocalLatticeVectorBase<dim>& other) const;
        VectorDimD cartesian() const;
        
    };
        
    template<int dim>
    LatticeVectorBase<dim> operator*(const typename LatticeVectorBase<dim>::IntScalarType& scalar, const LatticeVectorBase<dim>& L);
    

template <int dim>
class LatticeVector : public LatticeVectorBase<dim>
{
    
public:

    typedef typename LatticeCore<dim>::VectorDimD VectorDimD;
    typedef typename LatticeCore<dim>::MatrixDimD MatrixDimD;
    typedef typename LatticeCore<dim>::VectorDimI VectorDimI;
    typedef typename LatticeCore<dim>::MatrixDimI MatrixDimI;

    
    LatticeVector(const Lattice<dim>& lat) ;
    
    
    LatticeVector(const VectorDimD& d,
                  const Lattice<dim>& lat) ;
    
    
    template<typename OtherDerived>
    LatticeVector(const Eigen::MatrixBase<OtherDerived>& other,
                  const Lattice<dim>& lat) :
    /* init base */ LatticeVectorBase<dim>(other,lat)
    { /*!@param[in] d vector in real space
       * Constructs *this by mapping d to the lattice
       */
    }
    
};


template <>
class LatticeVector<3> : public LatticeVectorBase<3>
{
    
public:
    
    typedef typename LatticeCore<3>::VectorDimD VectorDimD;
    typedef typename LatticeCore<3>::MatrixDimD MatrixDimD;
    typedef typename LatticeCore<3>::VectorDimI VectorDimI;
    typedef typename LatticeCore<3>::MatrixDimI MatrixDimI;
    
    
    LatticeVector(const Lattice<3>& lat) ;
    
    
    LatticeVector(const VectorDimD& d,
                  const Lattice<3>& lat) ;
    
    
    template<typename OtherDerived>
    LatticeVector(const Eigen::MatrixBase<OtherDerived>& other,
                  const Lattice<3>& lat) :
    /* init base */ LatticeVectorBase<3>(other,lat)
    { /*!@param[in] d vector in real space
       * Constructs *this by mapping d to the lattice
       */
    }
    
    ReciprocalLatticeDirection<3> cross(const LatticeVector<3>& other) const;

    
};

    
} // end namespace
#endif
