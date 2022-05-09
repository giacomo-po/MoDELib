/* This file is part of model.
 *
 * model is distributed without any warranty under the MIT License.
 */

#ifndef model_ReciprocalLatticeVectorBase_h_
#define model_ReciprocalLatticeVectorBase_h_

#include <LatticeModule.h>

namespace model
{
    template <int dim>
    class ReciprocalLatticeVectorBase : public Eigen::Matrix<typename LatticeCore<dim>::IntScalarType,dim,1>
    {
        typedef Eigen::Matrix<typename LatticeCore<dim>::IntScalarType,dim,1> BaseType;
        BaseType& base();
        const BaseType& base() const;

    public:
        
        typedef typename LatticeCore<dim>::IntScalarType IntScalarType;
        typedef typename LatticeCore<dim>::VectorDimD VectorDimD;
        typedef typename LatticeCore<dim>::MatrixDimD MatrixDimD;
        typedef typename LatticeCore<dim>::VectorDimI VectorDimI;
        typedef typename LatticeCore<dim>::MatrixDimI MatrixDimI;

        const Lattice<dim>& lattice;
        
        ReciprocalLatticeVectorBase(const Lattice<dim>& lat);
        ReciprocalLatticeVectorBase(const VectorDimD& d, const Lattice<dim>& lat) ;
        template<typename OtherDerived>
        ReciprocalLatticeVectorBase(const Eigen::MatrixBase<OtherDerived>& other,
                                const Lattice<dim>& lat) :
        /* init base */ BaseType(other),
        /* init base */ lattice(lat)
        {/*!@param[in] d vector in real space
          * Constructs *this by mapping d to the lattice
          */
        }
        ReciprocalLatticeVectorBase(const ReciprocalLatticeVectorBase<dim>& other) = default;
        ReciprocalLatticeVectorBase(ReciprocalLatticeVectorBase<dim>&& other) =default;
        ReciprocalLatticeVectorBase<dim>& operator=(const ReciprocalLatticeVectorBase<dim>& other);
        ReciprocalLatticeVectorBase<dim>& operator=(ReciprocalLatticeVectorBase<dim>&& other);
        ReciprocalLatticeVectorBase<dim> operator+(const ReciprocalLatticeVectorBase<dim>& other) const;
        ReciprocalLatticeVectorBase<dim>& operator+=(const ReciprocalLatticeVectorBase<dim>& other);
        ReciprocalLatticeVectorBase<dim> operator-(const ReciprocalLatticeVectorBase<dim>& other) const;
        ReciprocalLatticeVectorBase<dim>& operator-=(const ReciprocalLatticeVectorBase<dim>& other);
        ReciprocalLatticeVectorBase<dim> operator*(const IntScalarType& scalar) const;
        IntScalarType dot(const LatticeVector<dim>& other) const;
        VectorDimD cartesian() const;
        double planeSpacing() const;
        VectorDimD interplaneVector() const;
        IntScalarType closestPlaneIndexOfPoint(const VectorDimD& P) const;
        IntScalarType planeIndexOfPoint(const VectorDimD& P) const ;
        IntScalarType planeIndexOfPoint(const LatticeVector<dim>& P) const;
        
    };
    
    template<int dim>
    ReciprocalLatticeVectorBase<dim> operator*(const typename ReciprocalLatticeVectorBase<dim>::IntScalarType& scalar, const ReciprocalLatticeVectorBase<dim>& L);
    


template <int dim>
class ReciprocalLatticeVector : public ReciprocalLatticeVectorBase<dim>
{

public:

    typedef typename LatticeCore<dim>::VectorDimD VectorDimD;
    typedef typename LatticeCore<dim>::MatrixDimD MatrixDimD;
    typedef typename LatticeCore<dim>::VectorDimI VectorDimI;
    typedef typename LatticeCore<dim>::MatrixDimI MatrixDimI;

    ReciprocalLatticeVector(const Lattice<dim>& lat) ;
    ReciprocalLatticeVector(const VectorDimD& d,const Lattice<dim>& lat) ;
    
    template<typename OtherDerived>
    ReciprocalLatticeVector(const Eigen::MatrixBase<OtherDerived>& other,
                  const Lattice<dim>& lat) :
    /* init base */ ReciprocalLatticeVectorBase<dim>(other,lat)
    { /*!@param[in] d vector in real space
       * Constructs *this by mapping d to the lattice
       */
    }
    
};


template <>
class ReciprocalLatticeVector<3> : public ReciprocalLatticeVectorBase<3>
{
    
public:

    
    typedef typename LatticeCore<3>::VectorDimD VectorDimD;
    typedef typename LatticeCore<3>::MatrixDimD MatrixDimD;
    typedef typename LatticeCore<3>::VectorDimI VectorDimI;
    typedef typename LatticeCore<3>::MatrixDimI MatrixDimI;

    
    ReciprocalLatticeVector(const Lattice<3>& lat) ;
    ReciprocalLatticeVector(const VectorDimD& d, const Lattice<3>& lat) ;
    
    template<typename OtherDerived>
    ReciprocalLatticeVector(const Eigen::MatrixBase<OtherDerived>& other,
                  const Lattice<3>& lat) :
    /* init base */ ReciprocalLatticeVectorBase<3>(other,lat)
    { /*!@param[in] d vector in real space
       * Constructs *this by mapping d to the lattice
       */
    }
    
    LatticeDirection<3> cross(const ReciprocalLatticeVectorBase<3>& other) const;

};


} // end namespace
#endif
