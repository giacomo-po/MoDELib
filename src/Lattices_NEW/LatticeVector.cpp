/* This file is part of model.
 *
 * model is distributed without any warranty under the MIT License.
 */


#ifndef model_LatticeVector_cpp_
#define model_LatticeVector_cpp_

#include<LatticeModule.h>

namespace model
{

    /**********************************************************************/
    template <int dim>
    typename LatticeVectorBase<dim>::BaseType& LatticeVectorBase<dim>::base()
    {
        return *this;
    }

    /**********************************************************************/
    template <int dim>
    const typename LatticeVectorBase<dim>::BaseType& LatticeVectorBase<dim>::base() const
    {
        return *this;
    }


    /**********************************************************************/
    template <int dim>
    LatticeVectorBase<dim>::LatticeVectorBase(const Lattice<dim> &lat) :
    /* init */ BaseType(VectorDimI::Zero()),
    /* init */ lattice(lat)
    ///* base init */ BaseType(LatticeBaseType::d2contra(d))
    { /*!@param[in] d vector in real space
          * Constructs *this by mapping d to the lattice
          */
    }

    /**********************************************************************/
    template <int dim>
    LatticeVectorBase<dim>::LatticeVectorBase(const VectorDimD &d,
                  const Lattice<dim> &lat) :
    /* init */ BaseType(LatticeCore<dim>::integerCoordinates(d,lat.reciprocalBasis.transpose())),
    /* init */ lattice(lat)
    { /*!@param[in] d vector in real space
          * Constructs *this by mapping d to the lattice
          */
    }

    /**********************************************************************/
    template <int dim>
    LatticeVectorBase<dim>& LatticeVectorBase<dim>::operator=(const LatticeVectorBase<dim> &other)
    {
        assert(&lattice == &other.lattice && "LatticeVectors belong to different Lattices.");
        base() = other.base();
        return *this;
    }

    /**********************************************************************/
    template <int dim>
    LatticeVectorBase<dim>& LatticeVectorBase<dim>::operator=(LatticeVectorBase<dim> &&other)
    {
        assert(&lattice == &other.lattice && "LatticeVectors belong to different Lattices.");
        base() = other.base();
        return *this;
    }

    /**********************************************************************/
    template <int dim>
    LatticeVectorBase<dim> LatticeVectorBase<dim>::operator+(const LatticeVectorBase<dim> &other) const
    {
        assert(&lattice == &other.lattice && "LatticeVectors belong to different Lattices.");
        return LatticeVectorBase<dim>(static_cast<VectorDimI>(*this) + static_cast<VectorDimI>(other), lattice);
    }

    /**********************************************************************/
    template <int dim>
    LatticeVectorBase<dim>& LatticeVectorBase<dim>::operator+=(const LatticeVectorBase<dim> &other)
    {
        assert(&lattice == &other.lattice && "LatticeVectors belong to different Lattices.");
        base() += other.base();
        return *this;
    }

    /**********************************************************************/
    template <int dim>
    LatticeVectorBase<dim> LatticeVectorBase<dim>::operator-(const LatticeVectorBase<dim> &other) const
    {
        assert(&lattice == &other.lattice && "LatticeVectors belong to different Lattices.");
        return LatticeVectorBase<dim>(static_cast<VectorDimI>(*this) - static_cast<VectorDimI>(other), lattice);
    }

    /**********************************************************************/
    template <int dim>
    LatticeVectorBase<dim>& LatticeVectorBase<dim>::operator-=(const LatticeVectorBase<dim> &other)
    {
        assert(&lattice == &other.lattice && "LatticeVectors belong to different Lattices.");
        base() -= other.base();
        return *this;
    }

    /**********************************************************************/
    template <int dim>
    typename LatticeVectorBase<dim>::IntScalarType LatticeVectorBase<dim>::dot(const ReciprocalLatticeVectorBase<dim> &other) const
    {
        assert(&lattice == &other.lattice && "LatticeVectors belong to different Lattices.");
        return static_cast<VectorDimI>(*this).dot(static_cast<VectorDimI>(other));
    }

    /**********************************************************************/
    template <int dim>
    typename LatticeVectorBase<dim>::VectorDimD LatticeVectorBase<dim>::cartesian() const
    {
        return lattice.latticeBasis * this->template cast<double>();
    }
    
    template<int dim>
    LatticeVectorBase<dim> operator*(const typename LatticeVectorBase<dim>::IntScalarType& scalar, const LatticeVectorBase<dim>& L)
    {
        return L*scalar;
    }

template <int dim>
LatticeVector<dim>::LatticeVector(const Lattice<dim> &lat) :
/* init base */ LatticeVectorBase<dim>(lat)
{ /*!@param[in] d vector in real space
      * Constructs *this by mapping d to the lattice
      */
}

template <int dim>
LatticeVector<dim>::LatticeVector(const VectorDimD &d,
              const Lattice<dim> &lat) :
/* init base */ LatticeVectorBase<dim>(d,lat)
{ /*!@param[in] d vector in real space
      * Constructs *this by mapping d to the lattice
      */
}

LatticeVector<3>::LatticeVector(const Lattice<3> &lat) :
/* init base */ LatticeVectorBase<3>(lat)
{ /*!@param[in] d vector in real space
      * Constructs *this by mapping d to the lattice
      */
}

LatticeVector<3>::LatticeVector(const VectorDimD &d,
              const Lattice<3> &lat) :
/* init base */ LatticeVectorBase<3>(d,lat)
{ /*!@param[in] d vector in real space
      * Constructs *this by mapping d to the lattice
      */
}

ReciprocalLatticeDirection<3> LatticeVector<3>::cross(const LatticeVector<3> &other) const
{
    assert(&lattice == &other.lattice && "LatticeVectors belong to different Lattices.");
    return ReciprocalLatticeDirection<3>(ReciprocalLatticeVector<3>(static_cast<VectorDimI>(*this).cross(static_cast<VectorDimI>(other)), lattice));
}

    template class LatticeVectorBase<1>;
    template LatticeVectorBase<1>operator*(const typename LatticeVectorBase<1>::IntScalarType& scalar, const LatticeVectorBase<1>& L);
    template class LatticeVectorBase<2>;
    template LatticeVectorBase<2>operator*(const typename LatticeVectorBase<2>::IntScalarType& scalar, const LatticeVectorBase<2>& L);
    template class LatticeVectorBase<3>;
    template LatticeVectorBase<3>operator*(const typename LatticeVectorBase<3>::IntScalarType& scalar, const LatticeVectorBase<3>& L);

    template class LatticeVector<1>;
    template class LatticeVector<2>;

} // end namespace
#endif
