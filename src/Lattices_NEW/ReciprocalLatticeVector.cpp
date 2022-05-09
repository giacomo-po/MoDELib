/* This file is part of model.
 *
 * model is distributed without any warranty under the MIT License.
 */


#ifndef model_ReciprocalLatticeVector_cpp_
#define model_ReciprocalLatticeVector_cpp_

#include <iostream>
#include <LatticeModule.h>

namespace model
{

    template <int dim>
    typename ReciprocalLatticeVectorBase<dim>::BaseType& ReciprocalLatticeVectorBase<dim>::base()
    {
        return *this;
    }

    template <int dim>
    const typename ReciprocalLatticeVectorBase<dim>::BaseType& ReciprocalLatticeVectorBase<dim>::base() const
    {
        return *this;
    }
    
    template <int dim>
    ReciprocalLatticeVectorBase<dim>::ReciprocalLatticeVectorBase(const Lattice<dim> &lat) :
    /* init */ BaseType(VectorDimI::Zero()),
    /* init */ lattice(lat)
    ///* base init */ BaseType(LatticeBaseType::d2contra(d))
    { /*!@param[in] d vector in real space
          * Constructs *this by mapping d to the lattice
          */
    }

    template <int dim>
    ReciprocalLatticeVectorBase<dim>::ReciprocalLatticeVectorBase(const VectorDimD &d,
                            /*                   */ const Lattice<dim> &lat) :
    /* init */ BaseType(LatticeCore<dim>::integerCoordinates(d, lat.latticeBasis.transpose())),
    /* init */ lattice(lat)
    { /*!@param[in] d vector in real space
          * Constructs *this by mapping d to the lattice
          */
    }

    template <int dim>
    ReciprocalLatticeVectorBase<dim>& ReciprocalLatticeVectorBase<dim>::operator=(const ReciprocalLatticeVectorBase<dim> &other)
    {
        assert(&lattice == &other.lattice && "ReciprocalLatticeVectorBase<dim> belong to different Lattices.");
        base() = other.base();
        return *this;
    }

    template <int dim>
    ReciprocalLatticeVectorBase<dim>& ReciprocalLatticeVectorBase<dim>::operator=(ReciprocalLatticeVectorBase<dim> &&other)
    {
        assert(&lattice == &other.lattice && "LatticeVectors belong to different Lattices.");
        base() = other.base();
        return *this;
    }

    template <int dim>
    ReciprocalLatticeVectorBase<dim> ReciprocalLatticeVectorBase<dim>::operator+(const ReciprocalLatticeVectorBase<dim> &other) const
    {
        assert(&lattice == &other.lattice && "LatticeVectors belong to different Lattices.");
        return ReciprocalLatticeVectorBase<dim>(static_cast<VectorDimI>(*this) + static_cast<VectorDimI>(other), lattice);
    }

    template <int dim>
    ReciprocalLatticeVectorBase<dim>& ReciprocalLatticeVectorBase<dim>::operator+=(const ReciprocalLatticeVectorBase<dim> &other)
    {
        assert(&lattice == &other.lattice && "LatticeVectors belong to different Lattices.");
        base() += other.base();
        return *this;
    }

    template <int dim>
    ReciprocalLatticeVectorBase<dim> ReciprocalLatticeVectorBase<dim>::operator-(const ReciprocalLatticeVectorBase<dim> &other) const
    {
        assert(&lattice == &other.lattice && "LatticeVectors belong to different Lattices.");
        return ReciprocalLatticeVectorBase<dim>(static_cast<VectorDimI>(*this) - static_cast<VectorDimI>(other), lattice);
    }

    template <int dim>
    ReciprocalLatticeVectorBase<dim>& ReciprocalLatticeVectorBase<dim>::operator-=(const ReciprocalLatticeVectorBase<dim> &other)
    {
        assert(&lattice == &other.lattice && "LatticeVectors belong to different Lattices.");
        base() -= other.base();
        return *this;
    }

    template <int dim>
    ReciprocalLatticeVectorBase<dim> ReciprocalLatticeVectorBase<dim>::operator*(const IntScalarType& scalar) const
    {
        return ReciprocalLatticeVectorBase<dim>(static_cast<VectorDimI>(*this) * scalar, lattice);
    }

    template <int dim>
    typename ReciprocalLatticeVectorBase<dim>::IntScalarType ReciprocalLatticeVectorBase<dim>::dot(const LatticeVector<dim> &other) const
    {
        assert(&lattice == &other.lattice && "LatticeVectors belong to different Lattices.");
        return static_cast<VectorDimI>(*this).dot(static_cast<VectorDimI>(other));
    }

    

    template <int dim>
    typename ReciprocalLatticeVectorBase<dim>::VectorDimD ReciprocalLatticeVectorBase<dim>::cartesian() const
    {
        return lattice.reciprocalBasis * this->template cast<double>();
    }
    
    template <int dim>
    double ReciprocalLatticeVectorBase<dim>::planeSpacing() const
    {
        return 1.0 / cartesian().norm();
    }

    template <int dim>
    typename ReciprocalLatticeVectorBase<dim>::VectorDimD ReciprocalLatticeVectorBase<dim>::interplaneVector() const
    {
        const VectorDimD c(cartesian());
        return c / c.squaredNorm();
    }

    template <int dim>
    typename ReciprocalLatticeVectorBase<dim>::IntScalarType ReciprocalLatticeVectorBase<dim>::closestPlaneIndexOfPoint(const VectorDimD &P) const
    {
        assert(this->squaredNorm() > 0 && "A null ReciprocalLatticeVector cannot be used to compute planeIndexOfPoint");
        const double hd(cartesian().dot(P));
        return std::lround(hd);
    }

    template <int dim>
    typename ReciprocalLatticeVectorBase<dim>::IntScalarType ReciprocalLatticeVectorBase<dim>::planeIndexOfPoint(const VectorDimD &P) const
    {
        assert(this->squaredNorm() > 0 && "A null ReciprocalLatticeVector cannot be used to compute planeIndexOfPoint");
        const double hd(cartesian().dot(P));
        const IntScalarType h(std::lround(hd));
        if (fabs(hd - h) > FLT_EPSILON)
        {
            std::cout << "P=" << P.transpose() << std::endl;
            std::cout << "r=" << this->cartesian().transpose() << std::endl;
            std::cout << "hd=" << std::setprecision(15) << std::scientific << hd << std::endl;
            std::cout << "h=" << h << std::endl;
            assert(0 && "P in not on a lattice plane.");
        }
        return h;
    }

    template <int dim>
    typename ReciprocalLatticeVectorBase<dim>::IntScalarType ReciprocalLatticeVectorBase<dim>::planeIndexOfPoint(const LatticeVector<dim> &P) const
    {
        assert(this->squaredNorm() > 0 && "A null ReciprocalLatticeVector cannot be used to compute planeIndexOfPoint");
        return dot(P);
    }

//    template <int dim>
//    typename ReciprocalLatticeVectorBase<dim>::VectorDimI ReciprocalLatticeVectorBase<dim>::d2cov(const VectorDimD &d,
//                            const Lattice<dim> &lat)
//    {
//        const VectorDimD nd(lat.latticeBasis.transpose() * d);
//        //            const VectorDimD rd(RoundEigen<double,dim>::round(nd));
//        const VectorDimD rd(nd.array().round());
//        if ((nd - rd).norm() > roundTol)
//        {
//            std::cout << "d2cov, nd=" << nd.transpose() << std::endl;
//            std::cout << "d2cov, rd=" << rd.transpose() << std::endl;
//            assert(0 && "Input vector is not a reciprocal lattice vector");
//        }
//        //            assert((nd-rd).norm()<roundTol && "Input vector is not a lattice vector");
//        return rd.template cast<IntScalarType>();
//    }

    //Operator starts here
    template <int dim>
    ReciprocalLatticeVectorBase<dim> operator*(const typename ReciprocalLatticeVectorBase<dim>::IntScalarType& scalar, const ReciprocalLatticeVectorBase<dim> &L)
    {
        return L * scalar;
    }


template <int dim>
ReciprocalLatticeVector<dim>::ReciprocalLatticeVector(const Lattice<dim> &lat) :
/* init base */ ReciprocalLatticeVectorBase<dim>(lat)
{ /*!@param[in] d vector in real space
      * Constructs *this by mapping d to the lattice
      */
}

template <int dim>
ReciprocalLatticeVector<dim>::ReciprocalLatticeVector(const VectorDimD &d,
                        /*                   */ const Lattice<dim> &lat) :
/* init base */ ReciprocalLatticeVectorBase<dim>(d,lat)
{ /*!@param[in] d vector in real space
      * Constructs *this by mapping d to the lattice
      */
}

ReciprocalLatticeVector<3>::ReciprocalLatticeVector(const Lattice<3> &lat) :
/* init base */ ReciprocalLatticeVectorBase<3>(lat)
{ /*!@param[in] d vector in real space
      * Constructs *this by mapping d to the lattice
      */
}

ReciprocalLatticeVector<3>::ReciprocalLatticeVector(const VectorDimD &d,
                        /*                   */ const Lattice<3> &lat) :
/* init base */ ReciprocalLatticeVectorBase<3>(d,lat)
{ /*!@param[in] d vector in real space
      * Constructs *this by mapping d to the lattice
      */
}

LatticeDirection<3> ReciprocalLatticeVector<3>::cross(const ReciprocalLatticeVectorBase<3> &other) const
{
    assert(&lattice == &other.lattice && "LatticeVectors belong to different Lattices.");
    return LatticeDirection<3>(LatticeVector<3>(static_cast<VectorDimI>(*this).cross(static_cast<VectorDimI>(other)), lattice));
}

template class ReciprocalLatticeVectorBase<1>;
template ReciprocalLatticeVectorBase<1> operator*(const typename ReciprocalLatticeVectorBase<1>::IntScalarType& scalar, const ReciprocalLatticeVectorBase<1> &L);
template class ReciprocalLatticeVectorBase<2>;
template ReciprocalLatticeVectorBase<2> operator*(const typename ReciprocalLatticeVectorBase<2>::IntScalarType&scalar, const ReciprocalLatticeVectorBase<2> &L);
    template class ReciprocalLatticeVectorBase<3>;
    template ReciprocalLatticeVectorBase<3> operator*(const typename ReciprocalLatticeVectorBase<3>::IntScalarType& scalar, const ReciprocalLatticeVectorBase<3> &L);

template class ReciprocalLatticeVector<1>;
template class ReciprocalLatticeVector<2>;

} // end namespace
#endif
