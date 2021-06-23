/* This file is part of PIL, the Particle Interaction Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LatticeVector_cpp_
#define model_LatticeVector_cpp_

#include<LatticeModule.h>

namespace model
{

    /**********************************************************************/
    template <int dim>
    typename LatticeVector<dim>::BaseType& LatticeVector<dim>::base()
    {
        return *this;
    }

    /**********************************************************************/
    template <int dim>
    const typename LatticeVector<dim>::BaseType& LatticeVector<dim>::base() const
    {
        return *this;
    }


    /**********************************************************************/
    template <int dim>
    LatticeVector<dim>::LatticeVector(const LatticeType &lat) :
    /* init base */ BaseType(VectorDimI::Zero()),
    /* init      */ lattice(lat)
    ///* base init */ BaseType(LatticeBaseType::d2contra(d))
    { /*!@param[in] d vector in real space
          * Constructs *this by mapping d to the lattice
          */
    }

    /**********************************************************************/
    template <int dim>
    LatticeVector<dim>::LatticeVector(const VectorDimD &d,
                  const LatticeType &lat) :
    /* init base */ BaseType(d2contra(d, lat)),
    /* init      */ lattice(lat)
    { /*!@param[in] d vector in real space
          * Constructs *this by mapping d to the lattice
          */
    }

//    /**********************************************************************/
//    template <int dim>
//    template <typename OtherDerived>
//    LatticeVector<dim>::LatticeVector(const Eigen::MatrixBase<OtherDerived> &other,
//                  const LatticeType &lat) :
//    /* init base */ BaseType(other),
//    /* init      */ lattice(lat)
//    { /*!@param[in] d vector in real space
//          * Constructs *this by mapping d to the lattice
//          */
//    }


    /**********************************************************************/
    template <int dim>
    typename LatticeVector<dim>::LatticeVectorType& LatticeVector<dim>::operator=(const LatticeVectorType &other)
    {
        assert(&lattice == &other.lattice && "LatticeVectors belong to different Lattices.");
        base() = other.base();
        return *this;
    }

    /**********************************************************************/
    template <int dim>
    typename LatticeVector<dim>::LatticeVectorType& LatticeVector<dim>::operator=(LatticeVectorType &&other)
    {
        assert(&lattice == &other.lattice && "LatticeVectors belong to different Lattices.");
        base() = other.base();
        return *this;
    }

    /**********************************************************************/
    template <int dim>
    typename LatticeVector<dim>::LatticeVectorType LatticeVector<dim>::operator+(const LatticeVectorType &other) const
    {
        assert(&lattice == &other.lattice && "LatticeVectors belong to different Lattices.");
        return LatticeVectorType(static_cast<VectorDimI>(*this) + static_cast<VectorDimI>(other), lattice);
    }

    /**********************************************************************/
    template <int dim>
    typename LatticeVector<dim>::LatticeVectorType& LatticeVector<dim>::operator+=(const LatticeVectorType &other)
    {
        assert(&lattice == &other.lattice && "LatticeVectors belong to different Lattices.");
        base() += other.base();
        return *this;
    }

    /**********************************************************************/
    template <int dim>
    typename LatticeVector<dim>::LatticeVectorType LatticeVector<dim>::operator-(const LatticeVectorType &other) const
    {
        assert(&lattice == &other.lattice && "LatticeVectors belong to different Lattices.");
        return LatticeVectorType(static_cast<VectorDimI>(*this) - static_cast<VectorDimI>(other), lattice);
    }

    /**********************************************************************/
    template <int dim>
    typename LatticeVector<dim>::LatticeVectorType& LatticeVector<dim>::operator-=(const LatticeVectorType &other)
    {
        assert(&lattice == &other.lattice && "LatticeVectors belong to different Lattices.");
        base() -= other.base();
        return *this;
    }

//    /**********************************************************************/
//    template <int dim>
//    typename LatticeVector<dim>::LatticeVectorType LatticeVector<dim>::operator*(const long int &scalar) const
//    {
//        return LatticeVectorType(static_cast<VectorDimI>(*this) * scalar, lattice);
//    }

    /**********************************************************************/
    template <int dim>
    long int LatticeVector<dim>::dot(const ReciprocalLatticeVectorType &other) const
    {
        assert(&lattice == &other.lattice && "LatticeVectors belong to different Lattices.");
        return static_cast<VectorDimI>(*this).dot(static_cast<VectorDimI>(other));
    }

    //        /**********************************************************************/
    //        ReciprocalLatticeVectorType cross(const LatticeVectorType& other) const
    //        {
    //            assert(&lattice==&other.lattice && "LatticeVectors belong to different Lattices.");
    //            return ReciprocalLatticeVectorType(static_cast<VectorDimI>(*this).cross(static_cast<VectorDimI>(other)),lattice);
    //        }

    /**********************************************************************/
    template <int dim>
    typename LatticeVector<dim>::ReciprocalLatticeDirectionType LatticeVector<dim>::cross(const LatticeVectorType &other) const
    {
        assert(&lattice == &other.lattice && "LatticeVectors belong to different Lattices.");
        return ReciprocalLatticeDirectionType(ReciprocalLatticeVectorType(static_cast<VectorDimI>(*this).cross(static_cast<VectorDimI>(other)), lattice));
    }

    /**********************************************************************/
    template <int dim>
    typename LatticeVector<dim>::VectorDimD LatticeVector<dim>::cartesian() const
    {
        return lattice.latticeBasis * this->template cast<double>();
    }

    /**********************************************************************/
    template <int dim>
    typename LatticeVector<dim>::VectorDimI LatticeVector<dim>::d2contra(const VectorDimD &d,
                               const LatticeType &lat)
    {
        const VectorDimD nd(lat.reciprocalBasis.transpose() * d);
        //            const VectorDimD rd(RoundEigen<double,dim>::round(nd));
        const VectorDimD rd(nd.array().round());
        if ((nd - rd).norm() > roundTol)
        {
            std::cout << "d2contra, nd=" << nd.transpose() << std::endl;
            std::cout << "d2contra, rd=" << rd.transpose() << std::endl;
            assert(0 && "Input vector is not a lattice vector");
        }
        return rd.template cast<long int>();
    }
    
    template<int dim>
    LatticeVector<dim> operator*(const long int& scalar, const LatticeVector<dim>& L)
    {
        return L*scalar;
    }

    template class LatticeVector<3>;
    template LatticeVector<3>operator*(const long int& scalar, const LatticeVector<3>& L);

} // end namespace
#endif

