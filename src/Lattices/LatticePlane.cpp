/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LatticePlane_cpp_
#define model_LatticePlane_cpp_

#include <LatticeModule.h>
namespace model
{


   /**************************************************************************/
    template <int dim>
    typename LatticePlaneKey<dim>::ArrayType LatticePlaneKey<dim>::hr2array(const VectorDimI &r, const LongIntType &h, const LongIntType &latticeID)
    {
        ArrayType temp;
        for (int k = 0; k < dim; ++k)
        {
            temp[k] = r(k);
        }
        temp[dim] = h;
        temp[dim + 1] = latticeID;
        return temp;
    }

    /**************************************************************************/
    template <int dim>
    int LatticePlaneKey<dim>::sgn(const LongIntType &val)
    {
        return (LongIntType(0) < val) - (val < LongIntType(0));
    }

    /**************************************************************************/
    template <int dim>
    typename LatticePlaneKey<dim>::ArrayType LatticePlaneKey<dim>::correct_h_sign(VectorDimI r, LongIntType h, const LongIntType &latticeID)
    {
        if(r.squaredNorm()==0)
        {
            throw std::runtime_error("LatticePlaneKey::correct_h_sign: A zero normal cannot be used as valid GlidePlane normal.");
        }
        for (int d = 0; d < dim; ++d)
        {
            const int sgnrd(sgn(r(d)));
            if (sgnrd != 0)
            {
                r *= sgnrd;
                h *= sgnrd;
                break;
            }
        }
        return hr2array(r, h, latticeID);
    }


    /**********************************************************************/
    template <int dim>
    LatticePlaneKey<dim>::LatticePlaneKey(const VectorDimI& r,
    const LongIntType& h,
    const LongIntType& latticeID) :
    /* init */ array(correct_h_sign(r,h,latticeID))
    {
    }
    
    /**********************************************************************/
    template <int dim>
    LatticePlaneKey<dim>::LatticePlaneKey(const VectorDimD& P,
    const ReciprocalLatticeDirection<dim>& r) :
    /* delegate */ LatticePlaneKey(r,r.planeIndexOfPoint(P),r.lattice.sID)
    {
    
    }
    
    /**********************************************************************/
    template <int dim>
    LatticePlaneKey<dim>::LatticePlaneKey(const LatticeVector<dim>& L,
    const ReciprocalLatticeDirection<dim>& r) :
    /* delegate */ LatticePlaneKey(r,r.planeIndexOfPoint(L),r.lattice.sID)
    {
    
    }
    
    /**********************************************************************/
    template <int dim>
    LatticePlaneKey<dim>::LatticePlaneKey(const long int& hin,
    const ReciprocalLatticeDirection<dim>& r) :
    /* init */ LatticePlaneKey(r,hin,r.lattice.sID)
    {
    
    }

    template <int dim>
    Eigen::Map<const typename LatticePlaneKey<dim>::VectorDimI> LatticePlaneKey<dim>::reciprocalDirectionComponents() const
    {
        return Eigen::Map<const VectorDimI>(&array.data()[0]);
    }

    template <int dim>
    const typename LatticePlaneKey<dim>::LongIntType& LatticePlaneKey<dim>::planeIndex() const
    {
        return array[dim];
    }

    template <int dim>
    const typename LatticePlaneKey<dim>::LongIntType& LatticePlaneKey<dim>::latticeID() const
    {
        return array[dim + 1];
    }

    template <int dim>
    bool LatticePlaneKey<dim>::operator<(const LatticePlaneKey<dim> &other) const
    { // allow use in sorted STL containers
        return array < other.array;
    }

    template struct LatticePlaneKey<3>;

//for the lattice plane class.
    /**********************************************************************/
    std::pair<bool, long int> LatticePlane::computeHeight(const ReciprocalLatticeDirection<dim> &r,
                                                   const VectorDimD &P)
    { /*! ????
          */
        if(r.squaredNorm()==0)
        {
            throw std::runtime_error("LatticePlane::computeHeight: A zero normal cannot be used as valid GlidePlane normal");
        }
        const double hd(P.dot(r.cartesian()));
        const long int h(std::lround(hd));
        return std::make_pair(fabs(hd - h) < FLT_EPSILON, h);
    }

    /**********************************************************************/
    std::pair<bool, long int> LatticePlane::computeHeight(const ReciprocalLatticeDirection<dim> &r,
                                                   const LatticeVector<dim> &L)
    { /*! ????
          */
        if(r.squaredNorm()==0)
        {
            throw std::runtime_error("LatticePlane::computeHeight: A zero normal cannot be used as valid GlidePlane normal");
        }
        return std::make_pair(true, L.dot(r));
    }

    /**********************************************************************/
    LatticePlane::LatticePlane(const VectorDimD& P,
    const ReciprocalLatticeDirection<dim>& r) :
    /* init */ key(P,r)
    /* init */,n(ReciprocalLatticeVector<dim>(key.reciprocalDirectionComponents(),r.lattice))
    /* init */,planeIndex(key.planeIndex())
    {
}
    
    /**********************************************************************/
    LatticePlane::LatticePlane(const LatticeVector<dim>& L,
    const ReciprocalLatticeDirection<dim>& r) :
    /* init */ key(L,r)
    /* init */,n(ReciprocalLatticeVector<dim>(key.reciprocalDirectionComponents(),r.lattice))
    /* init */,planeIndex(key.planeIndex())
    {
    
    }
    
    /**********************************************************************/
    LatticePlane::LatticePlane(const long int& hin,
    const ReciprocalLatticeDirection<dim>& r) :
    /* init */ key(hin,r)
    /* init */,n(ReciprocalLatticeVector<dim>(key.reciprocalDirectionComponents(),r.lattice))
    /* init */,planeIndex(key.planeIndex())
    {
    
    }
    /**********************************************************************/
    typename LatticePlane::VectorDimD LatticePlane::planeOrigin() const
    {
        return planeIndex * n.planeSpacing() * n.cartesian().normalized();
    }
}
#endif

