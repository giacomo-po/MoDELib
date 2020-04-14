/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LatticePlane_h_
#define model_LatticePlane_h_

#include <iomanip>
#include <LatticeVector.h>
#include <LatticePlaneBase.h>
#include <Plane.h>

namespace model
{
    
    
    template<int dim>
    struct LatticePlaneKey
    {// dim ints for direction, 1 int for heigth
        
        
        typedef long int LongIntType;
        typedef Eigen::Matrix<double,dim,1> VectorDimD;
        typedef Eigen::Matrix<LongIntType,dim,1> VectorDimI;
        typedef std::array<LongIntType,dim+2> ArrayType;
        //
        /**************************************************************************/
        static ArrayType hr2array(const VectorDimI& r, const LongIntType& h,const LongIntType& latticeID)
        {
            ArrayType temp;
            for(int k=0;k<dim;++k)
            {
                temp[k]=r(k);
            }
            temp[dim]=h;
            temp[dim+1]=latticeID;
            return temp;
        }
        
        /**************************************************************************/
        static int sgn(const LongIntType& val)
        {
            return (LongIntType(0) < val) - (val < LongIntType(0));
        }
        
        /**************************************************************************/
        static ArrayType correct_h_sign(VectorDimI r,LongIntType h,const LongIntType& latticeID)
        {
            assert(r.squaredNorm()>0 && "zero direction cannot be used as LatticePlaneKey");
            for(int d=0;d<dim;++d)
            {
                const int sgnrd(sgn(r(d)));
                if(sgnrd!=0)
                {
                    r*=sgnrd;
                    h*=sgnrd;
                    break;
                }
            }
            return hr2array(r,h,latticeID);
        }
        
        
        const ArrayType array;
        
        /**********************************************************************/
        LatticePlaneKey(const VectorDimI& r,
                        const LongIntType& h,
                        const LongIntType& latticeID) :
        /* init */ array(correct_h_sign(r,h,latticeID))
        {
        }
        
        /**********************************************************************/
        LatticePlaneKey(const VectorDimD& P,
                        const ReciprocalLatticeDirection<dim>& r) :
        /* delegate */ LatticePlaneKey(r,r.planeIndexOfPoint(P),r.lattice.sID)
        {
            
        }
        
        /**********************************************************************/
        LatticePlaneKey(const LatticeVector<dim>& L,
                        const ReciprocalLatticeDirection<dim>& r) :
        /* delegate */ LatticePlaneKey(r,r.planeIndexOfPoint(L),r.lattice.sID)
        {
            
        }
        
        /**********************************************************************/
        LatticePlaneKey(const long int& hin,
                        const ReciprocalLatticeDirection<dim>& r) :
        /* init */ LatticePlaneKey(r,hin,r.lattice.sID)
        {
            
        }

        
        Eigen::Map<const VectorDimI> reciprocalDirectionComponents() const
        {
            return Eigen::Map<const VectorDimI>(&array.data()[0]);
        }
        
        const LongIntType& planeIndex() const
        {
            return array[dim];
        }

        const LongIntType& latticeID() const
        {
            return array[dim+1];
        }
        
        bool operator<(const LatticePlaneKey<dim>& other) const
        {// allow use in sorted STL containers
            return array<other.array;
        }
        
        /**********************************************************************/
        template <class T>
        friend T& operator << (T& os, const LatticePlaneKey<dim>& key)
        {
            for(const auto& val : key.array)
            {
                os<<val<<" ";
            }
            return os;
        }
        
    };
    
    
    struct LatticePlane
    {
        static constexpr int dim=3;
        typedef Eigen::Matrix<double,dim,1> VectorDimD;
        typedef LatticeVector<dim>    LatticeVectorType;
        
        /**********************************************************************/
        static std::pair<bool,long int> computeHeight(const ReciprocalLatticeDirection<dim>& r,
                                                      const VectorDimD& P)
        {/*! ????
          */
            assert(r.squaredNorm()>0 && "A zero normal cannot be used as valid GlidePlane normal");
            const double hd(P.dot(r.cartesian()));
            const long int h(std::lround(hd));            
            return  std::make_pair(fabs(hd-h)<FLT_EPSILON,h);
        }
        
        /**********************************************************************/
        static std::pair<bool,long int> computeHeight(const ReciprocalLatticeDirection<dim>& r,
                                                      const LatticeVector<dim>& L)
        {/*! ????
          */
            assert(r.squaredNorm()>0 && "A zero normal cannot be used as valid GlidePlane normal");
            return  std::make_pair(true,L.dot(r));
        }
        
        const LatticePlaneKey<dim> key;
        const ReciprocalLatticeDirection<dim> n;
        const LatticePlaneKey<dim>::LongIntType planeIndex;
        
        /**********************************************************************/
        LatticePlane(const VectorDimD& P,
                     const ReciprocalLatticeDirection<dim>& r) :
        /* init */ key(P,r)
        /* init */,n(ReciprocalLatticeVector<dim>(key.reciprocalDirectionComponents(),r.lattice))
        /* init */,planeIndex(key.planeIndex())
        {

        }
        
        /**********************************************************************/
        LatticePlane(const LatticeVector<dim>& L,
                     const ReciprocalLatticeDirection<dim>& r) :
        /* init */ key(L,r)
        /* init */,n(ReciprocalLatticeVector<dim>(key.reciprocalDirectionComponents(),r.lattice))
        /* init */,planeIndex(key.planeIndex())
        {
            
        }
        
        /**********************************************************************/
        LatticePlane(const long int& hin,
                     const ReciprocalLatticeDirection<dim>& r) :
        /* init */ key(hin,r)
        /* init */,n(ReciprocalLatticeVector<dim>(key.reciprocalDirectionComponents(),r.lattice))
        /* init */,planeIndex(key.planeIndex())
        {
            
        }
        
        /**********************************************************************/
        VectorDimD planeOrigin() const
        {
            return planeIndex*n.planeSpacing()*n.cartesian().normalized();
        }
        

        
    };
    
}
#endif

