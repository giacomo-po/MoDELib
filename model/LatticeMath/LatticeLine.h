/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LatticeLine_h_
#define model_LatticeLine_h_

#include <model/LatticeMath/LatticeVector.h>
#include <model/LatticeMath/LatticeDirection.h>
#include <math.h>       /* round, floor, ceil, trunc */

namespace model
{
    struct LatticeLine
    {
        

        
        typedef LatticeVector<3>    LatticeVectorType;
        typedef LatticeDirection<3> LatticeDirectionType;
        typedef typename LatticeVectorType::VectorDimD VectorDimD;

        
        const LatticeVectorType P;
        const LatticeDirectionType d;
        
        LatticeLine(const LatticeVectorType& P_in,const LatticeVectorType& d_in) :
        /* init */ P(P_in),
        /* init */ d(d_in)
        {}
        
        
        LatticeLine(const LatticeVectorType& P_in,const LatticeDirectionType& d_in) :
        /* init */ P(P_in),
        /* init */ d(d_in)
        {}
        
        /**********************************************************************/
        VectorDimD snapToLattice(const VectorDimD& P0) const
        {
            const VectorDimD dc(d.cartesian());
            const VectorDimD Pc(P.cartesian());
            const double n=(P0-Pc).dot(dc)/dc.squaredNorm();
            return Pc+round(n)*dc;
        }
        
//        /**********************************************************************/
//        VectorDimD snapToDirection(const VectorDimD& DP) const
//        {
//            const VectorDimD dc(d.cartesian());
////            const double n=;
//            return round(DP.dot(dc)/dc.squaredNorm())*dc;
//        }
        
        /**********************************************************************/
        bool contains(const LatticeVectorType& P0) const
        {
            return LatticeDirectionType(LatticeVectorType(P0-P)).cross(d).squaredNorm()==0;
        }
        
    };
    
} // end namespace
#endif

//        /**********************************************************************/
//        PlaneLineIntersection intersectWith(const LatticePlane& plane) const
//        {
//            return PlaneLineIntersection(plane,*this);
//        }


//        /**********************************************************************/
//        LatticeVectorType closestPoint(const VectorDimD& P)
//        {
////            const double u();
//
//        }
