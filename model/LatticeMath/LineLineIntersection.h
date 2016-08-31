/* This file is part of PIL, the Particle Interaction Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LineLineIntersection_h_
#define model_LineLineIntersection_h_

#include <utility>      // std::pair, std::make_pair
#include <model/LatticeMath/LatticeVector.h>
#include <model/LatticeMath/LatticePlane.h>
#include <model/LatticeMath/LatticeLine.h>

namespace model
{
    
    
    struct LineLineIntersection
    {
        
        enum IntersectionType
        {
            parallelLines=0,
            coincidentLines=1,
            nonIntersectingLines=2,
            intersectingLines=3,
            offLattice=4
        };

        
        typedef Eigen::Matrix<long int,3,1> VectorDimI;
        typedef LatticeVector<3> LatticeVectorType;
        typedef LatticeDirection<3> LatticeDirectionType;
        typedef ReciprocalLatticeDirection<3> ReciprocalLatticeDirectionType;

        const LatticeVectorType P12; // d1xd2
        const ReciprocalLatticeDirectionType d1xd2; // d1xd2
        const ReciprocalLatticeDirectionType dPxd1; // (P1-P2)xd1
        const ReciprocalLatticeDirectionType dPxd2; // (P1-P2)xd2
        const IntersectionType intersectionType;
        const LatticeVector<3> P;
        
        
        static IntersectionType getIntersectionType(const LatticeVectorType& p12,
                                       const ReciprocalLatticeDirectionType& d1Xd2,
                                       const ReciprocalLatticeDirectionType& dPXd1,
                                       const ReciprocalLatticeDirectionType& dPXd2)
        {
            IntersectionType temp=parallelLines;
            if(d1Xd2.squaredNorm()==0) // parallel or coincident lines
            {

                if(dPXd1.squaredNorm()==0)
                {
                    assert(dPXd2.squaredNorm()==0 && "THIS SHOULD NEVER HAPPEN");
                    temp=coincidentLines;
                }
                else
                {
//                    temp.first=false;
                    temp=parallelLines;
                }
            }
            else // non-parallel, non-coincident
            {
                if(p12.dot(d1Xd2)!=0) // planarity condition fails
                {
//                    temp.first=false;
                    temp=nonIntersectingLines;
                }
                else // planarity condition ok
                {
                    if(LatticeGCD<3>::gcd(dPXd2.gCD,d1Xd2.gCD)==d1Xd2.gCD && LatticeGCD<3>::gcd(dPXd1.gCD,d1Xd2.gCD)==d1Xd2.gCD)
                    {
//                        temp.first=true;
                        temp=intersectingLines;
                    }
                    else
                    {
//                        temp.first=false;
                        temp=offLattice;
                    }
                }
            }
            return temp;
        }
        
        
        
        LineLineIntersection(const LatticeLine& l1, const LatticeLine& l2) :
        /* init */ P12(l1.P-l2.P),
        /* init */ d1xd2(l1.d,l2.d), // d1xd2
        /* init */ dPxd1(P12,l1.d), // (P1-P2)xd1
        /* init */ dPxd2(P12,l2.d), // (P1-P2)xd2
        /* init */ intersectionType(getIntersectionType(P12,d1xd2,dPxd1,dPxd2)),
        /* init */ P(intersectionType==coincidentLines? l1.P : (intersectionType==intersectingLines? (l1.P-dPxd2.gCD/d1xd2.gCD*(1-2*((dPxd2+d1xd2).squaredNorm()==0))*l1.d).eval() : VectorDimI::Zero()))
        {

        }
        
        friend std::ostream& operator << (std::ostream& os, const LineLineIntersection& lli)
        {
            os  <<lli.intersectionType<<"\n"
            /**/<<lli.P.cartesian().transpose();
            return os;
        }
        
    };
    
} // end namespace
#endif
