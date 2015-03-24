/* This file is part of PIL, the Particle Interaction Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LatticeLine_h_
#define model_LatticeLine_h_

#include <model/DislocationDynamics/LatticeMath/LatticeVector.h>
#include <model/DislocationDynamics/LatticeMath/LatticeDirection.h>
#include <model/DislocationDynamics/LatticeMath/LatticePlane.h>

namespace model
{
    struct LatticeLine
    {
        

        
        typedef LatticeVector<3>    LatticeVectorType;
        typedef LatticeDirection<3> LatticeDirectionType;
        typedef typename LatticeVectorType::VectorDimD VectorDimD;

        
        LatticeVectorType P;
        LatticeDirectionType d;
        
        LatticeLine(const LatticeVectorType& P_in,const LatticeVectorType d_in) :
        /* init */ P(P_in),
        /* init */ d(d_in)
        {}
        
        
        LatticeLine(const LatticeVectorType& P_in,const LatticeDirectionType d_in) :
        /* init */ P(P_in),
        /* init */ d(d_in)
        {}
        
        /**********************************************************************/
        std::pair<bool,LatticeVectorType> intersectWith(const LatticePlane& plane) const
        {
            std::pair<bool,LatticeVectorType> temp(false,LatticeVectorType());
            
            const long int num=plane.P.dot(plane.n)-P.dot(plane.n);
            
            if (num==0)
            {
                temp.first=true;
                temp.second=P;
            }
            else
            {
                const long int den=d.dot(plane.n);
                if(den!=0)
                {
                    const long int ri=num/den;
                    if(std::fabs(double(num)/den-ri)<LatticeVectorType::roundTol) // ri is a true integer
                    {
                        temp.first=true;
//                        temp.second=P+ri*d;
                        assert(0 && "FINISH HERE");
                    }
                }
            }
            return temp;
        }
        
        
        /**********************************************************************/
        LatticeVectorType closestPoint(const VectorDimD& P)
        {
//            const double u();
            
        }
        
        
    };
    
} // end namespace
#endif
