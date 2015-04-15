/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LatticePlaneBase_h_
#define model_LatticePlaneBase_h_

#include <model/LatticeMath/LatticeVector.h>
#include <model/LatticeMath/LatticeDirection.h>
#include <model/LatticeMath/ReciprocalLatticeDirection.h>

namespace model
{
    struct LatticePlaneBase : public ReciprocalLatticeDirection<3>
    {
        
        typedef LatticeVector<3>    LatticeVectorType;
        typedef LatticeDirection<3>    LatticeDirectionType;
        typedef ReciprocalLatticeDirection<3> ReciprocalLatticeDirectionType;
        
        const LatticeDirectionType d1;
        const LatticeDirectionType d2;
        
        LatticePlaneBase(const LatticeVectorType& v1_in, const LatticeVectorType& v2_in) :
        /* init base */ ReciprocalLatticeDirectionType(v1_in,v2_in),
        /* init */ d1(v1_in),
        /* init */ d2(v2_in)
        {
            assert(d1.squaredNorm()>0.0);
            assert(d2.squaredNorm()>0.0);
            assert(d1.cross(d2).squaredNorm()>0.0);
        }
        
        Eigen::Matrix<double,3,1> snapToLattice(const Eigen::Matrix<double,3,1>& P) const
        {
            const Eigen::Matrix<double,3,2> B((Eigen::Matrix<double,2,3>()<<d1.cartesian().transpose(),d2.cartesian().transpose()).finished().transpose());
            const Eigen::Matrix<double,2,1> nd=(B.transpose()*B).inverse()*B.transpose()*P;
            const Eigen::Matrix<double,3,1>  p=B*RoundEigen<double,2>::round(nd);
            return p;
        }

        
    };
    
} // end namespace
#endif
