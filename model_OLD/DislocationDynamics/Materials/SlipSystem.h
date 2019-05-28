/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_SLIPSYSTEM_H_
#define model_SLIPSYSTEM_H_

#include <memory>
#include <assert.h>
#include <LatticePlaneBase.h>
#include <LatticeVector.h>
#include <DislocationMobilityBase.h>

namespace model
{
    
    struct SlipSystem
    {
        

        const LatticePlaneBase n;
        const LatticeVector<3>  s;
        const Eigen::Matrix<double,3,1>  unitNormal;
        const std::shared_ptr<DislocationMobilityBase> mobility;
        
        
        SlipSystem(const LatticeVector<3>& a1,
                   const LatticeVector<3>& a2,
                   const LatticeVector<3>& slip_in,
                   const std::shared_ptr<DislocationMobilityBase>& mobility_in):
        /* init */ n(a1,a2)
        /* init */,s(slip_in)
        /* init */,unitNormal(n.cartesian().normalized())
        /* init */,mobility(mobility_in)
        {
            if(n.dot(s)!=0)
            {
                model::cout<<"PLANE NORMAL AND SLIP DIRECTION ARE NOT ORTHOGONAL. EXITING."<<std::endl;
                exit(EXIT_FAILURE);
            }
            if(!mobility)
            {
                model::cout<<"MOBILITY CANNOT BE A NULLPTR. EXITING."<<std::endl;
                exit(EXIT_FAILURE);
            }
        }
        
    };

}
#endif
