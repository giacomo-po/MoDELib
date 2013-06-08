/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2012 by Giacomo Po <gpo@ucla.edu>
 * Copyright (C) 2012 by Shao-Ching Huang <sch@ucla.edu>
 * Copyright (C) 2012 by Tajendra Singh <tvsingh@ucla.edu>
 * Copyright (C) 2012 by Tamer Crosby <tcrosby@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _ChargedParticle_h
#define _ChargedParticle_h

#include <tutorials/ParticleInteraction/ChargedParticles/ChargedCell.h>
#include <model/SpaceDecomposition/SpatialCellParticle.h>
#include <tutorials/ParticleInteraction/ChargedParticles/CoulombForce.h>

namespace model {
    
    
    
    
    class ChargedParticle :
    /* inheritance     */ public SpatialCellParticle<ChargedParticle,3>
    {
        
    public:
        
        enum{dim=3};
        typedef  SpatialCellParticle<ChargedParticle,3>::PositionType PositionType;
        typedef  SpatialCellParticle<ChargedParticle,3>::PositionType ForceType;
        
        typedef CoulombForce<ChargedParticle> CoulombForceInteraction;
        
    private:

        PositionType _p; // THIS SHOULD BE STORED IN SpatialCellParticle

    public:
        
        ForceType force;
        const double q; // the electric charge of the particle
        double energy;
        
        /*****************************************/
        ChargedParticle(const PositionType& pIN, const double& qIN) :
        /* init list */ SpatialCellParticle<ChargedParticle,3>::SpatialCellParticle(pIN), //  SpatialCellParticle must be constructed with initial position
        /* init list */ _p(pIN),
        /* init list */ force(ForceType::Zero()),
        /* init list */ q(qIN),
        /* init list */ energy(0.0)
        {/*!@param[in] pIN position of this ChargedParticle
          * @param[in] qIN charge of this ChargedParticle
          * 
          * Constructor with input position and charge
          */
        }
                
        /*****************************************/
        const PositionType& P() const
        {/*! The charge of this ChargedParticle
          */
            return _p;
        }
        
        /*****************************************/
        template <class T>
        friend T& operator << (T& os, const ChargedParticle& cP)
        {/*! operator << is used to output ChargedParticle info
          *  Example:
          *  ChargedParticle p;
          *  std::cout<<p;
          */
            os<<cP.sID<<"\t"
            <<cP.P().transpose()<<"\t"
            <<cP.q<<"\t"
//            <<cP.force.transpose()<<"\t"
            <<cP.get<CoulombForceInteraction>().transpose()<<"\t"
            <<cP.energy<<"\t"
            <<cP.mpiID<<"\t";
            return os;
        }
            
    };
            
}
#endif
