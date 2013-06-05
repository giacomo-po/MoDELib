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

namespace model {
    
    
    class ChargedParticle : public SpatialCellParticle<ChargedParticle,3>
    {
        
    public:
        
        enum{dim=3};
        typedef  SpatialCellParticle<ChargedParticle,3>::PositionType PositionType;
        typedef  SpatialCellParticle<ChargedParticle,3>::PositionType ForceType;
        

        
    private:

        PositionType _p;
        double _q; // the electric charge of the particle
        
    public:
        
        ForceType force;
        double energy;
        
        /*****************************************/
        ChargedParticle(const PositionType& pIN, const double& qIN) :
        /* init list */ SpatialCellParticle<ChargedParticle,3>::SpatialCellParticle(pIN), //  SpatialCellParticle must be constructed with initial position
        /* init list */ _p(pIN),
        /* init list */ _q(qIN),
        /* init list */ force(ForceType::Zero()),
        /* init list */ energy(0.0)
        {/*! Constructor with input position and charge
          */
            //std::cout<< "Creating particle ";
            //std::cout<<*this;
        }
        
        /*****************************************/
        const double& q() const
        {/*! The charge of this ChargedParticle
          */
            return _q;
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
            os  <<cP.sID<<" "
            <<cP.P().transpose()<<" "
            <<cP.q()<<" "
            <<cP.force.transpose()<<" "
            <<cP.energy;
            return os;
        }
            
    };
            
}
#endif
