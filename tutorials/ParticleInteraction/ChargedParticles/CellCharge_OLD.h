/* This file is part of PIL, the Particle Interaction Library.
 *
 * Copyright (C) 2012 by Giacomo Po <gpo@ucla.edu>
 * Copyright (C) 2012 by Shao-Ching Huang <sch@ucla.edu>
 * Copyright (C) 2012 by Tajendra Singh <tvsingh@ucla.edu>
 * Copyright (C) 2012 by Tamer Crosby <tcrosby@ucla.edu>
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_CellCharge_h
#define model_CellCharge_h

#include <iostream>
#include <model/SpaceDecomposition/SpatialCell.h>
//#include <tutorials/ParticleInteraction/ChargedParticles/ChargedParticle.h>



namespace model {
    
    template <typename ChargedParticleType>
    class CellCharge
    {
    
        typedef SpatialCell<ChargedParticleType,ChargedParticleType::dim> SpatialCellType;
 
    public:
        CellCharge(const SpatialCellType& sC)
        {
            double charge(0.0);
            
            for(typename SpatialCellType::ParticleContainerType::const_iterator pIter=sC.particleBegin();pIter!=sC.particleEnd();++pIter) // loop over neighbor particles
            {
                charge+=(*pIter)->q;
            }
            std::cout<<"Cell "<<sC.cellID.transpose()<<" has total charge: "<<charge<<std::endl;
        }
        
    };
    
}
#endif

