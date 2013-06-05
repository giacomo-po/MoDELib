/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _model_DislocationStress_h_
#define _model_DislocationStress_h_


#ifdef _MODEL_DD_MPI_
#include <model/ParticleInteraction/InteractionBase.h>
#endif

namespace model
{
	
	/**************************************************************************/
	/**************************************************************************/
	template<short unsigned int dim>
	class DislocationStress
#ifdef _MODEL_DD_MPI_
    /* inheritance */ : public InteractionBase<double,dim*(dim-1)>
#endif

    {
        
        
    public:
        
        DislocationStress()
        {/*!
          */
#ifdef _MODEL_DD_MPI_
            this->resultVector[cp1.rID*3+0]+=f(0);
            this->resultVector[cp1.rID*3+1]+=f(1);
            this->resultVector[cp1.rID*3+2]+=f(2);
#else

#endif
        }
				
	};
    
        
    /**************************************************************************/
    /**************************************************************************/
}	// close namespace
#endif
