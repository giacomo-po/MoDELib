/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SpatialCellProperty_h_
#define model_SpatialCellProperty_h_

#include <model/Utilities/CRTP.h>
#include <model/Utilities/StaticID.h>


namespace model {
	
    
    
	/********************************************************************************************/
	/********************************************************************************************/
	/*! \brief A dim-dimensional cell occupying the spatial region cellID<= x/cellSize < (cellID+1).
	 *  SpatialCell is aware off all ParticleType objects present inside it.
	 */
	template<typename Derived>
	struct SpatialCellProperty :
    /*                      */ private  CRTP<Derived>
    
    {
        
        
//        const Derived& get()
//        {
//            return this->derived();
//        }
        
	};
    
    
	
	////////////////////////////////////////////////////////////////////////////////
}	// close namespace model
#endif

