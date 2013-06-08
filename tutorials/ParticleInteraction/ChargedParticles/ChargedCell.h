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

#ifndef _ChargedCell_h
#define _ChargedCell_h

#include <Eigen/Dense>
#include <model/Utilities/TypeTraits.h>
#include <model/SpaceDecomposition/SpatialCell.h>


/**************************************************************************/
/* TRAITS *****************************************************************/
namespace model {
    
    class ChargedParticle;
    
    struct ChargedCell;
    
    template<>
    struct TypeTraits<ChargedParticle> {
        //enum {dim=3;}
        typedef ChargedParticle ParticleType;
        typedef ChargedCell     CellType;
    };
    
    template<>
    struct TypeTraits<ChargedCell> {
        typedef ChargedParticle ParticleType;
        typedef ChargedCell     CellType;
    };
    
    
    struct ChargedCell : public SpatialCell<ChargedCell,3>
    {
        
        typedef typename TypeTraits<ChargedCell>::ParticleType ParticleType;
        typedef SpatialCell<ChargedCell,3>::CellIdType CellIdType;
        
        
        /* Constructor *******************************************/
        ChargedCell(const CellIdType& cellID_in) :
        /* Base constructor */ SpatialCell<ChargedCell,3>::SpatialCell(cellID_in)
        {
        }
        
        
        /*****************************************/
        template <class T>
        friend T& operator<< (T& os, const ChargedCell& cC)
        {/*! Operator << use ParticleType specific << operator
          */
            os<<cC.cellID.transpose()<<"\t"<<cC.size()<<"\t"<<cC.neighborSize()<<"\t"<<cC.assignedRank;
            return os;
        }
        
    };
    
}
#endif
