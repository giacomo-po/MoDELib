/* This file is part of PIL, the Particle Interaction Library.
 *
 * Copyright (C) 2012 by Giacomo Po <gpo@ucla.edu>
 * Copyright (C) 2012 by Shao-Ching Huang <sch@ucla.edu>
 * Copyright (C) 2012 by Tajendra Singh <tvsingh@ucla.edu>
 * Copyright (C) 2012 by Tamer Crosby <tcrosby@ucla.edu>
 *
 * PIL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _SpatialParticle_h
#define _SpatialParticle_h

#include <boost/shared_ptr.hpp>
//#include <boost/utility.hpp> // defines shared_pointer<T>

#include <Eigen/Core>

#include <model/Utilities/StaticID.h>
#include <model/Utilities/CRTP.h>
#include <pil/SpatialCells/SpatialCellObserver.h>


namespace pil {

    template <typename Derived, int _dim>
    struct SpatialParticle :
    /* inheritance */ public model::StaticID<Derived>, // automatically generate unique ID(s) for SpatialParticle(s)
    /* inheritance */ public model::CRTP<Derived> // functions to cast this to derived type
    {
        
        enum{dim=_dim};        
        typedef Eigen::Matrix<double,dim,1> PositionType; // a vector of dim doubles
        
    protected:
        PositionType _p; // the spatial position of the particle
        
        
    public:
        
        typedef typename SpatialCellObserver<Derived,_dim>::SharedPtrType SharedPtrType;
        const SharedPtrType pCell;
        
        int rID;
        
        /*****************************************/
        SpatialParticle(const PositionType& posIN) :
        /* init list */ _p(posIN),
        /* init list */ pCell(SpatialCellObserver<Derived,_dim>::getCellByPosition(posIN)),
        /* init list */ rID(0)
        {/*! Constructor with input position
          * 1- initialize input position
          * 2- initialize shared pointer to SpatialCell (this either creates a new SpatialCell or points to an existing one)
          */
            pCell->addParticle(this->p_derived());
        }
        
        /* Destructor ********************************************/
		~SpatialParticle(){
			pCell->removeParticle(this->p_derived());
		}
        
        /*****************************************/
        const PositionType& P() const
        { /*! The position vector of this particle.
           */
            return _p;
        }
        
    };
    
}


#endif
