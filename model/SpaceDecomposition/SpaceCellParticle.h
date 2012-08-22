/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SPACECELLPARTICLE_H_
#define model_SPACECELLPARTICLE_H_

#include <math.h>
#include <boost/utility.hpp>
#include <Eigen/Dense>
#include <model/SpaceDecomposition/SpaceCellObserver.h>
#include <model/SpaceDecomposition/SpaceCell.h>

namespace model {
	

	
	/********************************************************************************************/
	/********************************************************************************************/
	template<typename Derived, short unsigned int dim, double & cellSize>
	struct SpaceCellParticle : boost::noncopyable,
	/*                      */ private SpaceCellObserver<SpaceCell<Derived,dim,cellSize>,dim,cellSize>,
	/*                      */ public  CRTP<Derived>{ 

		typedef SpaceCell<Derived,dim,cellSize> SpaceCellType;
		typedef typename SpaceCellType::ParticleContainerType ParticleContainerType;
		typedef SpaceCellObserver<SpaceCellType,dim,cellSize> SpaceCellObserverType;
		typedef typename SpaceCellObserverType::CellMapType  CellMapType;	
		typedef typename SpaceCellObserverType::VectorDimD  VectorDimD;	
		typedef typename SpaceCellObserverType::VectorDimI  VectorDimI;	
		typedef typename SpaceCellObserverType::SharedPtrType  SharedPtrType;	
		
				
	public:
		
		//! The cell ID
		const VectorDimI cellID;
		
		//! The pointer to the Cell
		const SharedPtrType pCell;
		
		/* Constructor **********************/
		SpaceCellParticle(const VectorDimD& P) : cellID(floorEigen<dim>(P/cellSize)),
		/*                                    */ pCell(this->getCellByID(cellID)){
			pCell->addParticle(this->p_derived());
		}
		
		/* Destructor ********************************************/
		~SpaceCellParticle(){			
			pCell->removeParticle(this->p_derived());
		}
		
//		/* neighborBegin() ******************************************/
//		typename ParticleContainerType::const_iterator neighborBegin() const {
//			return pCell->neighborParticleContainer.begin();
//		}
//
//		/* neighborEnd() ******************************************/
//		typename ParticleContainerType::const_iterator neighborEnd() const {
//			return pCell->neighborParticleContainer.end();
//		}
		
        /* nearCellsBegin ***************************************/
        typename CellMapType::const_iterator nearCellsBegin() const {
            return pCell->nearCellsBegin();
        }
        
        /* nearCellsEnd ***************************************/
        typename CellMapType::const_iterator nearCellsEnd() const {
            return pCell->nearCellsEnd();
        }
        
        /* farCellsBegin ***************************************/
        typename CellMapType::const_iterator farCellsBegin() const {
            return pCell->farCellsBegin();
        }
        
        /* farCellsEnd ***************************************/
        typename CellMapType::const_iterator farCellsEnd() const {
            return pCell->farCellsEnd();
        }
        
	};
		
	/********************************************************************************************/
	/********************************************************************************************/
}	// close namespace
#endif
