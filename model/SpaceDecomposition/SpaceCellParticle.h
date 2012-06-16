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
//#include <boost/shared_ptr.hpp>
#include <Eigen/Dense>
#include <model/SpaceDecomposition/SpaceCellObserver.h>
#include <model/SpaceDecomposition/SpaceCell.h>

namespace model {
	

	
	/********************************************************************************************/
	/********************************************************************************************/
	template<typename Derived, short unsigned int dim, double & cellSize>
	struct SpaceCellParticle : boost::noncopyable,
	/*                     */ private SpaceCellObserver<SpaceCell<Derived,dim,cellSize>,dim,cellSize>,
	/*                     */ public  CRTP<Derived>{ 

		typedef SpaceCell<Derived,dim,cellSize> SpaceCellType;
		typedef typename SpaceCellType::ParticleContainerType ParticleContainerType;

		typedef SpaceCellObserver<SpaceCellType,dim,cellSize> SpaceCellObserverType;
		typedef typename SpaceCellObserverType::CellMapType  CellMapType;	
		typedef typename SpaceCellObserverType::VectorDimD  VectorDimD;	
		typedef typename SpaceCellObserverType::VectorDimI  VectorDimI;	
		typedef typename SpaceCellObserverType::SharedPtrType  SharedPtrType;	
		
//		typedef boost::shared_ptr<SpaceCellType> SharedPtrType;

		
//		/* find_cell() **********************************************/ 
//		// THIS SHOULD BE A MEMBER FUNCTION OF SpaceCellObserver
//		SharedPtrType find_cell() const {
//			typename CellMapType::const_iterator iter(this->cellMap.find(cellID));
//			return (iter!=this->cellMap.end())? (*(iter->second->particleContainer.begin()))->pCell : SharedPtrType(new SpaceCellType(cellID));
//		}
		
		
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
		
		/* neighborBegin() ******************************************/
		typename ParticleContainerType::const_iterator neighborBegin() const {
			return pCell->neighborParticleContainer.begin();
		}

		/* neighborEnd() ******************************************/
		typename ParticleContainerType::const_iterator neighborEnd() const {
			return pCell->neighborParticleContainer.end();
		}
		
	};
		
	/********************************************************************************************/
	/********************************************************************************************/
}	// close namespace
#endif
