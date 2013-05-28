/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _model_SpatialCellParticle_h_
#define _model_SpatialCellParticle_h_

#include <math.h>
#include <boost/utility.hpp>
#include <Eigen/Dense>
#include <model/SpaceDecomposition/SpatialCellObserver.h>
#include <model/SpaceDecomposition/SpatialCell.h>
#include <model/Utilities/TypeTraits.h>
#include <model/Utilities/CRTP.h>
#include <model/Utilities/StaticID.h>


namespace model {
	

	
	/********************************************************************************************/
	/********************************************************************************************/
	template<typename Derived, short unsigned int _dim>
	struct SpatialCellParticle : boost::noncopyable,
	/*                      */ private SpatialCellObserver<typename TypeTraits<Derived>::CellType,_dim>,
	/*                      */ public  CRTP<Derived>,
    /*                      */ public  StaticID<Derived>
    {

        enum{dim=_dim};
		//typedef SpatialCell<Derived,dim,cellSize> SpatialCellType;
        typedef typename TypeTraits<Derived>::CellType SpatialCellType;
		typedef typename SpatialCellType::ParticleContainerType ParticleContainerType;
		typedef SpatialCellObserver<SpatialCellType,dim> SpatialCellObserverType;
		typedef typename SpatialCellObserverType::CellMapType  CellMapType;
        typedef typename SpatialCellObserverType::VectorDimD  PositionType;
		typedef typename SpatialCellObserverType::VectorDimD  VectorDimD;
//		typedef typename SpatialCellObserverType::VectorDimI  VectorDimI;
		typedef typename SpatialCellObserverType::SharedPtrType  SharedPtrType;	
		
				
	public:
		
		//! The cell ID
//		const VectorDimI cellID;
		
		//! The pointer to the Cell
		const SharedPtrType pCell;

#ifdef _MODEL_MPI_
        int rID;
#endif

		
		/* Constructor **********************/
		SpatialCellParticle(const VectorDimD& P) :
//        /* init list */ cellID(floorEigen<dim>(P/cellSize)),
//        /* init list */ cellID(getCellIDByPosition(P)),
//		/* init list */ pCell(this->getCellByID(cellID))
		/* init list */ pCell(this->getCellByPosition(P))
        {
			pCell->addParticle(this->p_derived());
		}
		
		/* Destructor ********************************************/
		~SpatialCellParticle(){			
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

        /* neighborCellsBegin ***************************************/
        typename CellMapType::const_iterator neighborCellsBegin() const {
            return pCell->neighborCellsBegin();
        }
        
        /* neighborCellsEnd ***************************************/
        typename CellMapType::const_iterator neighborCellsEnd() const {
            return pCell->neighborCellsEnd();
        }
        
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
