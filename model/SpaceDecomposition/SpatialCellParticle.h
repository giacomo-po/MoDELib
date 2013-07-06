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
//#include <model/SpaceDecomposition/SpatialCell.h>
//#include <model/Utilities/TypeTraits.h>
#include <model/Utilities/CRTP.h>
#include <model/Utilities/StaticID.h>


namespace model {
	

	
	/********************************************************************************************/
	/********************************************************************************************/
	template<typename Derived, short unsigned int _dim>
	struct SpatialCellParticle :
    /*                      */ boost::noncopyable,
	/*                      */ public  CRTP<Derived>,
    /*                      */ public  StaticID<Derived>//,
//    /*                      */ private SpatialCellObserver<Derived,_dim>
    {

        enum{dim=_dim};
		//typedef SpatialCell<Derived,dim,cellSize> SpatialCellType;
//        typedef typename TypeTraits<Derived>::CellType SpatialCellType;
//        typedef SpatialCell<Derived,_dim> SpatialCellType;
//		typedef typename SpatialCellType::ParticleContainerType ParticleContainerType;
		typedef SpatialCellObserver<Derived,dim> SpatialCellObserverType;
//		typedef SpatialCellObserverType::ParticleContainerType ParticleContainerType;
		typedef typename SpatialCellObserverType::CellMapType  CellMapType;
        typedef typename SpatialCellObserverType::VectorDimD  PositionType;
		typedef typename SpatialCellObserverType::VectorDimD  VectorDimD;
//		typedef typename SpatialCellObserverType::VectorDimI  VectorDimI;
		typedef typename SpatialCellObserverType::SharedPtrType  SharedPtrType;	

        typedef  SpatialCell<Derived,_dim> SpatialCellType;

        typedef typename SpatialCellType::ParticleContainerType ParticleContainerType;
		
				
	public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        
        //! The position vector
        const VectorDimD P;
        
		//! The pointer to the Cell
		const SharedPtrType pCell;

        /**********************************************************************/
		SpatialCellParticle(const VectorDimD& Pin) :
		/* init list */ P(Pin),
		/* init list */ pCell(SpatialCellObserverType::getCellByPosition(P))
        {/*!\param[in] P the position of this SpatialCellParticle 
          */
			pCell->addParticle(this->p_derived());
		}
		
        
        /**********************************************************************/
        ~SpatialCellParticle()
        {/*! Destructor removes this from SpatialCell
          */
			pCell->removeParticle(this->p_derived());
		}
        
        /**********************************************************************/
        void moveTo(const VectorDimD& Pin)
        {/*!
          */
            assert(0 && "NEED TO TEST");
//            const SharedPtrType temp(SpatialCellObserverType::getCellByPosition(Pin));
//            if (pCell!=temp)
//            {
//                pCell->removeParticle(this->p_derived());
//                pCell=temp;
//                pCell->addParticle(this->p_derived());
//            }
		}

		
        /* neighborCellsBegin ***************************************/
        typename CellMapType::const_iterator neighborCellsBegin() const
        {/*!\returns a const iterator to the first neighbor SpatialCell
          */
            return pCell->neighborCellsBegin();
        }
        
        /* neighborCellsEnd ***************************************/
        typename CellMapType::const_iterator neighborCellsEnd() const
        {/*!\returns a const iterator to the past-the-end neighbor SpatialCell
          */
            return pCell->neighborCellsEnd();
        }
        
//        /* nearCellsBegin ***************************************/
//        typename CellMapType::const_iterator nearCellsBegin() const
//        {
//            return pCell->nearCellsBegin();
//        }
//        
//        /* nearCellsEnd ***************************************/
//        typename CellMapType::const_iterator nearCellsEnd() const
//        {
//            return pCell->nearCellsEnd();
//        }
//        
//        /* farCellsBegin ***************************************/
//        typename CellMapType::const_iterator farCellsBegin() const
//        {
//            return pCell->farCellsBegin();
//        }
//        
//        /* farCellsEnd ***************************************/
//        typename CellMapType::const_iterator farCellsEnd() const
//        {
//            return pCell->farCellsEnd();
//        }
        
            
        
        
	};
		
	/********************************************************************************************/
	/********************************************************************************************/
}	// close namespace
#endif


//		/* neighborBegin() ******************************************/
//		typename ParticleContainerType::const_iterator neighborBegin() const {
//			return pCell->neighborParticleContainer.begin();
//		}
//
//		/* neighborEnd() ******************************************/
//		typename ParticleContainerType::const_iterator neighborEnd() const {
//			return pCell->neighborParticleContainer.end();
//		}
