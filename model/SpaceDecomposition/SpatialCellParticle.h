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
#include <Eigen/Dense>
#include <model/SpaceDecomposition/SpatialCellObserver.h>
#include <model/Utilities/CRTP.h>
#include <model/Utilities/StaticID.h>
#include <model/Utilities/NonCopyable.h>
#include <model/Math/CompileTimeMath/CTM.h>



namespace model
{
	/**************************************************************************/
	/**************************************************************************/
	template<typename Derived, short unsigned int _dim>
	struct SpatialCellParticle : public NonCopyable,
	/*                        */ public  CRTP<Derived>,
    /*                        */ public  StaticID<Derived>
    {

        enum{dim=_dim};
		typedef SpatialCellObserver<Derived,dim> SpatialCellObserverType;
		typedef typename SpatialCellObserverType::CellMapType  CellMapType;
        typedef typename SpatialCellObserverType::VectorDimD  PositionType;
		typedef typename SpatialCellObserverType::VectorDimD  VectorDimD;
		typedef typename SpatialCellObserverType::SharedPtrType  SharedPtrType;

        typedef  SpatialCell<Derived,_dim> SpatialCellType;

        typedef typename SpatialCellType::ParticleContainerType ParticleContainerType;
        typedef Eigen::Matrix<double,1,CTM::pow(2,dim)> VectorVerticesType;

				
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

        /* neighborCellsBegin *************************************************/
        typename CellMapType::const_iterator neighborCellsBegin() const
        {/*!\returns a const iterator to the first neighbor SpatialCell
          */
            return pCell->neighborCellsBegin();
        }
        
        /* neighborCellsEnd ***************************************************/
        typename CellMapType::const_iterator neighborCellsEnd() const
        {/*!\returns a const iterator to the past-the-end neighbor SpatialCell
          */
            return pCell->neighborCellsEnd();
        }
        
        /**********************************************************************/
        VectorVerticesType vertexWeigths() const
        {
            return pCell->vertexWeigths(P);
        }
        
	};
		
}	// close namespace
#endif




