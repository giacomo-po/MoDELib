/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SpatialCell_h_
#define model_SpatialCell_h_

#include <assert.h>
#include <math.h>
#include <set>
//#include <map>
//#include <boost/utility.hpp>
//#include <boost/shared_ptr.hpp>
//#include <memory> // std::shared_ptr (c++11)
#include <Eigen/Dense>
//#include <model/Math/CompileTimeMath/Pow.h>
#include <model/Utilities/CRTP.h>
#include <model/SpaceDecomposition/SpatialCellObserver.h>
//#include <model/SpaceDecomposition/NeighborShift.h>
#include <model/SpaceDecomposition/CellShift.h>
//#include <model/Utilities/CompareVectorsByComponent.h>


namespace model {
	
	/********************************************************************************************/
	/********************************************************************************************/
	/*! \brief A dim-dimensional cell occupying the spatial region cellID<= x/cellSize < (cellID+1).
	 *  SpatialCell is aware off all ParticleType objects present inside it.
	 */
	template<typename Derived, short unsigned int dim>
	struct SpatialCell : boost::noncopyable,
	/*              */ private SpatialCellObserver<Derived,dim>,
    /*              */ public  CRTP<Derived>{
        
        
 
		
        typedef typename TypeTraits<Derived>::ParticleType ParticleType;
        typedef Derived SpatialCellType;
		typedef SpatialCellObserver<SpatialCellType,dim> SpatialCellObserverType;
		typedef typename SpatialCellObserverType::CellMapType  CellMapType;
		typedef typename SpatialCellObserverType::VectorDimD  VectorDimD;
		typedef typename SpatialCellObserverType::CellIdType  CellIdType;
		typedef std::set<const ParticleType*> ParticleContainerType; // PTR COMPARE IS NOT NECESSARY
//		typedef std::set<ParticleType*> ParticleContainerType; // PTR COMPARE IS NOT NECESSARY

        
        enum{neighborLayer=1}; // = (1+2*1)^dim  cells = 27  cells in 3d
        enum{    nearLayer=2}; // = (1+2*3)^dim  cells = 343 cells in 3d
        
        typedef CellShift<dim,1> CellShiftType;
        typedef CellShift<dim,1>     NearShiftType;

        
        CellMapType  neighborCells;
        CellMapType      nearCells;
        CellMapType       farCells;
        
#ifdef _MODEL_MPI_
        int assignedRank;
#endif

        
        /* isNearCell *******************************************/
        bool isNearCell(const CellIdType& otherCellID) const {
            bool temp(false);
            for (int k=0;k<nearCellIDs.cols();++k){
                if(nearCellIDs.col(k)==otherCellID){ // cell is a neighbor
                    temp=true;
                    break;
                }
            }
            return temp;
        }
        
        /* isNeighborCell ***************************************/
        bool isNeighborCell(const CellIdType& otherCellID) const {
            bool temp(false);
            for (int k=0;k<neighborCellIDs.cols();++k){
                if(neighborCellIDs.col(k)==otherCellID){ // cell is a neighbor
                    temp=true;
                    break;
                }
            }
            return temp;
        }
        
        
	public:
        		
		//! The container of pointers to particles in this cell
		ParticleContainerType particleContainer; // TO DO: MAKE THIS A BASE CLASS
				
        //! The ID of this cell, defining the dim-dimensional spatial region cellID<= x/cellSize < (cellID+1).
		const CellIdType cellID;
        
        
        const VectorDimD center;
        
		//! The cellID(s) of the neighboring SpatialCell(s) (in column)
		const Eigen::Matrix<int,dim, CellShiftType::Nneighbors> neighborCellIDs;
		const Eigen::Matrix<int,dim,     NearShiftType::Nneighbors>     nearCellIDs;

        
		/* Constructor *******************************************/
		SpatialCell(const CellIdType& cellID_in) :
        /* init list */ cellID(cellID_in),
        /* init list */ center((cellID.template cast<double>().array()+0.5).matrix()*this->cellSize),
        /* init list */ neighborCellIDs(CellShiftType::neighborIDs(cellID)),
        /* init list */     nearCellIDs(    NearShiftType::neighborIDs(cellID))
        {/*! @param[in] cellID_in The ID of the current cell.
          * The constructor initializes cellID, center, neighborCellIDs, and
          * nearCellIDs. It then populates neighborCells, nearCells, and
          * farCells.
          */
            
			//! 1- Adds this to static SpatialCellObserver::cellMap
			assert(this->cellMap.insert(std::make_pair(cellID,this->p_derived())).second && "CANNOT INSERT Spatial CELL IN STATIC cellMap.");
            //! Populate nearCells and farCells

            for (int c=0;c<neighborCellIDs.cols();++c)
            {
                typename SpatialCellObserverType::isCellType isC(SpatialCellObserverType::isCell(neighborCellIDs.col(c)));
                if (isC.first)
                {
                    assert(neighborCells.insert(std::make_pair(isC.second->cellID,isC.second)).second && "CANNOT INSERT CELL IN NEIGHBORCELLS");
                    if(cellID!=isC.second->cellID){
                        assert(isC.second->neighborCells.insert(std::make_pair(cellID,this->p_derived())).second && "CANNOT INSERT THIS IN NEIGHBORCELLS");
                    }
                }
            }
            
//            for (typename CellMapType::const_iterator cellIter=this->cellBegin();cellIter!=this->cellEnd();++cellIter)
//            {
//                if (isNearCell(cellIter->second->cellID))
//                {
//                    if (isNeighborCell(cellIter->second->cellID))
//                    {
//                        assert(                  neighborCells.insert(std::make_pair(cellIter->first,cellIter->second)).second && "CANNOT INSERT CELL IN NEIGHBORCELLS");
//                        if(cellID!=cellIter->second->cellID)
//                        {
//                            assert(cellIter->second->neighborCells.insert(std::make_pair(cellID,this->p_derived())).second && "CANNOT INSERT THIS IN NEIGHBORCELLS");
//                        }
//                    }
//                    else
//                    {
//                        assert(                  nearCells.insert(std::make_pair(cellIter->first,cellIter->second)).second && "CANNOT INSERT CELL IN NEARCELLS");
//                        assert(cellIter->second->nearCells.insert(std::make_pair(cellID,this->p_derived())).second && "CANNOT INSERT THIS IN NEARCELLS");
//                    }
//                }
//                else
//                {
//                    assert(                  farCells.insert(std::make_pair(cellIter->first, cellIter->second)).second && "CANNOT INSERT CELL IN FARCELLS");
//                    assert(cellIter->second->farCells.insert(std::make_pair(cellID,this->p_derived())).second && "CANNOT INSERT THIS IN FARCELLS");
//                }
//            }
            
		}
		
		/* Destructor *******************************************/
		~SpatialCell()
        {/*! Removes this from static SpatialCellObserver::cellMap
          */
			this->cellMap.erase(cellID);
			assert(particleContainer.empty() && "DESTROYING NON-EMPTY Spatial CELL.");
		}
		
		/* addParticle *******************************************/
		void addParticle(const ParticleType* const pP)
        {/*! Adds pP to the particleContainer
          */
			assert(particleContainer.insert(pP).second && "CANNOT INSERT PARTICLE IN SpatialCELL");
		}
		
		/* removeParticle ****************************************/
		void removeParticle(const ParticleType* const pP)
        {/*! Removes pP from the particleContainer
          */
			assert(particleContainer.erase(pP)==1 && "CANNOT ERASE PARTICLE FROM particleContainer.");
		}

        /* neighborCellBegin ***************************************/
        typename CellMapType::const_iterator neighborCellsBegin() const {
            return neighborCells.begin();
        }
        
        /* neighborCellEnd ***************************************/
        typename CellMapType::const_iterator neighborCellsEnd() const {
            return neighborCells.end();
        }
        
        /* nearCellBegin ***************************************/
        typename CellMapType::const_iterator nearCellsBegin() const {
            return nearCells.begin();
        }
        
        /* nearCellEnd ***************************************/
        typename CellMapType::const_iterator nearCellsEnd() const {
            return nearCells.end();
        }
        
        /* nearCellBegin ***************************************/
        typename CellMapType::const_iterator farCellsBegin() const {
            return farCells.begin();
        }
        
        /* nearCellEnd ***************************************/
        typename CellMapType::const_iterator farCellsEnd() const {
            return farCells.end();
        }
        
        /* particleBegin ***************************************/
        typename ParticleContainerType::const_iterator particleBegin() const {
            return particleContainer.begin();
        }
        
        /* particleEnd *****************************************/
        typename ParticleContainerType::const_iterator particleEnd() const {
            return particleContainer.end();
        }
        
        /* size ************************************************/
        size_t size() const
        {
            return particleContainer.size();
        }
        
        /* neighborSize ****************************************/
        size_t neighborSize() const
        {
            size_t temp(0);
            for (typename CellMapType::const_iterator cIter=neighborCellsBegin(); cIter!=neighborCellsEnd(); ++cIter)
            {
                temp+=cIter->second->size();
            }
            return temp;
        }
        
        /* size ************************************************/
        size_t n2Weight() const
        {
            return size()*neighborSize();
        }
        
        
	};
    
	
	////////////////////////////////////////////////////////////////////////////////
}	// close namespace model
#endif

