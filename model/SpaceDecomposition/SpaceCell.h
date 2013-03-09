/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SPACECELL_H_
#define model_SPACECELL_H_

#include <assert.h>
#include <math.h>
#include <set>
#include <map>
//#include <boost/utility.hpp>
//#include <boost/shared_ptr.hpp>
//#include <memory> // std::shared_ptr (c++11)
#include <Eigen/Dense>
#include <model/Math/CompileTimeMath/Pow.h>
#include <model/SpaceDecomposition/SpaceCellObserver.h>
#include <model/SpaceDecomposition/NeighborShift.h>
#include <model/Utilities/CompareVectorsByComponent.h>


namespace model {
	
	/********************************************************************************************/
	/********************************************************************************************/
	/*! \brief A dim-dimensional cell occupying the spatial region cellID<= x/cellSize < (cellID+1).
	 *  SpaceCell is aware off all ParticleType objects present inside it.
	 */
	template<typename Derived, short unsigned int dim, double & cellSize>
	struct SpaceCell : boost::noncopyable,
	/*              */ private SpaceCellObserver<Derived,dim,cellSize>,
    /*              */ public  CRTP<Derived>{
        
        
 
		
        typedef typename TypeTraits<Derived>::ParticleType ParticleType;
        typedef Derived SpaceCellType;
		typedef SpaceCellObserver<SpaceCellType,dim,cellSize> SpaceCellObserverType;
		typedef typename SpaceCellObserverType::CellMapType  CellMapType;
		typedef typename SpaceCellObserverType::VectorDimD  VectorDimD;
		typedef typename SpaceCellObserverType::VectorDimI  VectorDimI;
		typedef std::set<const ParticleType*> ParticleContainerType; // PTR COMPARE IS NOT NECESSARY
        
        enum{neighborLayer=1}; // = (1+2*1)^dim  cells = 27  cells in 3d
        enum{    nearLayer=2}; // = (1+2*3)^dim  cells = 343 cells in 3d
        
        typedef NeighborShift<dim,1> NeighborShiftType;
        typedef NeighborShift<dim,1>     NearShiftType;

        
        CellMapType  neighborCells;
        CellMapType      nearCells;
        CellMapType       farCells;

        
        /* isNearCell *******************************************/
        bool isNearCell(const VectorDimI& otherCellID) const {
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
        bool isNeighborCell(const VectorDimI& otherCellID) const {
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
        
        //static int nearestNeighborOrder;
		
		//! The container of pointers to particles in this cell
		ParticleContainerType particleContainer;
				
        //! The ID of this cell, defining the dim-dimensional spatial region cellID<= x/cellSize < (cellID+1).
		const VectorDimI cellID;
        
        
        const VectorDimD center;
        
		//! The cellID(s) of the neighboring SpaceCell(s) (in column)
		const Eigen::Matrix<int,dim, NeighborShiftType::Nneighbors> neighborCellIDs;
		const Eigen::Matrix<int,dim,     NearShiftType::Nneighbors>     nearCellIDs;

        
		/* Constructor *******************************************/
		SpaceCell(const VectorDimI& cellID_in) :
        /* init list */ cellID(cellID_in),
        /* init list */ center((cellID.template cast<double>().array()+0.5).matrix()*cellSize),
        /* init list */ neighborCellIDs(NeighborShiftType::neighborIDs(cellID)),
        /* init list */     nearCellIDs(    NearShiftType::neighborIDs(cellID)){
            
			//! 1- Adds this to static SpaceCellObserver::cellMap
			assert(this->cellMap.insert(std::make_pair(cellID,this->p_derived())).second && "CANNOT INSERT SPACE CELL IN STATIC cellMap.");
            //! Populate nearCells and farCells
            for (typename CellMapType::const_iterator cellIter=this->begin();cellIter!=this->end();++cellIter){
                if (isNearCell(cellIter->second->cellID)){
                    if (isNeighborCell(cellIter->second->cellID)){
                        assert(                  neighborCells.insert(std::make_pair(cellIter->first,cellIter->second)).second && "CANNOT INSERT CELL IN NEIGHBORCELLS");
                        if(cellID!=cellIter->second->cellID){
                            assert(cellIter->second->neighborCells.insert(std::make_pair(cellID,this->p_derived())).second && "CANNOT INSERT THIS IN NEIGHBORCELLS");
                        }
                    }
                    else{
                        assert(                  nearCells.insert(std::make_pair(cellIter->first,cellIter->second)).second && "CANNOT INSERT CELL IN NEARCELLS");
                        assert(cellIter->second->nearCells.insert(std::make_pair(cellID,this->p_derived())).second && "CANNOT INSERT THIS IN NEARCELLS");
                    }
                }
                else{
                    assert(                  farCells.insert(std::make_pair(cellIter->first, cellIter->second)).second && "CANNOT INSERT CELL IN FARCELLS");
                    assert(cellIter->second->farCells.insert(std::make_pair(cellID,this->p_derived())).second && "CANNOT INSERT THIS IN FARCELLS");
                }
            }
            
		}
		
		/* Destructor *******************************************/
		~SpaceCell(){
			//! Removes this from static SpaceCellObserver::cellMap
			this->cellMap.erase(cellID);
			assert(particleContainer.empty() && "DESTROYING NON-EMPTY SPACE CELL.");
		}
		
		/* addParticle *******************************************/
		void addParticle(const ParticleType* const pP){
			//! 1- Adds pP to the particleContainer
			assert(particleContainer.insert(pP).second && "CANNOT INSERT PARTICLE IN SPACECELL");
		}
		
		/* removeParticle ****************************************/
		void removeParticle(const ParticleType* const pP){
			//! 1- Removes pP to the particleContainer
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
        
        
	};
    
    
    // Declare Static data
    //    template<typename ParticleType, short unsigned int dim, double & cellSize>
    //    int SpaceCell<ParticleType,dim,cellSize>::nearestNeighborOrder=1;
	
	////////////////////////////////////////////////////////////////////////////////
}	// close namespace model
#endif

