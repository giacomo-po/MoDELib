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

#ifndef _MODEL_SpatialCell_H_
#define _MODEL_SpatialCell_H_

#include <assert.h>
#include <math.h>
#include <set>
#include <map>
#include <boost/utility.hpp>
#include <boost/shared_ptr.hpp>
#include <Eigen/Dense>
#include <model/SpaceDecomposition/CellShift.h>


namespace model {
	
    // class pre-declaration
	template<typename ParticleType, int _dim>
	struct SpatialCellObserver;
    
	/********************************************************************************************/
	/********************************************************************************************/
	/*! \brief A dim-dimensional cell occupying the spatial region cellID<= x/cellSize < (cellID+1).
	 *  SpatialCell is aware off all ParticleType objects present inside it.
	 */
	template<typename ParticleType, int _dim>
	struct SpatialCell : boost::noncopyable,
	/*              */ private SpatialCellObserver<ParticleType,_dim>{
//    /*              */ public  CRTP<Derived>{
        
        
        enum{dim=_dim};
		
        //typedef typename TypeTraits<Derived>::ParticleType ParticleType;
        //typedef Derived SpatialCellType;
		typedef SpatialCellObserver<ParticleType,_dim> SpatialCellObserverType;
		typedef typename SpatialCellObserverType::CellMapType  CellMapType;
		typedef typename SpatialCellObserverType::VectorDimD  VectorDimD;
		typedef typename SpatialCellObserverType::VectorDimI  VectorDimI;
		typedef std::set<ParticleType*> ParticleContainerType; // PTR COMPARE IS NOT NECESSARY
        
        enum{neighborLayer=1}; // = (1+2*1)^dim  cells = 27  cells in 3d
        enum{    nearLayer=3}; // = (1+2*3)^dim  cells = 343 cells in 3d
        
        typedef CellShift<dim,1> CellShiftType;
        typedef CellShift<dim,1>     NearShiftType;

        
        CellMapType  neighborCells;
//        CellMapType      nearCells;
//        CellMapType       farCells;

        int assignedRank;
        
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
		
        //, double & cellSize
        
		//! The container of pointers to particles in this cell
		ParticleContainerType particleContainer;
				
        //! The ID of this cell, defining the dim-dimensional spatial region cellID<= x/cellSize < (cellID+1).
		const VectorDimI cellID;
        const double& cellSize;
        
        const VectorDimD center;
        
		//! The cellID(s) of the neighboring SpatialCell(s) (in column)
		const Eigen::Matrix<int,dim, CellShiftType::Nneighbors> neighborCellIDs;
		const Eigen::Matrix<int,dim,     NearShiftType::Nneighbors>     nearCellIDs;

        
		/* Constructor *******************************************/
		SpatialCell(const VectorDimI& cellID_in) :
        /* init list */ cellID(cellID_in),
        /* init list */ cellSize(SpatialCellObserver<ParticleType,dim>::cellSize),
        /* init list */ center((cellID.template cast<double>().array()+0.5).matrix()*cellSize),
        /* init list */ neighborCellIDs(CellShiftType::neighborIDs(cellID)),
        /* init list */     nearCellIDs(    NearShiftType::neighborIDs(cellID)),
        /* init list */ assignedRank(0){
                        
			//! 1- Adds this to static SpatialCellObserver::cellMap
			assert(this->cellMap.insert(std::make_pair(cellID,this)).second && "CANNOT INSERT SPACE CELL IN STATIC cellMap.");
            //! Populate nearCells and farCells
            
            
            for (int c=0;c<neighborCellIDs.cols();++c)
            {
                typename SpatialCellObserverType::isCellType isC(SpatialCellObserverType::isCell(neighborCellIDs.col(c)));
                if (isC.first)
                {
                    assert(neighborCells.insert(std::make_pair(isC.second->cellID,isC.second)).second && "CANNOT INSERT CELL IN NEIGHBORCELLS");
                    if(cellID!=isC.second->cellID){
                        assert(isC.second->neighborCells.insert(std::make_pair(cellID,this)).second && "CANNOT INSERT THIS IN NEIGHBORCELLS");
                    }
                }
            }
            
            
//            for (typename CellMapType::const_iterator cellIter=this->begin();cellIter!=this->end();++cellIter){
//              //  if (isNearCell(cellIter->second->cellID)){
//                    if (isNeighborCell(cellIter->second->cellID)){
//                        assert(                  neighborCells.insert(std::make_pair(cellIter->first,cellIter->second)).second && "CANNOT INSERT CELL IN NEIGHBORCELLS");
//                        if(cellID!=cellIter->second->cellID){
//                            assert(cellIter->second->neighborCells.insert(std::make_pair(cellID,this)).second && "CANNOT INSERT THIS IN NEIGHBORCELLS");
//                        }
//                    }
//               //     else{
//               //         assert(                  nearCells.insert(std::make_pair(cellIter->first,cellIter->second)).second && "CANNOT INSERT CELL IN NEARCELLS");
//               //         assert(cellIter->second->nearCells.insert(std::make_pair(cellID,this)).second && "CANNOT INSERT THIS IN NEARCELLS");
//               //     }
//               // }
//               // else{
//               //     assert(                  farCells.insert(std::make_pair(cellIter->first, cellIter->second)).second && "CANNOT INSERT CELL IN FARCELLS");
//               //     assert(cellIter->second->farCells.insert(std::make_pair(cellID,this)).second && "CANNOT INSERT THIS IN FARCELLS");
//               // }
//            }
            
		}
		
		/* Destructor *******************************************/
		~SpatialCell()
        {
			//! Removes this from static SpatialCellObserver::cellMap
			this->cellMap.erase(cellID);
			assert(particleContainer.empty() && "DESTROYING NON-EMPTY SPACE CELL.");
            
            //! Removes this from neighbor cells
            for (int c=0;c<neighborCellIDs.cols();++c)
            {
                typename SpatialCellObserverType::isCellType isC(SpatialCellObserverType::isCell(neighborCellIDs.col(c)));
                if (isC.first)
                {
//                    assert(                 neighborCells.insert(std::make_pair(isC.second->cellID,isC.second)).second && "CANNOT INSERT CELL IN NEIGHBORCELLS");
                    if(cellID!=isC.second->cellID){
                        assert(isC.second->neighborCells.erase(cellID)==1 && "CANNOT ERASE CURRENT CELL FROM NEIGHBOR CELL");
                    }
                }
            }
            
		}
		
		/* addParticle *******************************************/
		void addParticle(ParticleType* const pP)
        {
			//! 1- Adds pP to the particleContainer
			assert(particleContainer.insert(pP).second && "CANNOT INSERT PARTICLE IN SpatialCell");
		}
		
		/* removeParticle ****************************************/
		void removeParticle(ParticleType* const pP)
        {
			//! 1- Removes pP to the particleContainer
			assert(particleContainer.erase(pP)==1 && "CANNOT ERASE PARTICLE FROM particleContainer.");
		}

        /* neighborCellBegin ***************************************/
        typename CellMapType::const_iterator neighborCellsBegin() const
        {
            return neighborCells.begin();
        }
        
        /* neighborCellEnd ***************************************/
        typename CellMapType::const_iterator neighborCellsEnd() const
        {
            return neighborCells.end();
        }
        
//        /* nearCellBegin ***************************************/
//        typename CellMapType::const_iterator nearCellsBegin() const
//        {
//            return nearCells.begin();
//        }
//        
//        /* nearCellEnd ***************************************/
//        typename CellMapType::const_iterator nearCellsEnd() const
//        {
//            return nearCells.end();
//        }
//        
//        /* nearCellBegin ***************************************/
//        typename CellMapType::const_iterator farCellsBegin() const
//        {
//            return farCells.begin();
//        }
//        
//        /* nearCellEnd ***************************************/
//        typename CellMapType::const_iterator farCellsEnd() const
//        {
//            return farCells.end();
//        }
        
        /* particleBegin ***************************************/
        typename ParticleContainerType::const_iterator particleBegin() const
        {
            return particleContainer.begin();
        }
        
        /* particleEnd *****************************************/
        typename ParticleContainerType::const_iterator particleEnd() const
        {
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
        
//        /* nearSize ********************************************/
//        size_t nearSize() const {
//            size_t temp(0);
//            for (typename CellMapType::const_iterator cIter=nearCellsBegin(); cIter!=nearCellsEnd(); ++cIter)
//            {
//                temp+=cIter->second->size();
//            }
//            return temp;
//        }
//        
//        /* farSize ********************************************/
//        size_t farSize() const {
//            size_t temp(0);
//            for (typename CellMapType::const_iterator cIter=farCellsBegin(); cIter!=farCellsEnd(); ++cIter)
//            {
//                temp+=cIter->second->size();
//            }
//            return temp;
//        }
        
//        /*****************************************/
//        void setAssignedRank(const int& rank)
//        {
//            assignedRank=rank;
//        }
//        
//        /*****************************************/
//        void setAssignedRank(const int& rank)
//        {
//            assignedRank=rank;
//        }
        
        /*****************************************/
        template <class T>
        friend T& operator << (T& os, const SpatialCell<ParticleType,_dim>& sC)
        {/*! operator << is used to output SpatialCell info
          */
            os<<sC.cellID.transpose()<<" "<< sC.size()<<" "<<sC.neighborSize()<< " "<<sC.assignedRank;
            return os;
        }
        
        
	};
    
    
    // Declare Static data
    //    template<typename ParticleType, short unsigned int dim, double & cellSize>
    //    int SpatialCell<ParticleType,dim,cellSize>::nearestNeighborOrder=1;
	
	////////////////////////////////////////////////////////////////////////////////
}	// close namespace pil
#endif

