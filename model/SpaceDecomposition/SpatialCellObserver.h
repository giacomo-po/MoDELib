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

#ifndef model_SPACECELLOBSERVER_H_
#define model_SPACECELLOBSERVER_H_

#include <map>
#include <memory> // std::shared_ptr (c++11)
#include <model/Math/RoundEigen.h>
#include <model/SpaceDecomposition/SpatialCell.h>
#include <model/SpaceDecomposition/CellShift.h>

namespace model {
	
	/**************************************************************************/
	/**************************************************************************/
	template<typename ParticleType,short unsigned int dim>
	struct SpatialCellObserver
    {
        
        typedef Eigen::Matrix<double,dim,1> VectorDimD;
        typedef Eigen::Matrix<int,dim,1> CellIdType;
        typedef SpatialCell<ParticleType,dim> SpatialCellType;
		typedef std::map<CellIdType, SpatialCellType* const,CompareVectorsByComponent<int,dim> >  CellMapType;
		typedef std::shared_ptr<SpatialCellType> SharedPtrType;
        typedef std::pair<bool,SpatialCellType* const> isCellType;
        typedef CellShift<dim,1>    CellShiftType;
        
        
        
        
		/**********************************************************************/
		static typename CellMapType::const_iterator cellBegin()
        {/*! \returns A const_iterator to the first SpatialCellType cell.
          */
			return cellMap.begin();
		}
		
		/**********************************************************************/
		static typename CellMapType::const_iterator cellEnd()
        {/*! \returns A const_iterator to the past-the-last SpatialCellType cell.
          */
			return cellMap.end();
		}
        
        /**********************************************************************/
        static size_t totalCells()
        {/*! \returns The number of observed SpatialCellType cells.
          */
            return cellMap.size();
        }
        
		/**********************************************************************/
		static CellIdType getCellIDByPosition(const VectorDimD& P)
        {/*! \returns The CellIdType ID of the cell that contains P. The ID
          *  satisfies cellID <= P/cellSize < (cellID+1).
          */
            return floorEigen<dim>(P/_cellSize+VectorDimD::Constant(0.5));
		}
        
		/**********************************************************************/
		static SharedPtrType getCellByID(const CellIdType& cellID)
        {/*! @param[in] cellID The ID of the cell.
          *  \returns If a cell with ID=cellID exists, a shared-pointer to that
          *  cell is returned. Otherwise, a shared-pointer to a new cell with
          *  ID=cellID is returned.
          */
			typename CellMapType::const_iterator iter(cellMap.find(cellID));
			return (iter!=cellMap.end())? (*(iter->second->particleContainer.begin()))->pCell : SharedPtrType(new SpatialCellType(cellID));
		}
		
		/* getCellByPosition **************************************************/
		static SharedPtrType getCellByPosition(const VectorDimD& P)
        {/*! @param[in] P The position vector.
          *  \returns If a cell satisfying cellID <= P/cellSize < (cellID+1)
          *   exists, a shared-pointer to that cell is returned. Otherwise, a
          *   shared-pointer to a new cell satisfying cellID <= P/cellSize < (cellID+1)
          *   is returned.
          */
			return getCellByID(getCellIDByPosition(P));
		}
        
        /* isCell *************************************************************/
        static isCellType isCell(const CellIdType& cellID)
        {
            typename CellMapType::const_iterator iter(cellMap.find(cellID));
			return (iter!=cellMap.end())? std::make_pair(true,iter->second) : std::make_pair(false,( SpatialCellType*) NULL);
        }
        
        /**********************************************************************/
        static const CellMapType& cells()
        {
            return cellMap;
        }
        
        /**********************************************************************/
        static CellMapType neighborCells(const VectorDimD& P)
        {/*!@param[in] P the position vector
          *\returns a map of the SpatialCell(s) neighboring P.
          */
            const CellIdType cellID(getCellIDByPosition(P));
            const Eigen::Matrix<int,dim, CellShiftType::Nneighbors> neighborCellIDs(CellShiftType::neighborIDs(cellID));
            
            CellMapType temp;
            
            for (unsigned int c=0;c<CellShiftType::Nneighbors;++c)
            {
                isCellType isC(isCell(neighborCellIDs.col(c)));
                if (isC.first)
                {
                    const bool success=temp.emplace(isC.second->cellID,isC.second).second;
                    assert(success && "CANNOT INSERT CELL IN NEIGHBORCELLS");
                }
            }
            return temp;
        }
        
        /**********************************************************************/
        static CellMapType farCells(const VectorDimD& P)
        {/*!@param[in] P the position vector
          *\returns a map of the SpatialCell(s) not neighboring P.
          */
            const CellIdType cellID(getCellIDByPosition(P));
            const Eigen::Matrix<int,dim, CellShiftType::Nneighbors> neighborCellIDs(CellShiftType::neighborIDs(cellID));
            
            CellMapType temp(cellMap);
            
            for (unsigned int c=0;c<CellShiftType::Nneighbors;++c)
            {
                temp.erase(neighborCellIDs.col(c));
            }
            
            return temp;
        }
        
        /**********************************************************************/
        static void setCellSize(const double& newSize)
        {
            assert(newSize>0.0 && "cellSize MUST BE A POSITIVE double");
            assert(cellMap.empty() && "YOU ARE TRYING TO CHANGE cellSize WHILE SOME CELLS EXIST.");
            _cellSize=newSize;
        }
        
        /**********************************************************************/
        static double& cellSize()
        {
            return _cellSize;
        }
        
        /**********************************************************************/
        template <class T>
        friend T& operator<< (T& os, const SpatialCellObserver<SpatialCellType,dim>& sCO)
        {/*! Operator << uses ParticleType-specific operator <<
          */
            for (typename CellMapType::const_iterator cIter=sCO.begin();cIter!=sCO.end();++cIter)
            {
                os << (*cIter->second) << std::endl;
            }
            return os;
        }
        
        
    private:
        
        friend class SpatialCell<ParticleType,dim>;
        
        //! The container of SpatialCell pointers
        static  CellMapType cellMap;
        
        //! The size of a SpatialCell
        static double _cellSize;
        
    };
    
    // Static data member
    template <typename ParticleType,short unsigned int dim>
    double SpatialCellObserver<ParticleType,dim>::_cellSize=1.0;
    
    template <typename ParticleType,short unsigned int dim>
    std::map<Eigen::Matrix<int,dim,1>, SpatialCell<ParticleType,dim>* const,CompareVectorsByComponent<int,dim> > SpatialCellObserver<ParticleType,dim>::cellMap;
    
}
#endif
