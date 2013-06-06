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
#include <Eigen/Dense>
#include <model/Utilities/CompareVectorsByComponent.h>

namespace model {
	
	
	/********************************************************************************************/
	/********************************************************************************************/	
	template<typename SpatialCellType,short unsigned int dim>
	struct SpatialCellObserver
    {
		
		typedef Eigen::Matrix<double,dim,1> VectorDimD;
        typedef Eigen::Matrix<int,dim,1> CellIdType;
		typedef std::map<CellIdType, SpatialCellType* const,CompareVectorsByComponent<int,dim> >  CellMapType;
		typedef std::shared_ptr< SpatialCellType> SharedPtrType;
        typedef std::pair<bool, SpatialCellType* const> isCellType;


        static double cellSize;
        
		/* begin() ***************************************************/
		static typename CellMapType::const_iterator cellBegin()
        {/*! \returns A const_iterator to the first SpatialCellType cell. 
          */
			return cellMap.begin();
		}
		
		/* end() *****************************************************/		
		static typename CellMapType::const_iterator cellEnd()
        {/*! \returns A const_iterator to the past-the-last SpatialCellType cell.
          */
			return cellMap.end();
		}
        
        static size_t size()
        {/*! \returns The number of observed SpatialCellType cells.
          */
            return cellMap.size();
        }

        /* getCellIDByPosition ************************************************/
		static CellIdType getCellIDByPosition(const VectorDimD& P)
        {/*! \returns The CellIdType ID of the cell that contains P. The ID
          *  satisfies cellID <= P/cellSize < (cellID+1).
          */
			return floorEigen<dim>(P/cellSize);
		}
        
		/* getCellByID ********************************************************/
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
        
//        /*****************************************/
//        template <class T>
//        friend T& operator<< (T& os, const SpatialCellObserver<SpatialCellType,dim>& sCO)
//        {/*! Operator << use ParticleType specific << operator
//          */
//            for (typename CellMapType::const_iterator cIter=sCO.begin();cIter!=sCO.end();++cIter)
//            {
//                os << (*cIter->second) << std::endl;
//            }
//            return os;
//        }
		
		
	protected:
		static  CellMapType cellMap;		
	};

    /////////////////////////////
	// Declare static data member
	template <typename SpatialCellType,short unsigned int dim>
	double SpatialCellObserver< SpatialCellType,dim>::cellSize=1.0;

	template <typename SpatialCellType,short unsigned int dim>
	std::map<Eigen::Matrix<int,dim,1>, SpatialCellType* const,CompareVectorsByComponent<int,dim> > SpatialCellObserver< SpatialCellType,dim>::cellMap;
	
}	// close namespace
#endif

