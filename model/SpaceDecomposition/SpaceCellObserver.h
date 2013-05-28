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
	template<typename SpaceCellType,short unsigned int dim>
	struct SpaceCellObserver{
		
		typedef Eigen::Matrix<double,dim,1> VectorDimD;
		typedef Eigen::Matrix<int,dim,1>    VectorDimI;
		typedef std::map<VectorDimI,SpaceCellType* const,CompareVectorsByComponent<int,dim> >  CellMapType;
		typedef std::shared_ptr<SpaceCellType> SharedPtrType;
        typedef std::pair<bool,SpaceCellType* const> isCellType;


        static double cellSize;
        
		/* begin() ***************************************************/
		static typename CellMapType::const_iterator begin()
        {
			return cellMap.begin();
		}
		
		/* end() *****************************************************/		
		static typename CellMapType::const_iterator end()
        {
			return cellMap.end();
		}
        
        static size_t size()
        {/*! @param[out] The number of observed SpaceCellType cells.
          */
            return cellMap.size();
        }

        /* getCellIDByPosition ************************************************/
		static VectorDimI getCellIDByPosition(const VectorDimD& P)
        {
			return floorEigen<dim>(P/cellSize);
		}
        
		/* getCellByID ********************************************************/
		static SharedPtrType getCellByID(const VectorDimI& cellID)
        {
			typename CellMapType::const_iterator iter(cellMap.find(cellID));
			return (iter!=cellMap.end())? (*(iter->second->particleContainer.begin()))->pCell : SharedPtrType(new SpaceCellType(cellID));
		}
		
		/* getCellByPosition **************************************************/		
		static SharedPtrType getCellByPosition(const VectorDimD& P)
        {
//			return getCellByID(floorEigen<dim>(P/cellSize));
			return getCellByID(getCellIDByPosition(P));
		}
        
        /* isCell *************************************************************/
        static isCellType isCell(const VectorDimI& cellID)
        {
            typename CellMapType::const_iterator iter(cellMap.find(cellID));
			return (iter!=cellMap.end())? std::make_pair(true,iter->second) : std::make_pair(false,(SpaceCellType*) NULL);
        }
        
        /*****************************************/
        template <class T>
        friend T& operator<< (T& os, const SpaceCellObserver<SpaceCellType,dim>& sCO)
        {/*! Operator << use ParticleType specific << operator
          */
            for (typename CellMapType::const_iterator cIter=sCO.begin();cIter!=sCO.end();++cIter)
            {
                os << (*cIter->second) << std::endl;
            }
            return os;
        }
		
		
	protected:
		static  CellMapType cellMap;		
	};

    /////////////////////////////
	// Declare static data member
	template <typename SpaceCellType,short unsigned int dim>
	double SpaceCellObserver<SpaceCellType,dim>::cellSize=1.0;

	template <typename SpaceCellType,short unsigned int dim>
	std::map<Eigen::Matrix<int,dim,1>,SpaceCellType* const,CompareVectorsByComponent<int,dim> > SpaceCellObserver<SpaceCellType,dim>::cellMap;
	
}	// close namespace
#endif

