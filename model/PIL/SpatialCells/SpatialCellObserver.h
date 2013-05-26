/* This file is part of PIL, the Particle Interaction Library.
 *
 * Copyright (C) 2012 by Giacomo Po <gpo@ucla.edu>
 * Copyright (C) 2012 by Shao-Ching Huang <sch@ucla.edu>
 * Copyright (C) 2012 by Tajendra Singh <tvsingh@ucla.edu>
 * Copyright (C) 2012 by Tamer Crosby <tcrosby@ucla.edu>
 *
 * PIL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef pil_SpatialCellObserver_H_
#define pil_SpatialCellObserver_H_

#include <map>
#include <boost/shared_ptr.hpp>

#include <Eigen/Dense>
#include <model/Utilities/CompareVectorsByComponent.h>
#include <model/PIL/SpatialCells/SpatialCell.h>

namespace pil {
	
	template <short unsigned int dim>
	Eigen::Matrix<int,dim,1> floorEigen(const Eigen::Matrix<double,dim,1>& P) {
		return (Eigen::Matrix<int,dim,1>()<< (int)std::floor(P(0)), floorEigen<dim-1>(P.template segment<dim-1>(1))).finished();
	}
	
	template <>
	Eigen::Matrix<int,1,1> floorEigen<1>(const Eigen::Matrix<double,1,1>& P) {
		return (Eigen::Matrix<int,1,1>()<< (int)std::floor(P(0)) ).finished();
	}
	
	/********************************************************************************************/
	/********************************************************************************************/	
	template<typename ParticleType, int _dim>
	struct SpatialCellObserver{
		
        enum{dim=_dim}; // make dim available ouside class
		typedef Eigen::Matrix<double,dim,1> VectorDimD;
		typedef Eigen::Matrix<int,dim,1>    VectorDimI;
        //typedef SpatialCell<ParticleType> SpatialCellType;
        typedef Eigen::Matrix<int,_dim,1>    CellIdType;

		typedef std::map<Eigen::Matrix<int,_dim,1>,SpatialCell<ParticleType,_dim>* const,model::CompareVectorsByComponent<int,_dim> >  CellMapType;
		typedef boost::shared_ptr<SpatialCell<ParticleType,dim> > SharedPtrType;
        
        typedef std::pair<bool,SpatialCell<ParticleType,_dim>* const> isCellType;

//std::map<Eigen::Matrix<int,ParticleType::dim,1>,SpatialCell<ParticleType>* const,CompareVectorsByComponent<int,ParticleType::dim> >
//std::map<Eigen::Matrix<int,ParticleType::dim,1>,SpatialCell<ParticleType>* const,CompareVectorsByComponent<int,ParticleType::dim> >
        
        static double cellSize;
        
		/* begin() ***************************************************/
		static typename CellMapType::const_iterator begin() {
			return cellMap.begin();
		}
		
		/* end() *****************************************************/		
		static typename CellMapType::const_iterator end() {
			return cellMap.end();
		}
        
        static size_t size() {
            return cellMap.size();
        }
		
		/* getCellByID *****************************************************/
		static SharedPtrType getCellByID(const VectorDimI& cellID)  {
			typename CellMapType::const_iterator iter(cellMap.find(cellID));
			return (iter!=cellMap.end())? (*(iter->second->particleContainer.begin()))->pCell : SharedPtrType(new SpatialCell<ParticleType,dim>(cellID));
		}
		
		/* getCellByPosition **********************************************/		
		static SharedPtrType getCellByPosition(const VectorDimD& P)  {
			return getCellByID(floorEigen<dim>(P/cellSize));
		}
        
        /* getCellByID *****************************************************/
		static SpatialCell<ParticleType,dim>& getExistingCellByID(const VectorDimI& cellID)  {
			typename CellMapType::const_iterator iter(cellMap.find(cellID));
            assert(iter!=cellMap.end() && "REQUESTED CELL DOES NOT EXIST.");
            return *iter->second;
//			return (iter!=cellMap.end())? (*(iter->second->particleContainer.begin()))->pCell : SharedPtrType(new SpatialCell<ParticleType,dim>(cellID));
		}
        
//        /* getCellByID *****************************************************/
//		static const SpatialCell<ParticleType,dim>& getExistingCellByID(const VectorDimI& cellID)  {
//			typename CellMapType::const_iterator iter(cellMap.find(cellID));
//            assert(iter!=cellMap.end() && "REQUESTED CELL DOES NOT EXIST.");
//            return *iter->second;
//            //			return (iter!=cellMap.end())? (*(iter->second->particleContainer.begin()))->pCell : SharedPtrType(new SpatialCell<ParticleType,dim>(cellID));
//		}
        
        /* isCell ************************************************/
        static isCellType isCell(const VectorDimI& cellID)
        {
            typename CellMapType::const_iterator iter(cellMap.find(cellID));
			return (iter!=cellMap.end())? std::make_pair(true,iter->second) : std::make_pair(false,(SpatialCell<ParticleType,_dim>*) NULL);
        }
        
        /*****************************************/
        template <class T>
        friend T& operator<< (T& os, const SpatialCellObserver<ParticleType,dim>& sCO)
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
	template<typename ParticleType, int dim>
                std::map<Eigen::Matrix<int,dim,1>,SpatialCell<ParticleType,dim>* const,model::CompareVectorsByComponent<int,dim> > SpatialCellObserver<ParticleType,dim>::cellMap;

	template<typename ParticleType, int dim>
	double SpatialCellObserver<ParticleType,dim>::cellSize=1.0;

	
}	// close namespace
#endif

