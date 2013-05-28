/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SPACECELLOBSERVER_H_
#define model_SPACECELLOBSERVER_H_

#include <map>
//#include <boost/shared_ptr.hpp>
#include <memory> // std::shared_ptr (c++11)

#include <Eigen/Dense>
#include <model/Utilities/CompareVectorsByComponent.h>

namespace model {
	
//	template <short unsigned int dim>
//	Eigen::Matrix<int,dim,1> floorEigen(const Eigen::Matrix<double,dim,1>& P) {
//		return (Eigen::Matrix<int,dim,1>()<< (int)std::floor(P(0)), floorEigen<dim-1>(P.template segment<dim-1>(1))).finished();
//	}
//	
//	template <>
//	Eigen::Matrix<int,1,1> floorEigen<1>(const Eigen::Matrix<double,1,1>& P) {
//		return (Eigen::Matrix<int,1,1>()<< (int)std::floor(P(0)) ).finished();
//	}
	
	/********************************************************************************************/
	/********************************************************************************************/	
	template<typename SpaceCellType,short unsigned int dim,double & cellSize>
	struct SpaceCellObserver{
		
		typedef Eigen::Matrix<double,dim,1> VectorDimD;
		typedef Eigen::Matrix<int,dim,1>    VectorDimI;
		typedef std::map<VectorDimI,SpaceCellType* const,CompareVectorsByComponent<int,dim> >  CellMapType;
		typedef std::shared_ptr<SpaceCellType> SharedPtrType;


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
			return (iter!=cellMap.end())? (*(iter->second->particleContainer.begin()))->pCell : SharedPtrType(new SpaceCellType(cellID));
		}
		
		/* getCellByPosition **********************************************/		
		static SharedPtrType getCellByPosition(const VectorDimD& P)  {
			return getCellByID(floorEigen<dim>(P/cellSize));
		}
		
		
	protected:
		static  CellMapType cellMap;		
	};
	
	/////////////////////////////
	// Declare static data member
	template <typename SpaceCellType,short unsigned int dim, double & cellSize>
	std::map<Eigen::Matrix<int,dim,1>,SpaceCellType* const,CompareVectorsByComponent<int,dim> > SpaceCellObserver<SpaceCellType,dim,cellSize>::cellMap;
	
}	// close namespace
#endif

