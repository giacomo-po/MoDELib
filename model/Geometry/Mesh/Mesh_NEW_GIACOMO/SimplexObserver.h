// This file is part of mmdl, the C++ materials defect mechanics library.
//
// Copyright (C) 2011 by Giacomo Po <giacomopo@gmail.com>.
//
// mmdl is distributed without any warranty under the 
// GNU Lesser General Public License <http://www.gnu.org/licenses/>.


#ifndef mmdl_ADDRESSBOOK_H_
#define mmdl_ADDRESSBOOK_H_

#include <map>
#include <boost/shared_ptr.hpp>

#include <Eigen/Dense>
#include <mmdl/Utilities/CompareVectorsByComponent.h>

#include "mmdl/Geometry/Mesh/Mesh_NEW_GIACOMO/SimplexID.h"


namespace mmdl {
	
	
	template<short int ambientDim, short int sOrder>
	class Simplex;
	
	/********************************************************************************************/
	/********************************************************************************************/	
	template<short int ambientDim, short int sOrder>
	struct SimplexObserver{
		
//		typedef Eigen::Matrix<double,dim,1> VectorDimD;
//		typedef Eigen::Matrix<int,sOrder+1,1>    SimplexIdType;
		typedef std::map<SimplexID<sOrder>,Simplex<ambientDim,sOrder>* const,CompareVectorsByComponent<int,sOrder+1> >  simplexMapType; // This should be a set. 
		
		typedef boost::shared_ptr<Simplex<ambientDim,sOrder> > SharedPtrType;
		
		
		/* begin() ***************************************************/
		static typename simplexMapType::const_iterator begin() {
			return simplexMap.begin();
		}
		
		/* end() *****************************************************/		
		static typename simplexMapType::const_iterator end() {
			return simplexMap.end();
		}
		
//		/* getCellByID *****************************************************/		// THIS GOES IN MESH
//		static SharedPtrType getSimplex(const VectorDimI& simplexID)  {
//			typename simplexMapType::const_iterator iter(simplexMap.find(simplexID));
//			return (iter!=simplexMap.end())? (*(iter->second->particleContainer.begin()))->pCell : SharedPtrType(new SpaceCellType(cellID));
//		}
		
		/* getCellByPosition **********************************************/		
//		static SharedPtrType getCellByPosition(const VectorDimD& P)  {
//			//			typename simplexMapType::const_iterator iter(simplexMap.find(cellID));
//			return getCellByID(floorEigen<dim>(P));
//		}
		
		
	protected:
		static  simplexMapType simplexMap;		
	};
	
	/////////////////////////////
	// Declare static data member
	template <typename SpaceCellType,short unsigned int dim>
	std::map<Eigen::Matrix<int,dim,1>,SpaceCellType* const,CompareVectorsByComponent<int,dim> > SimplexObserver<SpaceCellType,dim,cellSize>::simplexMap;
	
}	// close namespace
#endif