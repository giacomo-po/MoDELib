// This file is part of mmdl, the C++ materials defect mechanics library.
//
// Copyright (C) 2011 by Giacomo Po <giacomopo@gmail.com>.
//
// mmdl is distributed without any warranty under the 
// GNU Lesser General Public License <http://www.gnu.org/licenses/>.


#ifndef mmdl_ADDRESSBOOK_H_
#define mmdl_ADDRESSBOOK_H_
#include <functional>	// for std::les
#include <map>
#include <iterator>		// for std::distance
#include <assert.h>		// for std::distance
#include <boost/utility.hpp>	// for noncopyable
#include <mmdl/Utilities/CRTP.h>
#include <mmdl/Utilities/StaticID.h>
#include <mmdl/Utilities/CompareVectorsByComponent.h>

namespace mmdl {
	
	
	template <typename KeyType, typename T, typename Compare>
	class ObservedByMap;
	
	
	//////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////
	template <typename KeyType, typename T, typename Compare>
	class ObservingMap{
		
	public:
		typedef std::map<KeyType,T* const,Compare> MapType;
		typedef typename MapType::iterator staticMapIteratorType;
		
	private:
		static MapType staticMap;		// the static map containing pointers to T objects
		friend class ObservedByMap<KeyType,T,Compare>;	// give ObservedByMap<T> access to staticMap
		
	public:
		//////////////////////////////////////////////////////////////
		// Naddresses
		size_t  size() const {
			//! Returns the number of pointers stored in AddressSet
			return this->staticMap.size();
		}	
		
		//////////////////////////////////////////////////////////////
		// address
		T* observed(const KeyType& key) const {
			// element type of std::map are of type pair<const Key,T>
			// return type of find is an iterator iter to a pair<const Key,T>
			// Dereferencing this iterator accesses the element's value, which is of type pair<const Key,T>
			return this->staticMap.find(key)->second;
		}
		
		staticMapIteratorType observedBegin() const {
			return this->staticMap.begin();
		}
		
		
		staticMapIteratorType observedEnd() const {
			return this->staticMap.end();
		}
	};
	
	/////////////////////////////
	// Declare static data member
	template <typename KeyType, typename T, typename Compare>
	typename ObservingMap<KeyType,T,Compare>::MapType ObservingMap<KeyType,T,Compare>::staticMap;
	
	
	
	

	
	
	
	//////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////
	// Template Specialization ObservingMap<Derived,1> (active)
//	template <typename Derived>
	template <typename KeyType, typename Derived, typename Compare>
	class ObservedByMap : boost::noncopyable,	// uniqueID cannot be copyed
	/*                 */ public CRTP<Derived>,
	/*			       */ public ObservingMap<KeyType,Derived,Compare>{
		
//		typedef Derived::KeyType KeyType;
//		typedef ObservingMapBase<Derived>::MapType MapType;
//		typedef typename MapType::iterator staticMapIteratorType;
		


		
		public:
				
		const KeyType uniqueID;
		
		
		//////////////////////////////////////////////////////////////
		// Constructor
		ObservedByMap(const KeyType& uniqueID_in) : uniqueID(uniqueID_in){
			//! Inserts the pointer to the derived class in staticMap
			assert(this->staticMap.insert(std::make_pair( uniqueID, this->p_derived() ) ).second && "COULD NOT INSERT OBSERVED TYPE POINTER INTO MAP.");
		}
		

				
		
		//////////////////////////////////////////////////////////////
		// Destructor
		~ObservedByMap(){
			//! Removes the pointer to the derived class from staticMap
			assert(this->staticMap.erase( uniqueID ) && "COULD NOT ERASE OBSERVED TYPE POINTER FROM MAP.");
		}
		
		//////////////////////////////////////////////////////////////
		// dID
		size_t dID() const {
			//! Returns the dynamic identifier of the derived object.
			return std::distance(this->staticMap.begin(),this->staticMap.find( uniqueID ) );
		}
		
		
		//		//////////////////////////////////////////////////////////////
		//		// Copy Constructor
		//		ObservingMap(const ObservingMap<Derived> & other){
		//			//! Inserts the pointer to the derived class in staticMap
		//			assert(this->staticMap.insert(std::make_pair( this->sID, this->p_derived() ) ).second && "COULD NOT INSERT OBSERVED TYPE POINTER INTO MAP.");
		//		}

		
	};
	
	

	
} // namespace mmdl
#endif
