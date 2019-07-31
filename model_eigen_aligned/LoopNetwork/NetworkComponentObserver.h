/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * mmdl is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_NetworkComponentObserver_H_
#define model_NetworkComponentObserver_H_

#include <iostream>
#include <map>
#include <iterator>		// for std::distance
#include <CRTP.h>
#include <StaticID.h>

namespace model
{

    
    /*!\brief NetworkComponentObserver is a class template that stores a static map of pointers to
     * Derived class objects using the CRTP pattern. The key in the map is the staticID.
     * The derived object gains access to the pointers of all the other objects of derived type that are instatiated.
     * The static map is automatically and dynamically updated during construction/destruction of derived objects.
     * A dynamic ID into the static vector is also provided. The typical implementation is:
     *
     * \code
     class Derived : public NetworkComponentObserver<Derived>{
     ...
     };
     *
     * \endcode
     *
     * For multiple inheritance it's usefulto use the following pattern:
     *
     * \code
     template <typename Derived2>
     class Derived1 : public NetworkComponentObserver<Derived2>{
     ...
     };
     *
     class Derived2 : public Derived1<Derived2>{
     ...
     };
     *
     * \endcode
     */
	template <typename NetworkComponentType>
	class NetworkComponentObserver
    {
		
    public:
		typedef std::map<size_t,NetworkComponentType* const> NetworkComponentContainerType;
		
	private:
		static NetworkComponentContainerType networkComponentMap;
		
	public:
        
        static void addComponent(NetworkComponentType* const pNC)
        {
            const bool success=networkComponentMap.insert(std::make_pair(pNC->sID,pNC)).second;
            assert(success && "COULD NOT INSERT NetworkComponent");
        }
        
        static void removeComponent(NetworkComponentType* const pNC)
        {
            const int erased=networkComponentMap.erase(pNC->sID);
            assert(erased==1 && "COULD NOT ERASE NetworkComponent");
        }
		
        static const NetworkComponentContainerType& components()
        {
            return networkComponentMap;
        }
        
	};
	
	/////////////////////////////
	// Declare static data member
	template <typename NetworkComponentType>
	std::map<size_t,NetworkComponentType* const> NetworkComponentObserver<NetworkComponentType>::networkComponentMap;
	
	
	
	
//
//	
//	template <typename Derived, bool active = 1>
//	class NetworkComponentObserver{};
//	
//	
//	
//	//////////////////////////////////////////////////////////////
//	//////////////////////////////////////////////////////////////
//	// Template Specialization NetworkComponentObserver<Derived,1> (active)
//	template <typename Derived>
//	class NetworkComponentObserver<Derived,1> :	private NetworkComponentObserverBase<Derived>,
//	/*							*/  public  StaticID<Derived>,
//	/*							*/  public  CRTP<Derived>{
//		
//		typedef std::map<size_t,Derived* const> AddressMapType;
//		typedef typename AddressMapType::iterator AddressMapIteratorType;
//		
//
//		
//		
//		public:
//				
//		AddressMapIteratorType ABbegin() const {
//			return this->AddressMap.begin();
//		}
//
//		
//		AddressMapIteratorType ABend() const {
//			return this->AddressMap.end();
//		}
//		
//		
//		//////////////////////////////////////////////////////////////
//		// Constructor
//		NetworkComponentObserver(){
//			//! 1- Insert the pointer to the derived class in AddressMap
//			this->AddressMap.insert(std::make_pair( this->sID, this->p_derived() ) );
//			//newInstance();
//		}
//		
//		//////////////////////////////////////////////////////////////
//		// Copy Constructor
//		NetworkComponentObserver(const NetworkComponentObserver<Derived> & I){
//			//! 1- Insert the pointer to the derived class in AddressMap
//			this->AddressMap.insert(std::make_pair( this->sID, this->p_derived() ) );
//			//newInstance();
//		}
//				
//		
//		//////////////////////////////////////////////////////////////
//		// Destructor
//		~NetworkComponentObserver(){
//			this->AddressMap.erase( this->sID );
//		}
//		
//		//////////////////////////////////////////////////////////////
//		// dID
//		size_t dID() const {
//			//! Returns the dynamic identifier of the derived object.
//			return std::distance(this->AddressMap.begin(),this->AddressMap.find( this->sID ) );
//		}
//		
//		//////////////////////////////////////////////////////////////
//		// Naddresses
//		size_t  Naddresses() const {
//			//! Returns the number of pointers stored in AddressSet
//			return this->AddressMap.size();
//		}	
//		
//		//////////////////////////////////////////////////////////////
//		// address
//		Derived* address(const size_t& k) const {
//			// element type of std::map are of type pair<const Key,T>
//			// return type of find is an iterator iter to a pair<const Key,T>
//			// Dereferencing this iterator accesses the element's value, which is of type pair<const Key,T>
//			return this->AddressMap.find(k)->second;
//		}
//		
//	};
//	
//	
//	//////////////////////////////////////////////////////////////
//	//////////////////////////////////////////////////////////////
//	// Template Specialization NetworkComponentObserver<Derived,0> (pasive)
//	template <typename T>
//	class NetworkComponentObserver<T,0> : private NetworkComponentObserverBase<T>{
//		
//		typedef std::map<size_t,T* const> AddressMapType;
//		typedef typename AddressMapType::iterator AddressMapIteratorType;
//				
//	public:
//		//////////////////////////////////////////////////////////////
//		// Naddresses
//		size_t  Naddresses() const {
//			//! Returns the number of pointers stored in AddressSet
//			return this->AddressMap.size();
//		}	
//		
//		//////////////////////////////////////////////////////////////
//		// address
//		T* address(const size_t& k) const {
//			// element type of std::map are of type pair<const Key,T>
//			// return type of find is an iterator iter to a pair<const Key,T>
//			// Dereferencing this iterator accesses the element's value, which is of type pair<const Key,T>
//			return this->AddressMap.find(k)->second;
//		}
//		
//		AddressMapIteratorType ABbegin() const {
//			return this->AddressMap.begin();
//		}
//		
//		
//		AddressMapIteratorType ABend() const {
//			return this->AddressMap.end();
//		}
//		
//	};
	
} // namespace model
#endif

