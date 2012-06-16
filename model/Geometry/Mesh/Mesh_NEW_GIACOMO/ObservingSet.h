// This file is part of mmdl, the C++ materials defect mechanics library.
//
// Copyright (C) 2011 by Giacomo Po <giacomopo@gmail.com>.
//
// mmdl is distributed without any warranty under the 
// GNU Lesser General Public License <http://www.gnu.org/licenses/>.


#ifndef mmdl_OBSERVINGSET_H_
#define mmdl_OBSERVINGSET_H_
#include <set>
#include <iterator>
#include <assert.h>		
#include <mmdl/Utilities/CRTP.h>

namespace mmdl {
	
	//////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////
	/*!  
	 *
	 * \brief ObservingSet is a class template that stores a static map of pointers to
	 * Derived class objects using the CRTP pattern. The key in the map is the staticID.
	 * The derived object gains access to the pointers of all the other objects of derived type that are instatiated. 
	 * The static map is automatically and dynamically updated during construction/destruction of derived objects. 
	 * A dynamic ID into the static vector is also provided. The typical implementation is:
	 *
	 * \code
	 class Derived : public ObservingSet<Derived>{
	 ...
	 };
	 *
	 * \endcode
	 *
	 * For multiple inheritance it's usefulto use the following pattern:
	 *
	 * \code
	 template <typename Derived2>
	 class Derived1 : public ObservingSet<Derived2>{
	 ...
	 };
	 *
	 class Derived2 : public Derived1<Derived2>{
	 ...
	 };
	 *
	 * \endcode
	 */
	
	
	/****************************************************************/
	/****************************************************************/
	template <typename Derived>
	class ObservedBySet;	// class predeclaration
	
	
	/****************************************************************/
	/****************************************************************/
	template <typename T>
	struct ObservingSet{
		
		typedef std::set<T*> SetType;	// std::set stores const objects so addresses will be const
		typedef typename SetType::const_iterator SetIteratorType;
		friend class ObservedBySet<T>;	// give ObservedBySet<T> access to staticSet

		/* size **********************************************/
		size_t  size() const {
			//! Returns the number of pointers stored in AddressSet
			return this->staticSet.size();
		}	
		
		/* observedBegin *************************************/
		SetIteratorType observedBegin() const {
			return this->staticSet.begin();
		}
		
		/* observedEnd ***************************************/
		SetIteratorType observedEnd() const {
			return this->staticSet.end();
		}
		
		/* observed ******************************************/
		T& observed(const unsigned int& k){
			assert(k<staticSet.size() && "REQUESTED ELEMENT OUT OF BOUND.");
			SetIteratorType it(staticSet.begin());
			std::advance (it,k);
			return **it;
		}
		
	private:
		static SetType staticSet; // the static map containing pointers to T objects

	};
	
	/////////////////////////////
	// Declare static data member
	template <typename T>
	std::set<T*> ObservingSet<T>::staticSet;
	
	
	
	/****************************************************************/
	/****************************************************************/
	template <typename Derived>
	struct ObservedBySet : public CRTP<Derived>,
	/*			       */  public ObservingSet<Derived>{
		
		/* Constructor **********************************************/
		ObservedBySet() {
			//! Inserts the pointer to the derived class in staticSet
			assert(this->staticSet.insert( this->p_derived() ).second && "COULD NOT INSERT OBSERVED TYPE POINTER INTO SET.");
		}
		
		/* Copy Constructor *****************************************/
		ObservedBySet(const ObservedBySet<Derived> & other){
			//! Inserts the pointer to the derived class in staticSet
			assert(this->staticSet.insert( this->p_derived() ).second && "COULD NOT INSERT OBSERVED TYPE POINTER INTO SET.");
		}
		
		/* Destructor ***********************************************/
		~ObservedBySet(){
			//! Removes the pointer to the derived class from staticSet
			assert(this->staticSet.erase( this->p_derived() ) && "COULD NOT ERASE OBSERVED TYPE POINTER FROM SET.");
		}
		
		/* dID ******************************************************/
		size_t dID() const {
			//! The dynamic ID of the derived object, i.e. the position in the ObservingSet
			return std::distance(this->staticSet.begin(),this->staticSet.find( this->p_derived() ) );
		}
		
	};
	
	
} // namespace mmdl
#endif

