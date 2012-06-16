/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_STATICID_H_
#define model_STATICID_H_

#include <assert.h> 

namespace model {
	
	/************************************************************/	
	/************************************************************/
	/*! \brief 
	 *  A class template that implements a counter of the number of instances of Derived type that are created at runtime.
	 *  It also provides a unique increasing static identifier for each instance.
	 *
	 * \code
	 * #include <model/Utilities/StaticID.h>
	 * 
	 * class MyClass : public model::StaticID<MyClass>{};
	 *
	 * int main(){
	 *
	 * MyClass a;		// with default constructor
	 * MyClass b(a);	// with copy constructor
	 *
	 * std::cout<<a.sID<<std::endl;
	 * std::cout<<b.sID<<std::endl;
	 * }
	 * \endcode
	 */
	
	template<typename Derived>
	class StaticID {

		// The increment
		static int increment;
		
		// The incremental counters
		static size_t count;
		static size_t nextCount;
		
	public:
		
		//! The static ID of this
		const  size_t sID;
		
		/* Default Constructor *****************************/
		StaticID() : sID(count) {
			count=nextCount;
			nextCount+=increment;
		}
		
		/* Copy Constructor ********************************/
		StaticID(const StaticID&) : sID(count) {
			count=nextCount;
			nextCount+=increment;
		}
		
		/* set_count  **************************************/

		static void set_count(const size_t& newCount){
			assert(newCount>=count && "YOU ARE TRYING TO SET THE COUNTER TO A LOWER VALUE THAN THE CURRENT ONE.");
			count =  newCount;
			nextCount=count+increment;
		}
		
		/* set_increment  **********************************/
		static void set_increment(const int& newIncrement){
			assert(newIncrement>=1 && "newIncrement MUST BE >=1.");
			nextCount+=(newIncrement-increment);	// now next time count will be used it will be correct
			increment = newIncrement;
		}
		
	};
	
	/* Declare static data members  *****************************/
	template<typename Derived>
	int StaticID<Derived>::increment = 1;

	template<typename Derived>
	size_t StaticID<Derived>::count = 0;

	template<typename Derived>
	size_t StaticID<Derived>::nextCount = 1;
	
	/************************************************************/	
	/************************************************************/
} // namespace model
#endif
