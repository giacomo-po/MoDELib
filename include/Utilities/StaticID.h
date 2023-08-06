/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_STATICID_H_
#define model_STATICID_H_


namespace model
{
	
	/************************************************************/	
	/************************************************************/
	/*! \brief A class template that implements a counter of the number of 
     * instances of Derived type that are created at runtime. It also provides a 
     * unique increasing static identifier (sID) for each instance.
	 *
     * Example:
     * \include test/test_StaticID/main.cpp 
     * Output:
     * \include test/test_StaticID/output.txt
	 */
	template<typename Derived>
	class StaticID
    {

		// The increment
		static size_t increment;
		
		// The incremental counters
		static size_t count;
        static bool count_used;
		
	public:
		
		//! The static ID of this
		const  size_t sID;
		
        /**********************************************************************/
		StaticID() : sID(count)
        {
            count_used=true;
            count+=increment;
		}
		
        /**********************************************************************/
		StaticID(const StaticID&) : sID(count)
        {
            count_used=true;
            count+=increment;
		}
        
        /**********************************************************************/
        static size_t nextID()
        {
            return count;
        }
        
        static size_t& get_count()
        {
            return count;
        }
		
        /**********************************************************************/
		static void set_count(const size_t& newCount)
        {
//			model_checkInput(newCount>=count && "YOU ARE TRYING TO SET THE COUNTER TO A LOWER VALUE THAN THE CURRENT ONE.");
            assert(newCount>=count && "YOU ARE TRYING TO SET THE COUNTER TO A LOWER VALUE THAN THE CURRENT ONE.");
			count =  newCount;
            count_used=false;
		}
        
        static void force_count(const size_t& newCount)
        {
            count =  newCount;
            count_used=false;
        }
		
        /**********************************************************************/
		static void set_increment(const size_t& newIncrement)
        {
//			model_checkInput(newIncrement>=1 && "newIncrement MUST BE >=1.");
            assert(newIncrement>=1 && "newIncrement MUST BE >=1.");
            if(count_used)
            {
                count-=increment;
                count+=newIncrement;
            }
            increment = newIncrement;
		}
		
	};
	
	/* Static data members  *****************************/
	template<typename Derived>
	size_t StaticID<Derived>::increment = 1;

	template<typename Derived>
	size_t StaticID<Derived>::count = 0;

    template<typename Derived>
    bool StaticID<Derived>::count_used = false;
	
} // namespace model
#endif
