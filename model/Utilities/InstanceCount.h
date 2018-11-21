/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_InstanceCount_H_
#define model_InstanceCount_H_

//#include <modelMacros.h> // model_checkInput

namespace model
{
	
	/************************************************************/	
	/************************************************************/
	/*! \brief A class template that implements a counter of the number of 
     * instances of Derived type that are created at runtime. It also provides a 
     * unique increasing static identifier (sID) for each instance.
	 *
     * Example:
     * \include test/test_InstanceCount/main.cpp 
     * Output:
     * \include test/test_InstanceCount/output.txt
	 */
	template<typename Derived>
	class InstanceCount
    {

		
		// The incremental counters
		static size_t count;
		
	public:
		

		
        /**********************************************************************/
		InstanceCount()
        {
            count+=1;
		}
		
        /**********************************************************************/
		InstanceCount(const InstanceCount&) 
        {
            count+=1;
		}
        
        /**********************************************************************/
        ~InstanceCount()
        {
            count-=1;
        }
        
        /**********************************************************************/
        static size_t& counter()
        {
            return count;
        }
		

		

		
	};
	
	template<typename Derived>
	size_t InstanceCount<Derived>::count = 0;
	
} // namespace model
#endif
