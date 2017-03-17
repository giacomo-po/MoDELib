/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_ARESAMETYPE_H_
#define model_ARESAMETYPE_H_


namespace model {
	
    template <typename A, typename B>
    struct AreSameType{
        enum{value=0};
    } __attribute__ ((deprecated));
    
    
    template <typename A>
    struct AreSameType<A,A>{
        enum{value=1};
    } __attribute__ ((deprecated));
	
	//////////////////////////////////////////////////////////////
} // namespace model
#endif


