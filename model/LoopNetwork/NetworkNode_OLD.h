/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_NetworkNode_H_
#define model_NetworkNode_H_

//#include <iostream>
//#include <list>
//#include <map>
#include <memory>
//#include <iterator>

#include <model/Utilities/StaticID.h>

namespace model
{
    template<typename Derived>
    class NetworkNode : public StaticID<Derived>
    {
        
    public:
        
                NetworkNode(const NetworkNode&) =delete;
//        const std::shared_ptr<LoopType>& pLoop;
    
//    NetworkNode()
//    {
////        std::cout<<"Constructing LoopNode "<<this->sID<<std::endl;
//    
//    }
    
    };
    
    
}
#endif
