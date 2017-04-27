/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LoopNode_H_
#define model_LoopNode_H_

#include <iostream>
//#include <list>
//#include <map>
#include <memory>
//#include <iterator>

#include <model/Utilities/StaticID.h>
#include <model/Utilities/CRTP.h>
//#include <iomanip>
//#include <vector>
//#include <Eigen/Dense>
//#include <Eigen/StdVector>
#include <model/LoopNetwork/NodeObserver.h>

namespace model
{
    template<typename Derived>
    class LoopNode : public StaticID<Derived>,
    /*            */ public CRTP<Derived>
    {
        
    public:
        
                LoopNode(const LoopNode&) =delete;
//        const std::shared_ptr<LoopType>& pLoop;
    
    LoopNode()
    {
        std::cout<<"Constructing LoopNode "<<this->sID<<std::endl;
        NodeObserver<Derived>::addNode(this->p_derived());
    }
        
        ~LoopNode()
        {
//            std::cout<<"Constructing LoopNode "<<this->sID<<std::endl;
            NodeObserver<Derived>::removeNode(this->p_derived());
        }
    
    };
    
    
}
#endif
