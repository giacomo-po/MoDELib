/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PeriodicDislocationLoop_H_
#define model_PeriodicDislocationLoop_H_


#include <memory>
#include <string>
#include <list>






namespace model
{

    template<typename DislocationNetworkType>
    class PeriodicDislocationNode : private std::map<Eigen::Matrix<double,dim,1>,std::shared_ptr<typename DislocationNetworkType::NodeType>,CompareVectorsByComponent<double,dim,float>>
    
    
    typedef std::map<Eigen::Matrix<double,dim,1>,std::shared_ptr<typename DislocationNetworkType::NodeType>,CompareVectorsByComponent<double,dim,float>> NodeContainerType;
    
    public:
    
    PeriodicDislocationNode()
    {
        
    }
    
    const NodeContainerType& nodes() const
    {
        return *this;
    }
    
    NodeContainerType& nodes()
    {
        return *this;
    }
    
}
#endif
