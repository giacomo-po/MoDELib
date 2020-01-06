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



#include <StaticID.h>
#include <PeriodicLoopObserver.h>

namespace model
{

//    template<typename DislocationNetworkType>
//    class PeriodicDislocationNode : private std::map<Eigen::Matrix<double,dim,1>,std::shared_ptr<typename DislocationNetworkType::NodeType>,CompareVectorsByComponent<double,dim,float>>
//
//
//    typedef std::map<Eigen::Matrix<double,dim,1>,std::shared_ptr<typename DislocationNetworkType::NodeType>,CompareVectorsByComponent<double,dim,float>> NodeContainerType;
//
//    public:
//
//    PeriodicDislocationNode()
//    {
//
//    }
//
//    const NodeContainerType& nodes() const
//    {
//        return *this;
//    }
//
//    NodeContainerType& nodes()
//    {
//        return *this;
//    }
    
    template<typename DislocationNetworkType>
    class PeriodicDislocationLoop : public StaticID<PeriodicDislocationLoop<DislocationNetworkType>>
    /*                           */,public std::map<size_t,typename DislocationNetworkType::LoopType*>
    {
        
        typedef typename DislocationNetworkType::LoopType LoopType;
        typedef std::map<size_t,LoopType*> LoopContainerType;
        typedef PeriodicDislocationLoop<DislocationNetworkType> PeriodicDislocationLoopType;
        typedef PeriodicLoopObserver<PeriodicDislocationLoopType> PeriodicDislocationLoopObserverType;

        // PeriodicGlidePlane<dim> periodicPlane;
        
        PeriodicLoopObserver<PeriodicDislocationLoopType>* const observer;
        
        public:
        
        PeriodicDislocationLoop(PeriodicDislocationLoopObserverType* const obs) :
        /* init */ observer(obs)
//        periodicPlane(...)
        {
            observer->addPeriodicLoop(this);
        }
        
        ~PeriodicDislocationLoop()
        {
            observer->removePeriodicLoop(this);
        }
        
        const LoopContainerType& loops() const
        {
            return *this;
        }
        
        //        LoopContainerType& loops()
        //        {
        //            return *this;
        //        }
        
        void addLoop(const LoopType* )
        {
            // add to LoopContainerType
            // update outLoop reconstrcuction
//            periodicPlane.addPatch(loop->periodicShift);
        }
        
        void removeLoop(const LoopType* )
        {
            // femove from LoopContainerType
            // update outLoop reconstrcuction
        }
        

        
        
    };

    
}
#endif
