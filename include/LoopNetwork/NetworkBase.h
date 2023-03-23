/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_NetworkBase_H_
#define model_NetworkBase_H_

#include <TerminalColors.h>
#include <map>
#include <memory>


namespace model
{
    
    template<typename Derived,typename _KeyType>
    class NetworkBase
    {

        typedef typename TypeTraits<Derived>::LoopNetworkType LoopNetworkType;
        typedef std::map<_KeyType,const std::weak_ptr<Derived>> WeakPtrContainerType;
        LoopNetworkType* const _network;
        WeakPtrContainerType* const weakPtrContainer;
        
    public:
        
        typedef _KeyType KeyType;
        
        const KeyType key;
        NetworkBase(const NetworkBase<LoopNetworkType,KeyType>& ) =delete;
        NetworkBase& operator=(const NetworkBase<LoopNetworkType,KeyType>& other) =delete;

        
        /**********************************************************************/
        NetworkBase(LoopNetworkType* const metwork_in,
                    WeakPtrContainerType* const wpc,
                    const KeyType& key_in) :
        /* init */ _network(metwork_in)
        /* init */,weakPtrContainer(wpc)
        /* init */,key(key_in)
        {
        }
        
        /**********************************************************************/
        ~NetworkBase()
        {
            weakPtrContainer->erase(key);
//            size_t erased(weakPtrContainer->erase(key));
//            assert(erased==1 && "Could not erase key");
        }
        
        
        /**********************************************************************/
        const LoopNetworkType& network() const
        {
            return *_network;
        }
        
        LoopNetworkType& network()
        {
            return *_network;
        }
        
        /**********************************************************************/
        LoopNetworkType* p_network() const
        {
            return _network;
        }
        
        size_t networkID() const
        {
            return std::distance(weakPtrContainer->begin(),weakPtrContainer->find(key));
        }

        
//        LoopNetworkType* const p_network() const
//        {
//            return _network;
//        }
        
        
    };
    
    
    
    
}
#endif
