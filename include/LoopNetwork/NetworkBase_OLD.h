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


namespace model
{
    
    template<typename NetworkType,typename _KeyType>
    class NetworkBase
    {

        NetworkType* const _network;
        
    public:
        
        typedef _KeyType KeyType;
        
        const KeyType key;
        NetworkBase(const NetworkBase<NetworkType,KeyType>& ) =delete;
        NetworkBase& operator=(const NetworkBase<NetworkType,KeyType>& other) =delete;

        
        /**********************************************************************/
        NetworkBase(NetworkType* const metwork_in,const KeyType& key_in) :
        /* init */ _network(metwork_in)
        /* init */,key(key_in)
        {
        }
        
        /**********************************************************************/
        ~NetworkBase()
        {
            size_t erased(_network->erase(key));
            assert(erased==1 && "Could not erase key");
        }
        
        
        /**********************************************************************/
        const NetworkType& network() const
        {
            return *_network;
        }
        
        NetworkType& network()
        {
            return *_network;
        }
        
        /**********************************************************************/
        NetworkType* p_network() const
        {
            return _network;
        }
        
//        NetworkType* const p_network() const
//        {
//            return _network;
//        }
        
        
    };
    
    
    
    
}
#endif
