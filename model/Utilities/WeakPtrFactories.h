/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef  model_WeakPtrFactories_H_
#define  model_WeakPtrFactories_H_

#include <map>
#include <memory>

#include <CRTP.h>

namespace model
{
    
    template<typename Derived,typename ValueType>
    struct WeakPtrFactory : public std::map<typename ValueType::KeyType,
    /*                                   */ const std::weak_ptr<ValueType>>
    {
        
        typedef typename ValueType::KeyType   KeyType;
        typedef std::weak_ptr<ValueType> WeakPtrType;
        typedef std::shared_ptr<ValueType> SharedPtrType;
        
        
        template<typename...Args>
        SharedPtrType create(const Args&...args)
        {
            SharedPtrType temp(new ValueType(static_cast<Derived*>(this),args...));
            
            const auto iter(this->find(temp->key));
            if(iter!=this->end())
            {// key found
                if(iter->second.expired())
                {// shared_ptr was deleted elsewhere
                    this->erase(iter);
                    return this->emplace(temp->key,WeakPtrType(temp)).first->second.lock();
                }
                else
                {// shared_ptr exists
                    assert(false && "key already existing");
                    return SharedPtrType(nullptr);
                }
            }
            else
            {// key not found found
                return this->emplace(temp->key,WeakPtrType(temp)).first->second.lock();
            }
            
        }
        
        /**************************************************************************/
        SharedPtrType get(const KeyType& key) const
        {
            const auto iter(this->find(key));
            if(iter!=this->end())
            {// key found
                if(iter->second.expired())
                {// shared_ptr was deleted elsewhere
                    return SharedPtrType(nullptr);
                }
                else
                {// shared_ptr exists
                    return iter->second.lock();
                }
            }
            else
            {// key not found found
                return SharedPtrType(nullptr);
            }
        }
        
        /**************************************************************************/
        template<typename...Args>
        SharedPtrType clone(const KeyType& key,const Args&... args)
        {
            const auto iter(this->find(key));
            if(iter!=this->end())
            {// key found
                if(iter->second.expired())
                {// shared_ptr was deleted elsewhere
                    return nullptr;
                }
                else
                {// shared_ptr exists
                    const auto temp(iter->second.lock()->clone(args...));
                    return this->emplace(temp->key,WeakPtrType(temp)).first->second.lock();
                }
            }
            else
            {// key not found found
                return nullptr;
            }
        }
        
    };

    template<typename Derived,typename ValueType>
    struct KeyConstructableWeakPtrFactory : public std::map<typename  ValueType::KeyType,
    /*                                                   */ const std::weak_ptr<ValueType>>
    {
        typedef typename ValueType::KeyType   KeyType;
//        typedef typename TypeTraits<Derived>::ValueType ValueType;
        typedef std::weak_ptr<ValueType> WeakPtrType;
        typedef std::shared_ptr<ValueType> SharedPtrType;
        typedef std::map<KeyType,const SharedPtrType> MapType;
        typedef typename MapType::size_type SizeType;
        
        /**************************************************************************/
        template<typename...Args>
        SharedPtrType get(const Args&...args)
        {
            const KeyType key(ValueType::getKey(args...));
            const auto iter(this->find(key));
            if(iter!=this->end())
            {// key found
                if(iter->second.expired())
                {// shared_ptr was deleted elsewhere
                    this->erase(iter);
                    return this->emplace(key,SharedPtrType(new ValueType(static_cast<Derived*>(this),args...))).first->second.lock();
                }
                else
                {// shared_ptr exists
                    return iter->second.lock();
                }
            }
            else
            {// key not found found
//                return this->emplace(key,new ValueType(this->derived(),key)).first->second;
                return this->emplace(key,SharedPtrType(new ValueType(static_cast<Derived*>(this),args...))).first->second.lock();
            }
            
        }
        
    };
    
		
}
#endif
