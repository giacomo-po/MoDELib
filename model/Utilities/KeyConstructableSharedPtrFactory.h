/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef  model_KeyConstructableSharedPtrFactory_H_
#define model_KeyConstructableSharedPtrFactory_H_

#include <map>
#include <memory>

#include <CRTP.h>

namespace model
{

    template<typename Derived>
    struct KeyConstructableWeakPtrFactory : public CRTP<Derived>
    /*                                     */,public std::map<typename TypeTraits<Derived>::KeyType,
    /*                                                     */ const std::weak_ptr<typename TypeTraits<Derived>::ValueType>,
    /*                                                     */ typename TypeTraits<Derived>::CompareType>
    {
        
        typedef typename TypeTraits<Derived>::KeyType   KeyType;
        typedef typename TypeTraits<Derived>::ValueType ValueType;
        typedef std::weak_ptr<ValueType> WeakPtrType;
        typedef std::shared_ptr<ValueType> SharedPtrType;
        typedef std::map<KeyType,const SharedPtrType> MapType;
        typedef typename MapType::size_type SizeType;
        
        /**************************************************************************/
        SharedPtrType get(const KeyType& key)
        {
            const auto iter(this->find(key));
            if(iter!=this->end())
            {// key found
                if(iter->second.expired())
                {// shared_ptr was deleted elsewhere
                    this->erase(iter);
                    return this->emplace(key,SharedPtrType(new ValueType(this->derived(),key))).first->second.lock();
                }
                else
                {// shared_ptr exists
                    return iter->second.lock();
                }
            }
            else
            {// key not found found
//                return this->emplace(key,new ValueType(this->derived(),key)).first->second;
                return this->emplace(key,SharedPtrType(new ValueType(this->derived(),key))).first->second.lock();
            }
            
        }
        
    };
    
    template<typename Derived>
    struct KeyConstructableSharedPtrFactory : public CRTP<Derived>
    /*                                     */,public std::map<typename TypeTraits<Derived>::KeyType,
    /*                                                      */ const std::shared_ptr<typename TypeTraits<Derived>::ValueType>,
    /*                                                      */ typename TypeTraits<Derived>::CompareType>
    {

        typedef typename TypeTraits<Derived>::KeyType   KeyType;
        typedef typename TypeTraits<Derived>::ValueType ValueType;
        typedef std::shared_ptr<ValueType> SharedPtrType;
        typedef std::map<KeyType,const SharedPtrType> MapType;
        typedef typename MapType::size_type SizeType;

        /**************************************************************************/
        SharedPtrType get(const KeyType& key)
        {
            const auto iter(this->find(key));
            if(iter!=this->end())
            {// key found
                return iter->second;
            }
            else
            {
                return this->emplace(key,new ValueType(this->derived(),key)).first->second;
            }

        }


        //        /**************************************************************************/
        //        SharedPtrType get(const KeyType& key)
        //        {
        //            return MapType::emplace(std::piecewise_construct,
        //                                   std::forward_as_tuple(key),
        //                                   std::forward_as_tuple(this->derived(),key)
        //                                   ).first->second;
        //        }

        //        /**************************************************************************/
        //        SizeType erase(const SharedPtrType& shared) FINISH THIS PART
        //        {
        //            const auto iter(MapType::find(shared->key()));
        //            if(iter!=MapType::end())
        //            {
        //                if(iter->second->use_count()<=2)
        //                {// current map
        //
        //                }
        //            }
        //            else
        //            {
        //                return 0;
        //            }
        //        }

    };
		
}
#endif
