/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef  model_SharedPtrFactories_h_
#define  model_SharedPtrFactories_h_

#include <map>
#include <memory>

#include <CRTP.h>

namespace model
{
    
    
    template<typename Derived,typename ValueType>
    struct KeyConstructableSharedPtrFactory : public CRTP<Derived>
    /*                                     */,public std::map<typename ValueType::KeyType,
    /*                                                      */ const std::shared_ptr<ValueType>,
    /*                                                      */ typename ValueType::CompareType>
    {

        typedef typename ValueType::KeyType   KeyType;
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
