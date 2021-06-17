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
#include <assert.h>

#include <CRTP.h>

namespace model
{
    
    
    template<typename Derived,typename ValueType,typename CompareType=std::less<typename  ValueType::KeyType>>
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

    template<typename Derived,typename ValueType,typename CompareType=std::less<typename  ValueType::KeyType>>
    struct KeyConstructableWeakPtrFactory : public std::map<typename  ValueType::KeyType,
    /*                                                   */ const std::weak_ptr<ValueType>,
    /*                                                   */ CompareType>
    {
        typedef typename ValueType::KeyType   KeyType;
        typedef std::weak_ptr<ValueType> WeakPtrType;
        typedef std::shared_ptr<ValueType> SharedPtrType;
        
        /**************************************************************************/
        SharedPtrType getFromKey(const KeyType& key)
        {
            const auto iter(this->find(key));
            if(iter!=this->end())
            {// key found
                if(iter->second.expired())
                {// shared_ptr was deleted elsewhere
                    this->erase(iter);
                    return this->emplace(key,SharedPtrType(new ValueType(static_cast<Derived*>(this),key))).first->second.lock();
                }
                else
                {// shared_ptr exists
                    return iter->second.lock();
                }
            }
            else
            {// key not found found
                //                return this->emplace(key,new ValueType(this->derived(),key)).first->second;
                return this->emplace(key,SharedPtrType(new ValueType(static_cast<Derived*>(this),key))).first->second.lock();
            }
            
        }
        
//        size_t keyID(const KeyType& key) const
//        {
//            return std::distance(this->find(key),this->begin());
//        }
        
        
        //        /**************************************************************************/
        //        template<typename...Args>
        //        SharedPtrType get(const Args&...args)
        //        {
        //            const KeyType key(ValueType::getKey(args...));
        //            const auto iter(this->find(key));
        //            if(iter!=this->end())
        //            {// key found
        //                if(iter->second.expired())
        //                {// shared_ptr was deleted elsewhere
        //                    this->erase(iter);
        //                    return this->emplace(key,SharedPtrType(new ValueType(static_cast<Derived*>(this),args...))).first->second.lock();
        //                }
        //                else
        //                {// shared_ptr exists
        //                    return iter->second.lock();
        //                }
        //            }
        //            else
        //            {// key not found found
        ////                return this->emplace(key,new ValueType(this->derived(),key)).first->second;
        //                return this->emplace(key,SharedPtrType(new ValueType(static_cast<Derived*>(this),args...))).first->second.lock();
        //            }
        //
        //        }
        
    };
    
		
}
#endif


//
//        struct iterator : public BaseType::iterator
//        {
//
//            typedef typename BaseType::iterator BaseIteratorType;
//
//            const BaseType* const p_map;
//
//
//            iterator(const BaseType* const p_map_in, const BaseIteratorType& iter) :
//            /* init */ BaseType::iterator(iter)
//            /* init */,p_map(p_map_in)
//            {
//
//            }
//
//            iterator& operator++()
//            {
//                iterBase().operator++();
//                while(iterBase()!=p_map->end())
//                {
//                    if(iterBase()->second.expired())
//                    {
//                        iterBase().operator++();
//                    }
//                    else
//                    {
//                        break;
//                    }
//
//                }
//                return *this;
//            }
//
//            bool operator==(const iterator& other) const
//            {
//                return iterBase()==other.iterBase();
//            }
//
//            bool operator!=(const iterator& other) const
//            {
//                return iterBase()!=other.iterBase();
//            }
//
//            BaseIteratorType& iterBase()
//            {
//                return *this;
//            }
//
//            const BaseIteratorType& iterBase() const
//            {
//                return *this;
//            }
//
//        };
//
//
//        struct const_iterator : public BaseType::const_iterator
//        {
//
//            typedef typename BaseType::const_iterator BaseIteratorType;
//
//            const BaseType* const p_map;
//
//
//            const_iterator(const BaseType* const p_map_in, const BaseIteratorType& iter) :
//            /* init */ BaseType::const_iterator(iter)
//            /* init */,p_map(p_map_in)
//            {
//
//            }
//
//            const_iterator& operator++()
//            {
//                iterBase().operator++();
//                while(iterBase()!=p_map->end())
//                {
//                    if(iterBase()->second.expired())
//                    {
//                        iterBase().operator++();
//                    }
//                    else
//                    {
//                        break;
//                    }
//
//                }
//                return *this;
//            }
//
//            bool operator==(const const_iterator& other) const
//            {
//                return iterBase()==other.iterBase();
//            }
//
//            bool operator!=(const const_iterator& other) const
//            {
//                return iterBase()!=other.iterBase();
//            }
//
//            BaseIteratorType& iterBase()
//            {
//                return *this;
//            }
//
//            const BaseIteratorType& iterBase() const
//            {
//                return *this;
//            }
//
//        };
//
//        iterator begin()
//        {
//            auto baseIter(base().begin());
//            while(baseIter!=base().end())
//            {
//                if(baseIter->second.expired())
//                {
//                    baseIter++;
//                }
//                else
//                {
//                    break;
//                }
//            }
//            return iterator(this,baseIter);
//        }
//
//        const_iterator begin() const
//        {
//            auto baseIter(base().begin());
//            while(baseIter!=base().end())
//            {
//                if(baseIter->second.expired())
//                {
//                    baseIter++;
//                }
//                else
//                {
//                    break;
//                }
//            }
//            return const_iterator(this,baseIter);
//        }
//
//        iterator end()
//        {
//            return iterator(this,base().end());
//        }
//
//        const_iterator end() const
//        {
//            return const_iterator(this,base().end());
//        }
//
//        BaseType& base()
//        {
//            return *this;
//        }
//
//        const BaseType& base() const
//        {
//            return *this;
//        }


//    template<typename ValueType,typename CompareType>
//    struct WeakPtrFactoryBase : public std::map<typename ValueType::KeyType,
//    /*                                       */ const std::weak_ptr<ValueType>,
//    /*                                       */ CompareType>
//    {
//
//        typedef typename ValueType::KeyType   KeyType;
//        typedef std::weak_ptr<ValueType> WeakPtrType;
//        typedef std::shared_ptr<ValueType> SharedPtrType;
//        typedef std::map<KeyType, const WeakPtrType> BaseType;
//
//
//        class iterator : public BaseType::iterator
//        {
//
//            typedef typename BaseType::iterator BaseIteratorType;
//
//            const BaseType* const p_map;
//
//        public:
//
//            iterator(const BaseType* const p_map_in, const BaseIteratorType& iter) :
//            /* init */ BaseType::iterator(iter)
//            /* init */,p_map(p_map_in)
//            {
//
//            }
//
//            iterator& operator++()
//            {
//                iterBase().operator++();
//                while(iterBase()!=p_map->end())
//                {
//                    if(iterBase()->second.expired())
//                    {
//                        iterBase().operator++();
//                    }
//                    else
//                    {
//                        break;
//                    }
//
//                }
//                return *this;
//            }
//
//            bool operator==(const iterator& other) const
//            {
//                return iterBase()==other.iterBase();
//            }
//
//            bool operator!=(const iterator& other) const
//            {
//                return iterBase()!=other.iterBase();
//            }
//
//            BaseIteratorType& iterBase()
//            {
//                return *this;
//            }
//
//            const BaseIteratorType& iterBase() const
//            {
//                return *this;
//            }
//
//        };
//
//        class const_iterator : public BaseType::const_iterator
//        {
//
//            typedef typename BaseType::const_iterator BaseIteratorType;
//
//            const BaseType* const p_map;
//
//            public:
//
//            const_iterator(const BaseType* const p_map_in, const BaseIteratorType& iter) :
//            /* init */ BaseType::const_iterator(iter)
//            /* init */,p_map(p_map_in)
//            {
//
//            }
//
//            const_iterator& operator++()
//            {
//                iterBase().operator++();
//                while(iterBase()!=p_map->end())
//                {
//                    if(iterBase()->second.expired())
//                    {
//                        iterBase().operator++();
//                    }
//                    else
//                    {
//                        break;
//                    }
//
//                }
//                return *this;
//            }
//
//            bool operator==(const const_iterator& other) const
//            {
//                return iterBase()==other.iterBase();
//            }
//
//            bool operator!=(const const_iterator& other) const
//            {
//                return iterBase()!=other.iterBase();
//            }
//
//            BaseIteratorType& iterBase()
//            {
//                return *this;
//            }
//
//            const BaseIteratorType& iterBase() const
//            {
//                return *this;
//            }
//
//        };
//
//        iterator begin()
//        {
//            auto baseIter(base().begin());
//            while(baseIter!=base().end())
//            {
//                if(baseIter->second.expired())
//                {
//                    baseIter++;
//                }
//                else
//                {
//                    break;
//                }
//            }
//            return iterator(this,baseIter);
//        }
//
//        const_iterator begin() const
//        {
//            auto baseIter(base().begin());
//            while(baseIter!=base().end())
//            {
//                if(baseIter->second.expired())
//                {
//                    baseIter++;
//                }
//                else
//                {
//                    break;
//                }
//            }
//            return const_iterator(this,baseIter);
//        }
//
//        iterator end()
//        {
//            return iterator(this,base().end());
//        }
//
//        const_iterator end() const
//        {
//            return const_iterator(this,base().end());
//        }
//
//        BaseType& base()
//        {
//            return *this;
//        }
//
//        const BaseType& base() const
//        {
//            return *this;
//        }
//
//        size_t keyID(const KeyType& key) const
//        {
//            size_t temp(0);
//            for(const_iterator pair=begin();pair!=end();++pair)
//            {
//                if(pair->first==key)
//                {
//                    return temp;
//                }
//                temp++;
//            }
//            assert(false && "key nof found");
//            return temp;
//        }
//    };
