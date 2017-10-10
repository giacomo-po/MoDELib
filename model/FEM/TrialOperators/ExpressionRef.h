/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef  model_ExpressionRef_h_
#define  model_ExpressionRef_h_

#include <memory> // std::unique_ptr
#include <typeinfo>  //for 'typeid' to work

namespace model
{
    
    //    namespace internal
    //    {
    template<typename T>
    class ExpressionRef
    {
        
        //            const T* const temp;
        //            std::unique_ptr<T> temp;
        const std::shared_ptr<const T> temp;
        
        const T& ref;
        
        
    public:
        
        ExpressionRef(const T& t):
        /* init */ temp(nullptr),
        /* init */ ref(t)
        {
            //                std::cout<<typeid(t).name()<<std::endl;
            //                std::cout<<"ExpressionRef from const&"<<std::endl;
        }
        
        ExpressionRef(T&& t):
        //            /* init */ temp(std::shared_ptr<const T>(std::move(t))),
        /* init */ temp(std::make_shared<T>(std::move(t))),
        //          /* init */ temp((std::move(in))),
        /* init */ ref(*temp)
        {
            //                std::cout<<typeid(*temp).name()<<std::endl;
            //                std::cout<<"ExpressionRef from &&"<<std::endl;
        }
        
        //            ExpressionRef(const ExpressionRef<T>& other)=delete;
        //            ExpressionRef(ExpressionRef<T>&& other)=default;
        
        const T& operator()() const
        {
            return ref;
        }
        
        //            T& operator()()
        //            {
        //                return ref;
        //            }
        
        //            T& operator()()
        //            {
        //                return ref;
        //            }
        
        //            ~ExpressionRef()
        //            {
        //                delete temp;
        //            }
        
    };
    //    }
}
#endif
