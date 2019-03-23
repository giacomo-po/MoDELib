/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef  model_NonCopyable_H_
#define  model_NonCopyable_H_

namespace model
{
    /*! \brief A non-copyable base class using c++11.
     */
    class NonCopyable
    {
    protected:
        NonCopyable() = default;   //default constructor
        ~NonCopyable() = default;  //default destructor
        
//        NonCopyable(NonCopyable &&) = default; //default move constructor
//        NonCopyable& operator=(NonCopyable &&) = default; //default move assignment
        
    public:
        NonCopyable(const NonCopyable&) = delete; // non construction-copyable
        const NonCopyable& operator=(const NonCopyable&) = delete; // non assignable
    };
    
}
#endif
