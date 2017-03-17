/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <map>
#include <Model/Utilities/NonCopyable.h>

using namespace model;

class A : public model::NonCopyable
{
    
    
public:
    
    const char i;
    const double d;
    

    
    
    A(const char& j, const double& dd) : i(j), d(dd)
    {
        std::cout<<"Creating object A "<<i<<std::endl;
    }
    
    ~A()
    {
        std::cout<<"Destroying object A "<<i<<std::endl;
    }
    
    
};



int main (int argc, char * const argv[])
{
    //A a('a');
    
    std::map<int,A> m;
    m.emplace(std::piecewise_construct, std::make_tuple(1), std::make_tuple('a', 2.0));
    
    return 0;
}


