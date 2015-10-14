#include <iostream>
#include <model/Utilities/StaticID.h>

class MyClass : public model::StaticID<MyClass>{};

int main()
{
    
    MyClass a;		// with default constructor
    MyClass b(a);	// with copy constructor
    MyClass c=b;	// with assigment operator
    
    std::cout<<a.sID<<std::endl;
    std::cout<<b.sID<<std::endl;
    std::cout<<c.sID<<std::endl;
    
    return 0;
}