#include <iostream>
#include <model/Utilities/StaticID.h>

class MyClass : public model::StaticID<MyClass>{};

int main()
{
    
    MyClass a;		// with default constructor
    MyClass b(a);	// with copy constructor
    MyClass c=b;	// with assigment operator

    model::StaticID<MyClass>::set_increment(3);
    MyClass d;

    model::StaticID<MyClass>::set_count(10);
    MyClass e;
    MyClass f;

    
    std::cout<<a.sID<<std::endl;
    std::cout<<b.sID<<std::endl;
    std::cout<<c.sID<<std::endl;
    std::cout<<d.sID<<std::endl;
    std::cout<<e.sID<<std::endl;
    std::cout<<f.sID<<std::endl;
    
    return 0;
}