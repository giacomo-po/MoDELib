#include <MyClassTemplate.h>
#include <iostream>

namespace internal
{
    template<int dim>
    void MyClassTemplate<dim>::print() const
    {
        std::cout<<"dim="<<dim<<std::endl;
    }

    template class MyClassTemplate<3>; // explicit instantiation
    
}


