
#include <iostream>
#include <array>




int main (int argc, char * const argv[])
{
    
    std::array<int,3> a{1,5,4};
    std::array<int,3> b{1,2,4};
    
    std::cout<<(a<b)<<std::endl;
    
    return 0;
}

