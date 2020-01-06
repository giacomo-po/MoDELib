
#ifndef model_DDcuda_h_
#define model_DDcuda_h_

#include <iostream>
#include <vector>

namespace model
{
    class DDcuda
    {
    
    
        
    public:
        
        DDcuda();
        
        template <typename T1,typename T2>
        void computeStressAtQuadrature(const std::vector<T1>&,
                                     std::vector<T2>&);
        
        
    };
}
#endif
