#include <iostream>
#include <model/Math/CompileTimeMath/PermutationWithRepetition.h>

using namespace model;
int main()
{
    Eigen::Matrix<int,1,3> pool;
    pool<<0,1,2;
    
    std::cout<<PermutationWithRepetition<5>::permute(pool)<<std::endl;
    
    return 0;
}