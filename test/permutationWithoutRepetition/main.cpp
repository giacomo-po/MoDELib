#include <iostream>
#include <Eigen/Dense>
#include <model/Math/CompileTimeMath/PermutationWithoutRepetition.h>

using namespace model;
int main()
{
    Eigen::Matrix<int,1,3> pool;
    pool<<0,1,2;
    
    std::cout<<PermutationWithoutRepetition<2>::permute(pool)<<std::endl;
    
    return 0;
}