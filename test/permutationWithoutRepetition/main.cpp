#include <iostream>
#include <Eigen/Dense>
#include <model/Math/CompileTimeMath/PermutationWithoutRepetition.h>

using namespace model;
int main()
{
    Eigen::Matrix<int,1,3> pool;
    pool<<1,2,3;
    
    std::cout<<PermutationWithoutRepetition<3>::permute(pool)<<std::endl<<std::endl;

    std::cout<<PermutationWithoutRepetition<3>::permuteWithPlusMinusSign(pool)<<std::endl;

    
    return 0;
}