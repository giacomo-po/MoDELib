#include <iostream>
#include <Eigen/Geometry>
#include <model/Math/CompileTimeMath/PermutationWithRepetition.h>
#include <model/Math/CompileTimeMath/PermutationWithOutRepetition.h>

using namespace model;
int main()
{
    Eigen::Matrix<int,1,3> pool;
    pool<<1,1,1;
    
    std::cout<<PermutationWithRepetition<3>::permute(pool)<<std::endl;
    
    return 0;
}
