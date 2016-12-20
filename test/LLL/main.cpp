#include <iostream>
#include <chrono>
#include <Eigen/Dense>
#include <model/LatticeMath/LLL.h>

using namespace model;



int main()
{

    Eigen::Matrix<int,3,4> B;
    B<< 1,-1, 3, 1,
    /**/1, 0, 5, 1,
    /**/1, 2, 6, 1;
    
    LLL lll(B);
    
      return 0;
}
