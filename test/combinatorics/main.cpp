#include <iostream>
#include <CTM.h>
#include <PermutationWithRepetition.h>
#include <PermutationWithOutRepetition.h>
#include <CombinationWithRepetition.h>

using namespace model;

template<int N,int k>
void printPermuation(const Eigen::Matrix<int,1,N>& pool)
{
    std::cout<<"### k="<<k<<" ###"<<std::endl;
    Eigen::MatrixXi c=CombinationWithRepetition<N,k>::combine(pool);
    for(int i=0;i<c.rows();++i)
    {
        Eigen::Matrix<int,1,k> row(c.row(i));
//        std::cout<<"##"<<row<<std::endl;
        Eigen::MatrixXi p=PermutationWithoutRepetition<k>::permute(row);
//        std::cout<<p<<std::endl;

        for(int j=0;j<p.cols();++j)
        {
            std::cout<<p.col(j).transpose()<<" ";
        }
        std::cout<<std::endl;
//        std::cout<<row<<std::endl;
    }
}

int main()
{
    Eigen::Matrix<int,1,3> pool;
    pool<<1,2,3;
    
    printPermuation<3,1>(pool);
    printPermuation<3,2>(pool);
    printPermuation<3,3>(pool);
    //std::cout<<PermutationWithRepetition<3>::permute(pool)<<std::endl;
    
//    std::cout<<"### k=1 ###"<<std::endl;
//    Eigen::MatrixXi c=CombinationWithRepetition<3,1>::combine(pool);
//    for(int r=0;r<c.rows();++r)
//    {
//        Eigen::MatrixXi p=PermutationWithRepetition<1>(c.row(r));
//    }
//    std::cout<<<<std::endl;
//
//    std::cout<<"### k=2 ###"<<std::endl;
//    std::cout<<CombinationWithRepetition<3,2>::combine(pool)<<std::endl;
//
//    std::cout<<"### k=3 ###"<<std::endl;
//    std::cout<<CombinationWithRepetition<3,3>::combine(pool)<<std::endl;
    
//    std::cout<<"### k=4 ###"<<std::endl;
//    std::cout<<CombinationWithRepetition<3,4>::combine(pool)<<std::endl;
//
//    std::cout<<"### k=5 ###"<<std::endl;
//    std::cout<<CombinationWithRepetition<3,5>::combine(pool)<<std::endl;
    
    return 0;
}
