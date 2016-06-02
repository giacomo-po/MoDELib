#include <iostream>
#include <random>
#include <chrono>
#include <omp.h>
#include <deque>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/PardisoSupport>

int main()
{
    
    typedef Eigen::SparseMatrix<double> SparseMatType;
    typedef Eigen::Triplet<double> TripletType;
    
    // Generate random symmetric sparse matrix
    const int n=10000;  // matrix size
    const double fillRatio=0.2; // sparse fillRatio in [0 1]
    std::random_device rd; // random number generator
    std::uniform_int_distribution<int> dist(0, 100); // distribution for rd
    
    std::deque<TripletType> deqT;
    
    std::cout<<"Filling triplets"<<std::endl;
    for (int i=0;i<n;++i)
    {
        for(int j=i;j<n;++j)
        {
            const double v(dist(rd)/100.0); // random double in [0 1], also used as value in matrix
            if(v<fillRatio)
            {
                deqT.emplace_back(i,j,v);
                deqT.emplace_back(j,i,v); // make symmetric
            }
        }
    }
    
    std::cout<<"Creating matrix from triplets"<<std::endl;
    SparseMatType A(n,n);
    A.setFromTriplets(deqT.begin(),deqT.end());
    std::cout<<"non-zeros="<<A.nonZeros()<<std::endl;
    std::cout<<"fill ratio="<<((double)A.nonZeros())/n/n<<std::endl;

    std::cout<<"PardisoLDLT factorizing"<<std::endl;
    Eigen::PardisoLDLT<SparseMatType> solver(A);

    std::cout<<"PardisoLDLT solving"<<std::endl;
    Eigen::VectorXd b=Eigen::VectorXd::Random(n);
    Eigen::VectorXd x=solver.solve(b);
    
    return 0;
}