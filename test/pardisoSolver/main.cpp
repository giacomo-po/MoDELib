#include <iostream>
#include <random>
#include <chrono>
#include <omp.h>
#include <deque>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/PardisoSupport>


template<typename SolverType>
void testSolver(const Eigen::SparseMatrix<double>& A,
                const Eigen::VectorXd& b)
{
    std::cout<<"    factorizing"<<std::flush;
    const auto t0= std::chrono::system_clock::now();
    SolverType solver(A);
    std::cout<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<std::endl;
    std::cout<<"    factorizing succesful? "<< (solver.info()==Eigen::Success)<<std::endl;
    
    std::cout<<"    solving"<<std::flush;
    const auto t1= std::chrono::system_clock::now();
    Eigen::VectorXd x=solver.solve(b);
    std::cout<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t1)).count()<<" sec]"<<std::endl;
    std::cout<<"    relative error: "<<(A*x-b).norm()/b.norm()<<std::endl;
}


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
    const Eigen::VectorXd b=Eigen::VectorXd::Random(n); // lhs
    
    
    std::cout<<"Filling triplets..."<<std::flush;
    const auto t0= std::chrono::system_clock::now();
    for (int i=0;i<n;++i)
    {
        for(int j=i;j<n;++j)
        {
            if(j==i)
            {
                deqT.emplace_back(i,j,n); // force SPD
            }
            else
            {
                double v(dist(rd)/100.0); // uniform random double in [0 1], also used as value in matrix
                if(v<fillRatio)
                {
                    deqT.emplace_back(i,j,v);
                    deqT.emplace_back(j,i,v); // make symmetric
                }
            }
            
        }
    }
    std::cout<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<std::endl;
    
    
    std::cout<<"Creating matrix from triplets..."<<std::flush;
    const auto t1= std::chrono::system_clock::now();
    SparseMatType A(n,n);
    A.setFromTriplets(deqT.begin(),deqT.end());
    std::cout<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t1)).count()<<" sec]"<<std::endl;
    std::cout<<"matrix size="<<n<<"x"<<n<<std::endl;
    std::cout<<"non-zeros="<<A.nonZeros()<<std::endl;
    std::cout<<"fill ratio="<<((double)A.nonZeros())/n/n<<std::endl;
    
    // Testing various solvers
    std::cout<<"PardisoLLT"<<std::endl;
    testSolver<Eigen::PardisoLLT<SparseMatType>>(A,b);
    
    std::cout<<"PardisoLDLT"<<std::endl;
    testSolver<Eigen::PardisoLDLT<SparseMatType>>(A,b);

    std::cout<<"SimplicialLLT"<<std::endl;
    testSolver<Eigen::SimplicialLLT<SparseMatType>>(A,b);
    
    std::cout<<"SimplicialLDLT"<<std::endl;
    testSolver<Eigen::SimplicialLDLT<SparseMatType>>(A,b);
    
    std::cout<<"SparseLU"<<std::endl;
    testSolver<Eigen::SparseLU<SparseMatType>>(A,b);
    

    
    return 0;
}