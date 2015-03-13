//g++ main.cpp -o test -O3 -std=c++11 -fopenmp -I/usr/local/include/ -I.


#include <iostream>
#include <chrono>
#include <omp.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/PardisoSupport>

//#include <Eigen/SparseQR>
//#include <Eigen/OrderingMethods>
//
//#include <model/Math/SparseNullSpace.h>

// http://stackoverflow.com/questions/2181418/computing-the-null-space-of-a-matrix-as-fast-as-possible

using namespace model;


int main()
{

    typedef Eigen::SparseMatrix<double> SparseMatType; // declares a column-major sparse matrix type of double
    typedef Eigen::Triplet<double> TripletType;


    int r=4;
    int c=6;
    
    SparseMatType A(r,c);
    A.insert(0,0) =  1.0;
    A.insert(0,2) = -3.0;
    A.insert(0,4) =  2.0;
    A.insert(0,5) = -8.0;
    A.insert(1,1) =  1.0;
    A.insert(1,2) =  5.0;
    A.insert(1,4) = -1.0;
    A.insert(1,5) =  4.0;
    A.insert(2,3) =  1.0;
    A.insert(2,4) =  7.0;
    A.insert(2,5) = -9.0;
    A.makeCompressed();                        // optional

//    SparseNullSpace<SparseMatType> ns(A);
    

//    SparseMatType AZ;
//    AZ=A*ns.matrixZ();

//    std::cout<<"(A*Z).norm()="<<(A*ns.matrixZ()).norm()<<std::endl;
//    std::cout<<"(A*Y).norm()="<<(A*ns.matrixY()).norm()<<std::endl;

    Eigen::PardisoDLDT<SparseMatType> solver(A);
    
    return 0;
}