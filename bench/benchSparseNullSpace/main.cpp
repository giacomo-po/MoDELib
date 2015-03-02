//g++ main.cpp -o test -O3 -std=c++11 -fopenmp -I/usr/local/include/ -I.


#include <iostream>
#include <chrono>
#include <omp.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseQR>
#include <Eigen/OrderingMethods>

#include <SparseNullSpace.h>

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

    SparseNullSpace<SparseMatType> ns(A);
    
//    SparseMatType AT=A.transpose();
//    Eigen::SparseQR<SparseMatType,Eigen::COLAMDOrdering<int> > qr(AT);
//    SparseMatType Q;
//    Q=qr.matrixQ();
//    SparseMatType R;
//    R=qr.matrixR();
//    
//    std::cout<<Eigen::MatrixXd(Q)<<std::endl<<std::endl;
//
//    std::cout<<Eigen::MatrixXd(R)<<std::endl<<std::endl;
//
//
////    std::cout<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<std::endl;
//
//    SparseMatType Z;
//    Z=Q;
////
////    std::cout<<Eigen::MatrixXd(Z)<<std::endl<<std::endl;
//
    SparseMatType AZ;
    AZ=A*ns.matrixZ();
//    
    std::cout<<Eigen::MatrixXd(AZ)<<std::endl;

    
    return 0;
}