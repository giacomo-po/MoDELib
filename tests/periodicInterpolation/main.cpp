
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <set>
#include <math.h>       /* fmod */

#include <memory>
#include <vector>
#include <Eigen/Dense>
#include <TextFileParser.h>
#include <PeriodicLatticeInterpolant.h>



using namespace model;



int main(int argc, char * argv[])
{
 
    

    const Eigen::Matrix<double,2,2> A(TextFileParser("input.txt").readMatrix<double>("A",2,2,true));
    const Eigen::Matrix<double,2,1> N(TextFileParser("input.txt").readMatrix<double>("N",2,1,true));
    const Eigen::Matrix<double,2,1> D(Eigen::Matrix<double,2,1>::Ones());
    const Eigen::Matrix<double,Eigen::Dynamic,3> f(TextFileParser("input.txt").readMatrixCols<double>("f",3,true));
    const Eigen::Matrix<double,Eigen::Dynamic,5> df(TextFileParser("input.txt").readMatrixCols<double>("df",5,true));

    
    PeriodicLatticeInterpolant<2> pli(A,N.template cast<size_t>(),D.template cast<size_t>(),f,df);

    std::ofstream ofs("output.txt");
    for(int k=0;k<pli.waveVectors.rows();++k)
    {
        ofs<<std::setprecision(15)<<std::scientific<<pli.waveVectors.row(k)<<" "<<pli.sinCosCoeffs.row(k)<<"\n";
    }

    return 0;
}



//
//    const Eigen::Matrix<double,2,2> B(2.0*M_PI*A.inverse().transpose());
//
//    std::vector<Eigen::Matrix<double,2,1>> kPoints;
//
//    for(int i=0;i<N(0);++i)
//    {
//        for(int j=0;j<N(1);++j)
//        {
//            kPoints.push_back(B*(Eigen::Matrix<double,2,1>()<<double(i),double(j)).finished());
//        }
//    }
//
//    Eigen::MatrixXd M1(Eigen::MatrixXd::Zero(f.rows(),2*kPoints.size()));
//    Eigen::MatrixXd V1(Eigen::MatrixXd::Zero(f.rows(),1));
//    for(int i=0;i<f.rows();i++)
//    {
//        for(size_t j=0;j<kPoints.size();j++)
//        {
//            M1(i,2*j)=sin(kPoints[j].dot(f.block<1,2>(i,0)));
//            M1(i,2*j+1)=cos(kPoints[j].dot(f.block<1,2>(i,0)));
//        }
//
//        V1(i)=f(i,2);
//    }
//
//    std::cout<<"M1=\n"<<M1<<std::endl;
//    std::cout<<"V1=\n"<<V1<<std::endl;
//
//    Eigen::MatrixXd M2(Eigen::MatrixXd::Zero(df.rows(),2*kPoints.size()));
//    Eigen::MatrixXd V2(Eigen::MatrixXd::Zero(df.rows(),1));
//    for(int i=0;i<df.rows();i++)
//    {
//        for(size_t j=0;j<kPoints.size();j++)
//        {
//            M2(i,2*j)=kPoints[j].dot(df.block<1,2>(i,2))*cos(kPoints[j].dot(df.block<1,2>(i,0)));
//            M2(i,2*j+1)=-kPoints[j].dot(df.block<1,2>(i,2))*sin(kPoints[j].dot(df.block<1,2>(i,0)));
//        }
//
//        V2(i)=df(i,4);
//    }
//
//    std::cout<<"M2=\n"<<M2<<std::endl;
//    std::cout<<"V2=\n"<<V2<<std::endl;
//
//
//    Eigen::MatrixXd M(M1.rows()+M2.rows(),2*kPoints.size()-1);
//    M<<M1.block(0,1,M1.rows(),M1.cols()-1),M2.block(0,1,M2.rows(),M2.cols()-1);
//    Eigen::MatrixXd V(V1.rows()+V2.rows(),1);
//    V<<V1,V2;
//
//    std::cout<<M<<std::endl;
//    std::cout<<V<<std::endl;
//
//    Eigen::MatrixXd x((M.transpose()*M).llt().solve(M.transpose()*V));
//    Eigen::MatrixXd x1(x.rows()+1,1);
//    x1<<0.0,x;
//    std::cout<<x<<std::endl;
//
//        std::ofstream ofs("output.txt");
//        for(size_t k=0;k<kPoints.size();++k)
//        {
//            ofs<<std::setprecision(15)<<std::scientific<<x1(2*k)<<" "<<x1(2*k+1)<<" "<<kPoints[k].transpose()<<"\n";
//        }

