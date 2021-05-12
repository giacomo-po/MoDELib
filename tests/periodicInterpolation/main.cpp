
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


void output(const PeriodicLatticeInterpolant<2>& pli,std::vector<Eigen::Matrix<double,2,1>> pts)
{
    std::ofstream outputFile("output.txt");
    for(int k=0;k<pli.waveVectors.rows();++k)
    {
        outputFile<<std::setprecision(15)<<std::scientific<<pli.waveVectors.row(k)<<" "<<pli.sinCosCoeffs.row(k)<<"\n";
    }
    
    std::ofstream dataFile("data.txt");
    for(const auto& pt : pts)
    {
        dataFile<<pt(0)<<" "<<pt(1)<<" "<<pli(pt)<<"\n";
    }
}


int main(int argc, char * argv[])
{
    
    const int rotSymm(argc > 1 ? atoi(argv[1]) : 1);
    const bool fullSurface(argc > 2 ? atoi(argv[2]) : true);

    const Eigen::Matrix<double,2,2> A(TextFileParser("input.txt").readMatrix<double>("A",2,2,true));
    const Eigen::Matrix<double,Eigen::Dynamic,3> f(TextFileParser("input.txt").readMatrixCols<double>("f",3,true));
    const Eigen::Matrix<double,Eigen::Dynamic,2> m(TextFileParser("input.txt").readMatrixCols<double>("m",2,true));
    const Eigen::Matrix<double,Eigen::Dynamic,2> waveVectorBasisCoordinates(TextFileParser("input.txt").readMatrixCols<double>("waveVecInt",2,true));
    
    std::vector<Eigen::Matrix<double,2,1>> mirSymm;
    for(int i=0;i<m.rows();i++)
    {
        mirSymm.push_back(m.row(i).transpose());
    }
    
    
    PeriodicLatticeInterpolant<2> pli(A,waveVectorBasisCoordinates,f,rotSymm,mirSymm);
    
    
    std::vector<Eigen::Matrix<double,2,1>> pts;
    const Eigen::VectorXd x(TextFileParser("points.txt").readMatrixCols<double>("x",1,false));
    const Eigen::VectorXd y(TextFileParser("points.txt").readMatrixCols<double>("y",1,false));

    if(fullSurface)
    {

        for(int i=0;i<x.size();++i)
        {
            for(int j=0;j<y.size();++j)
            {
                pts.push_back((Eigen::Matrix<double,2,1>()<<x(i),y(j)).finished());
            }
        }
    }
    else
    {
        assert(x.size()==y.size());
        for(int i=0;i<x.size();++i)
        {
            pts.push_back((Eigen::Matrix<double,2,1>()<<x(i),y(i)).finished());
        }
    }
    
    output(pli,pts);


    return 0;
}

