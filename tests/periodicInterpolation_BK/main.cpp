
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


void output(const PeriodicLatticeInterpolant<2>& pli)
{
    std::ofstream ofs("output.txt");
    for(int k=0;k<pli.waveVectors.rows();++k)
    {
        ofs<<std::setprecision(15)<<std::scientific<<pli.waveVectors.row(k)<<" "<<pli.sinCosCoeffs.row(k)<<"\n";
    }
}


int main(int argc, char * argv[])
{
    
    const bool readWaveVectors(argc > 1 ? atoi(argv[1]) : 0);
    
    const Eigen::Matrix<double,2,2> A(TextFileParser("input.txt").readMatrix<double>("A",2,2,true));
    const Eigen::Matrix<double,Eigen::Dynamic,3> f(TextFileParser("input.txt").readMatrixCols<double>("f",3,true));
    const Eigen::Matrix<double,Eigen::Dynamic,5> df(TextFileParser("input.txt").readMatrixCols<double>("df",5,true));

    if(readWaveVectors)
    {
        const Eigen::Matrix<double,Eigen::Dynamic,2> waveVectorBasisCoordinates(TextFileParser("input.txt").readMatrixCols<double>("waveVecInt",2,true));
        PeriodicLatticeInterpolant<2> pli(A,waveVectorBasisCoordinates,f,df);
        output(pli);
    }
    else
    {
        const Eigen::Matrix<double,2,1> N(TextFileParser("input.txt").readMatrix<double>("N",2,1,true));
        const Eigen::Matrix<double,2,1> D(TextFileParser("input.txt").readMatrix<double>("D",2,1,true));
        PeriodicLatticeInterpolant<2> pli(A,N.template cast<size_t>(),D.template cast<size_t>(),f,df);
        output(pli);
    }
    
    return 0;
}

