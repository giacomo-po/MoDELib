#include <iostream>
#include <iomanip>      // std::setprecision
#include <fstream>
#include <chrono>
#include <deque>
#include <string>     // std::string, std::to_string
#include <Eigen/Dense>
#include <model/Math/SmithDecomposition.h>
#include <model/LatticeMath/CSL.h>
#include <model/LatticeMath/DSCL.h>
#include <model/DislocationDynamics/Materials/FCCcrystal.h>
#include <model/DislocationDynamics/Materials/SimpleCubiccrystal.h>
#include <model/DislocationDynamics/Materials/BCCcrystal.h>

using namespace model;


void fillMap(std::deque<Eigen::Matrix<double,3,3>>& gbMap)
{
    // Rotation matrices from Grimmer, Acta Cryst. (1974). A30, 197, table 1

    // [100]
    gbMap.push_back((Eigen::Matrix<double,3,3>()<<2,-1,2,2,2,-1,-1,2,2).finished()/3);
    gbMap.push_back((Eigen::Matrix<double,3,3>()<<5,0,0,0,4,-3,0,3,4).finished()/5);
    gbMap.push_back((Eigen::Matrix<double,3,3>()<<13,0,0,0,12,-5,0,5,12).finished()/13);
    
}

template<typename LatticeType>
void out_csl_dscl(const Eigen::Matrix<double,3,3>& R)
{
    const int dim=3;
    Lattice<dim> A;
    Lattice<dim> B;
    A.setLatticeBasis(LatticeType::template getLatticeBasis<dim>());
    B.setLatticeBasis(R*LatticeType::template getLatticeBasis<dim>());
    CSL<dim> csl(A,B);
    DSCL<dim> dscl(A,B);
    std::ofstream file (std::string(LatticeType::name) + "_sigma_" + std::to_string(csl.sigma()) + ".txt", std::ofstream::out);
    file<<std::setprecision(15)<<std::scientific<<A.covBasis()<<std::endl;
    file<<std::setprecision(15)<<std::scientific<<B.covBasis()<<std::endl;
    file<<std::setprecision(15)<<std::scientific<<csl.covBasis()<<std::endl;
    file<<std::setprecision(15)<<std::scientific<<dscl.covBasis()<<std::endl;
}

int main()
{
    
    
    
    std::deque<Eigen::Matrix<double,3,3>> gbMap;
    fillMap(gbMap);
    
    for (const auto& R : gbMap)
    {
        
//        // SimpleCubic
        out_csl_dscl<SimpleCubic>(R);
        
        // FCC
        out_csl_dscl<FCC>(R);
        
        // BCC
        out_csl_dscl<BCC>(R);
        
    }
    
    
    
    return 0;
}
