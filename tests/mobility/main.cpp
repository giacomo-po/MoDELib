/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <fstream>
#include <deque>
#include <memory>
//#include <model/DislocationDynamics/Materials/PeriodicElement.h>
#include <DislocationMobilityFCC.h>
#include <DislocationMobilityBCC.h>
#include <DislocationMobilityHEXprismatic.h>
#include <DislocationMobilityHEXbasal.h>
#include <DislocatedMaterialBase.h>


using namespace model;

int main (int argc, char * const argv[])
{
    
    const std::string materialFile= argc > 1 ? std::string(argv[1]) : "../../tutorials/DislocationDynamics/MaterialsLibrary/Zr.txt";
    
    //    PeriodicElement<74,Isotropic> W;
    
    DislocatedMaterialBase material(materialFile);
//    TextFileParser parser(material.materialFile); // CREATE UNI_PTR OF THIS AND PUT IN IF STATEMENT
    
    Eigen::Matrix<double,3,1> s ;  // Burgers vector
    Eigen::Matrix<double,3,1> n ;  // plane normal
    Eigen::Matrix<double,3,1> m;
    
    std::unique_ptr<DislocationMobilityBase> mobility;
    
    if(material.crystalStructure=="BCC")
    {
        s = (Eigen::Matrix<double,3,1>()<< 1.0,1.0,1.0).finished().normalized();
        n = (Eigen::Matrix<double,3,1>()<<-1.0,0.0,1.0).finished().normalized();
        m = (Eigen::Matrix<double,3,1>()<<0.0,0.0,1.0).finished().normalized();
        mobility=std::make_unique<DislocationMobilityBCC>(material);
    }
    else if(material.crystalStructure=="FCC")
    {
        s = (Eigen::Matrix<double,3,1>()<<-1.0,0.0,1.0).finished().normalized();
        n = (Eigen::Matrix<double,3,1>()<< 1.0,1.0,1.0).finished().normalized();
        m = (Eigen::Matrix<double,3,1>()<<0.0,0.0,1.0).finished().normalized();
        mobility=std::make_unique<DislocationMobilityFCC>(material);
    }
    else if(material.crystalStructure=="HEX")
    {
        s = (Eigen::Matrix<double,3,1>()<< 1.0,0.0,0.0).finished().normalized();
        m = (Eigen::Matrix<double,3,1>()<<1.0,2.0,3.0).finished().normalized();
        bool isBasal(true);
        if(isBasal)
        {
            n = (Eigen::Matrix<double,3,1>()<< 0.0,0.0,1.0).finished().normalized();
            mobility=std::make_unique<DislocationMobilityHEXbasal>(material);
        }
        else
        {
            n = (Eigen::Matrix<double,3,1>()<< 0.0,1.0,0.0).finished().normalized();
            mobility=std::make_unique<DislocationMobilityHEXprismatic>(material);
        }
    }
    else
    {
        assert(0 && "TEST NOT SUPPORTED FOR MATERIAL TYPE");
    }
    
    
    // Construct array containing stress values
    std::deque<double> TAU;
    TAU.push_back(material.mu_SI*pow(10.0,-5));
    const int nS=100;
    for(int e=-5;e<=-3;++e)
    {
        for(int k=1;k<nS;++k)
        {
            TAU.push_back(material.mu_SI*pow(10.0,e+double(k)/nS));
        }
    }
    
    // Construct array containing temperature values
    const int nT=101;
    std::deque<double> T;
    for(int n=0;n<nT;++n)
    {
        T.push_back(material.Tm*n/nT);
    }
    
    
    const Eigen::Matrix<double,3,1> xiS = s; // screw
    const Eigen::Matrix<double,3,1> xiE = n.cross(s).normalized(); // screw
    
    std::ofstream vS_file("velocityS.txt");
    std::ofstream vE_file("velocityE.txt");

    for (const auto& t : T)
    {
        for(const auto& tau : TAU)
        {
            const double sf=m.dot(s)*m.dot(n);
            const Eigen::Matrix<double,3,3> S=tau/sf*m*m.transpose();
            
            vS_file<<std::setprecision(15)<<std::scientific<<tau/material.mu_SI<<" "<<t/material.Tm<<" "<<mobility->velocity(S/material.mu_SI,s,xiS,n,t,0.0,0.0,false)<<std::endl;
            vE_file<<std::setprecision(15)<<std::scientific<<tau/material.mu_SI<<" "<<t/material.Tm<<" "<<mobility->velocity(S/material.mu_SI,s,xiE,n,t,0.0,0.0,false)<<std::endl;
        }
    }
    
    return 0;
}


