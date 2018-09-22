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
//#include <model/DislocationDynamics/Materials/PeriodicElement.h>
#include <model/DislocationDynamics/MobilityLaws/DislocationMobility.h>
#include <model/DislocationDynamics/Materials/Material.h>


using namespace model;

int main (int argc, char * const argv[])
{
    

    
    
//    PeriodicElement<74,Isotropic> W;
  
    Material<3,Isotropic> material("./W.txt");
    
     Eigen::Matrix<double,3,1> s ;  // Burgers vector
     Eigen::Matrix<double,3,1> n ;  // plane normal
    
    if(material.crystalStructure=="BCC")
    {
         s = (Eigen::Matrix<double,3,1>()<< 1.0,1.0,1.0).finished().normalized();
         n = (Eigen::Matrix<double,3,1>()<<-1.0,0.0,1.0).finished().normalized();
    }
    else if(material.crystalStructure=="FCC")
    {
         s = (Eigen::Matrix<double,3,1>()<<-1.0, 0.0,1.0).finished().normalized();
         n = (Eigen::Matrix<double,3,1>()<< 1.0,1.0,1.0).finished().normalized();
    }
    else
    {
        assert(0 && "TEST NOT SUPPORTED FOR MATERIAL TYPE");
    }
    
    const Eigen::Matrix<double,3,1> m = (Eigen::Matrix<double,3,1>()<<0.0, 0.0,1.0).finished().normalized(); // load direction

    const Eigen::Matrix<double,3,1> xi = s; // screw

    
    TextFileParser parser(material.materialFile);
    const double B0e=parser.readScalar<double>("B0e_SI",true);
    const double B1e=parser.readScalar<double>("B1e_SI",true);
    const double B0s=parser.readScalar<double>("B0s_SI",true);
    const double B1s=parser.readScalar<double>("B1s_SI",true);
    const double Bk=parser.readScalar<double>("Bk_SI",true);
    const double dH=parser.readScalar<double>("dH0_eV",true);
    const double p=parser.readScalar<double>("p",true);
    const double q=parser.readScalar<double>("q",true);
    const double Tm=parser.readScalar<double>("Tm",true);
    const double Tf=parser.readScalar<double>("Tf",true);
    const double tauC=parser.readScalar<double>("tauC_SI",true);
    const double a0=parser.readScalar<double>("a0",true);
    const double a1=parser.readScalar<double>("a1",true);
    const double a2=parser.readScalar<double>("a2",true);
    const double a3=parser.readScalar<double>("a3",true);
    
    DislocationMobility<BCClattice<3>> mobility(material.b_SI,
                                       material.mu_SI,
                                       material.cs_SI,
                                       B0e,B1e,
                                       B0s,B1s,
                                       Bk,
                                       dH,
                                       p,q,
                                       Tm*Tf,
                                       tauC,
                                       a0,a1,a2,a3);
    
    std::deque<double> TAU;
    TAU.push_back(material.mu_SI*pow(10.0,-5));
    for(int e=-5;e<=-2;++e)
    {
        for(int k=2;k<=10;++k)
        {
            TAU.push_back(material.mu_SI*k*pow(10.0,e));
        }
    }
    
    const int nT=101;
    std::deque<double> T;
    for(int n=0;n<nT;++n)
    {
        T.push_back(Tm*n/nT);
    }
    
    
    
    std::ofstream v_file("velocity.txt");

    for (const auto& t : T)
    {
        for(const auto& tau : TAU)
        {
            const double sf=m.dot(s)*m.dot(n);
            const Eigen::Matrix<double,3,3> S=tau/sf*m*m.transpose();

            v_file<<std::setprecision(15)<<std::scientific<<tau/material.mu_SI<<" "<<t/Tm<<" "<<mobility.velocity(S/material.mu_SI,s,xi,n,t,0.0,0.0,false)<<std::endl;
        }
    }
    
    return 0;
}


