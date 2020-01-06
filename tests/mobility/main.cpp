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
#include <DislocationMobility.h>
#include <Material.h>


using namespace model;

int main (int argc, char * const argv[])
{
    
    const std::string materialFile= argc > 1 ? std::string(argv[1]) : "../../tutorials/DislocationDynamics/MaterialsLibrary/W.txt";
    
    //    PeriodicElement<74,Isotropic> W;
    
    Material<3,Isotropic> material(materialFile);
    TextFileParser parser(material.materialFile); // CREATE UNI_PTR OF THIS AND PUT IN IF STATEMENT
    
    Eigen::Matrix<double,3,1> s ;  // Burgers vector
    Eigen::Matrix<double,3,1> n ;  // plane normal
    
    std::unique_ptr<DislocationMobilityBase> mobility;
    
    if(material.crystalStructure=="BCC")
    {
        s = (Eigen::Matrix<double,3,1>()<< 1.0,1.0,1.0).finished().normalized();
        n = (Eigen::Matrix<double,3,1>()<<-1.0,0.0,1.0).finished().normalized();
        mobility=std::make_unique<DislocationMobility<BCClattice<3>>>(material);
    }
    else if(material.crystalStructure=="FCC")
    {
        s = (Eigen::Matrix<double,3,1>()<<-1.0, 0.0,1.0).finished().normalized();
        n = (Eigen::Matrix<double,3,1>()<< 1.0,1.0,1.0).finished().normalized();
        mobility=std::make_unique<DislocationMobility<FCClattice<3>>>(material);
    }
    else
    {
        assert(0 && "TEST NOT SUPPORTED FOR MATERIAL TYPE");
    }
    
    const Eigen::Matrix<double,3,1> m = (Eigen::Matrix<double,3,1>()<<0.0, 0.0,1.0).finished().normalized(); // load direction
    
    const Eigen::Matrix<double,3,1> xi = s; // screw
    
    
    // Construct array containing stress values
    std::deque<double> TAU;
    TAU.push_back(material.mu_SI*pow(10.0,-5));
    for(int e=-5;e<=-2;++e)
    {
        for(int k=2;k<=10;++k)
        {
            TAU.push_back(material.mu_SI*k*pow(10.0,e));
        }
    }
    
    // Construct array containing temperature values
    const int nT=101;
    std::deque<double> T;
    for(int n=0;n<nT;++n)
    {
        T.push_back(material.Tm*n/nT);
    }
    
    
    
    std::ofstream v_file("velocity.txt");
    
    for (const auto& t : T)
    {
        for(const auto& tau : TAU)
        {
            const double sf=m.dot(s)*m.dot(n);
            const Eigen::Matrix<double,3,3> S=tau/sf*m*m.transpose();
            
            v_file<<std::setprecision(15)<<std::scientific<<tau/material.mu_SI<<" "<<t/material.Tm<<" "<<mobility->velocity(S/material.mu_SI,s,xi,n,t,0.0,0.0,false)<<std::endl;
        }
    }
    
    return 0;
}


