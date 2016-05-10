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
#include <model/DislocationDynamics/Materials/PeriodicElement.h>
#include <model/DislocationDynamics/MobilityLaws/DislocationMobility.h>


using namespace model;

int main (int argc, char * const argv[])
{
    
    const Eigen::Matrix<double,3,1> m = (Eigen::Matrix<double,3,1>()<<-0.1706, 0.3711,0.9128).finished().normalized();
    const Eigen::Matrix<double,3,1> s = (Eigen::Matrix<double,3,1>()<<-1.0, 0.0,1.0).finished().normalized();
    const Eigen::Matrix<double,3,1> n = (Eigen::Matrix<double,3,1>()<< 1.0,1.0,1.0).finished().normalized();
    const Eigen::Matrix<double,3,1> xi = s; // screw
    
    
    PeriodicElement<74,Isotropic> W;
    
    std::deque<double> TAU;
    TAU.push_back(W.mu*pow(10.0,-5));
    for(int e=-5;e<=-2;++e)
    {
        for(int k=2;k<=10;++k)
        {
            TAU.push_back(W.mu*k*pow(10.0,e));
        }
    }
    
    const int nT=101;
    std::deque<double> T;
    for(int n=0;n<nT;++n)
    {
        T.push_back(W.Tm*n/nT);
    }
    
    
    
    std::ofstream v_file("velocity.txt");

    for (const auto& t : T)
    {
        for(const auto& tau : TAU)
        {
            const double sf=m.dot(s)*m.dot(n);
            const Eigen::Matrix<double,3,3> S=tau/sf*m*m.transpose();

            v_file<<std::setprecision(15)<<std::scientific<<tau/W.mu<<" "<<t/W.Tm<<" "<<W.dm.velocity(S/W.mu,s,xi,n,t)<<std::endl;
        }
    }
    
    return 0;
}


