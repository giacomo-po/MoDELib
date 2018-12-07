/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <random>
#include <vector>
#include <set>



int main (int argc, char * const argv[])
{
 
    const double m = argc > 1 ? atof(argv[1]) : 0.0;
    const double s = argc > 2 ? atof(argv[2]) : 1.0;
    
    std::cout<<"mean="<<m<<std::endl;
    std::cout<<"std="<<s<<std::endl;
    
    const int nrolls=1000000;  // number of experiments
    
    std::default_random_engine generator;
    std::lognormal_distribution<double> distribution(m,s);
    
    std::set<double> values;
    for (int i=0; i<nrolls; ++i)
    {
        values.insert(distribution(generator));
    }
    std::cout<<values.size()<<std::endl;
    
    std::ofstream file("probability.txt");

    for(const double& val : values)
    {
        file<<std::setprecision(15)<<std::scientific<<val<<"\n";
    }
    

    
    return 0;
}


