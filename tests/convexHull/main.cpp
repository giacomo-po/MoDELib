/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <fstream>

#include <IDreader.h>
#include <ConvexHull.h>


using namespace model;

int main (int argc, char * const argv[])
{
    
    typedef model::HullPoint<2,double> HullPointType;
    typedef model::ConvexHull<2,double> ConvexHullType;
    
    std::ifstream inFile("inPoints.txt");
    if(inFile.is_open())
    {

        HullPointType p;
        ConvexHullType hull;
        while (inFile >> p[0] >> p[1])
        {
//            hull.push_back(p); // &x[0] is unused
            hull.insert(p); // &x[0] is unused

        }
        
        const auto outPoints=hull.getPoints();
        std::ofstream outFile("outPoints.txt");
        for(const auto& point : outPoints)
        {
            outFile<<point[0]<<" "<<point[1]<<"\n";
        }
    }
    else
    {
        std::cout<<"Cannot open file inPoint.txt"<<std::endl;
    }

    
    
    

    
    return 0;
}


