/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>

#include <model/DislocationDynamics/Visualization/DDreader.h>
#include <model/Network/Readers/VertexReader.h>
#include <model/Network/Readers/EdgeReader.h>

using namespace model;

int main(int argc, char** argv)
{
    
    
    
    DDreader ddr;
    ddr.list();

    for (const auto& file : ddr)
    {
        const size_t fileID=file.first;
        
        VertexReader<'V',10,double> vReader;
        if(vReader.isGood(fileID,0)) // V_fileID.bin is good to read
        {
            vReader.read(fileID,0);
            std::ostringstream filestream;
            filestream <<"V/V" << "_" << fileID << ".txt";
            std::ofstream outFile(filestream.str().c_str(), std::ios::out);
            for (const auto& v : vReader)
            {
                outFile<<v.first<<" "<<std::setprecision(15)<<std::scientific<<v.second<<"\n";
            }
        }
        
        EdgeReader  <'E',11,double> eReader;
        if(eReader.isGood(fileID,0)) // E_fileID.bin is good to read
        {
            eReader.read(fileID,0);
            std::ostringstream filestream;
            filestream <<"E/E" << "_" << fileID << ".txt";
            std::ofstream outFile(filestream.str().c_str(), std::ios::out);
            for (const auto& e : eReader)
            {
                outFile<<e.first.first<<" "<<e.first.second<<" "<<std::setprecision(15)<<std::scientific<<e.second<<"\n";
            }
        }
        
        

    }
    
    return 0;
    
}


