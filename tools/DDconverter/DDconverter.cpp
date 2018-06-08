/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

//#include <iostream>
//#include <fstream>
#include <string>
//#include <sstream>
//#include <iomanip>

#include <model/DislocationDynamics/IO/EVLio.h>
//#include <model/Network/Readers/VertexReader.h>
//#include <model/Network/Readers/EdgeReader.h>
#include <experimental/filesystem>

using namespace model;

namespace fs = std::experimental::filesystem;


int main(int argc, char** argv)
{
    const std::string path = "./evl";

    const std::string subS1="evl_";
    const std::string subS2=".bin";

    
    for (const auto & p : fs::directory_iterator(path))
    {
        
        const size_t found1=p.path().string().find(subS1);
        const size_t found2=p.path().string().find(subS2);
        
        if (   found1 != std::string::npos
            && found2 != std::string::npos)
        {
            const std::string stringID=p.path().string().substr(found1+subS1.length(),found2-(found1+subS1.length()));
            const size_t runID=std::stol(stringID);
            
            EVLio<3> evlio;
            evlio.bin2txt(runID);
        }
    }

    
    
//    DDreader ddr;
//    ddr.list();
//
//    for (const auto& file : ddr)
//    {
//        const size_t fileID=file.first;
//        
//        VertexReader<'V',9,double> vReader;
//        if(vReader.isGood(fileID,0)) // V_fileID.bin is good to read
//        {
//            vReader.read(fileID,0);
//            std::ostringstream filestream;
//            filestream <<"V/V" << "_" << fileID << ".txt";
//            std::ofstream outFile(filestream.str().c_str(), std::ios::out);
//            for (const auto& v : vReader)
//            {
//                outFile<<v.first<<" "<<std::setprecision(15)<<std::scientific<<v.second<<"\n";
//            }
//        }
//        
//        EdgeReader  <'E',11,double> eReader;
//        if(eReader.isGood(fileID,0)) // E_fileID.bin is good to read
//        {
//            eReader.read(fileID,0);
//            std::ostringstream filestream;
//            filestream <<"E/E" << "_" << fileID << ".txt";
//            std::ofstream outFile(filestream.str().c_str(), std::ios::out);
//            for (const auto& e : eReader)
//            {
//                outFile<<e.first.first<<" "<<e.first.second<<" "<<std::setprecision(15)<<std::scientific<<e.second<<"\n";
//            }
//        }
//        
//        
//
//    }
    
    return 0;
    
}


