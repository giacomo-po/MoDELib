/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DDtraitsIO_cpp_
#define model_DDtraitsIO_cpp_

#include <DDtraitsIO.h>
#include <TextFileParser.h>

namespace model
{
    
    

        DDtraitsIO::DDtraitsIO(const std::string& folderName) :
        /* init */ simulationFolder(folderName)
        /* init */,inputFilesFolder(simulationFolder+"/inputFiles")
        /* init */,evlFolder(simulationFolder+"/evl")
        /* init */,auxFolder(simulationFolder+"/evl")
        /* init */,fFolder(simulationFolder+"/F")
        /* init */,ddFile(inputFilesFolder+"/DD.txt")
        /* init */,fFile(fFolder+"/F_0.txt")
        /* init */,flabFile(fFolder+"/F_labels.txt")
        /* init */,polyFile(inputFilesFolder+"/polycrystal.txt")
        /* init */,materialFile(inputFilesFolder+"/"+TextFileParser(polyFile).readString("materialFile",false))
        /* init */,microstructureFile(inputFilesFolder+"/initialMicrostructure.txt")
        /* init */,meshFile(inputFilesFolder+"/"+TextFileParser(polyFile).readString("meshFile",false))
        {
            std::cout<<"simulationFolder="<<simulationFolder<<std::endl;
            std::cout<<"inputFilesFolder="<<inputFilesFolder<<std::endl;
            std::cout<<"evlFolder="<<evlFolder<<std::endl;
            std::cout<<"auxFolder="<<auxFolder<<std::endl;
            std::cout<<"fFolder="<<fFolder<<std::endl;
            std::cout<<"ddFile="<<ddFile<<std::endl;
            std::cout<<"fFile="<<fFile<<std::endl;
            std::cout<<"flabFile="<<flabFile<<std::endl;
            std::cout<<"polyFile="<<polyFile<<std::endl;
            std::cout<<"materialFile="<<materialFile<<std::endl;
            std::cout<<"microstructureFile="<<microstructureFile<<std::endl;
            std::cout<<"meshFile="<<meshFile<<std::endl;
        }
        
        

    
}
#endif

