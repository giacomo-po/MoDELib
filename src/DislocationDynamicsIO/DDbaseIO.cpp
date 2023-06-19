/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DDaux_cpp_
#define model_DDaux_cpp_


#include <iostream>     // std::cout
#include <fstream>      // std::ifstream
#include <string>      // std::ifstream
#include <DDbaseIO.h>

namespace model
{
    
        
        DDbaseIO::DDbaseIO(const std::string& outFolderName_in,
                 const std::string& fileNamePrefix_in,
                 const std::string& outFolderSuffix_in) :
        /* init */ outFolderName(outFolderName_in)
        /* init */,fileNamePrefix(fileNamePrefix_in)
        /* init */,outFolderSuffix(outFolderSuffix_in)
        {}
        

        std::string DDbaseIO::getBinFilename(const size_t& runID)
        {
            return outFolderName+outFolderSuffix+"/"+fileNamePrefix+"_"+std::to_string(runID)+".bin";
        }
        
        std::string DDbaseIO::getTxtFilename(const size_t& runID)
        {
            return outFolderName+outFolderSuffix+"/"+fileNamePrefix+"_"+std::to_string(runID)+".txt";
        }
        
        std::string DDbaseIO::getTxtSegmentFilename(const size_t& runID)
        {
            return outFolderName+outFolderSuffix+"/"+"segments"+std::to_string(runID)+".txt";
        }
        
        bool DDbaseIO::isBinGood(const size_t& frameID)
        {
            return std::ifstream(getBinFilename(frameID).c_str(), std::ios::in|std::ios::binary).good();

        }

        bool DDbaseIO::isTxtGood(const size_t& frameID)
        {
            return std::ifstream(getTxtFilename(frameID).c_str(), std::ios::in).good();
            
        }

    
}
#endif

