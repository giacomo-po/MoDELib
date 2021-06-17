/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DDaux_H_
#define model_DDaux_H_

#include <iostream>     // std::cout
#include <fstream>      // std::ifstream
#include <string>      // std::ifstream


namespace model
{
    
    struct DDbaseIO
    {
        
        const std::string outFolderName;
        const std::string fileNamePrefix;
        const std::string outFolderSuffix;
        
        DDbaseIO(const std::string& outFolderName_in,
                 const std::string& fileNamePrefix_in,
                 const std::string& outFolderSuffix_in="") :
        /* init */ outFolderName(outFolderName_in)
        /* init */,fileNamePrefix(fileNamePrefix_in)
        /* init */,outFolderSuffix(outFolderSuffix_in)
        {}
        
        /**********************************************************************/
        template<typename T>
        static void binWrite(std::ofstream& of,const T& o)
        {
            of.write((char *) &o, (sizeof o));
        }
        
        /**********************************************************************/
        template<typename T>
        static void binRead(std::ifstream& file,
                            T*& memblock,
                            const size_t& arraySize)
        {
            memblock = new T [arraySize];
            file.read (reinterpret_cast<char*>(memblock), arraySize*sizeof (T));
        }
        
        /**********************************************************************/
        std::string getBinFilename(const size_t& runID)
        {
            return outFolderName+outFolderSuffix+"/"+fileNamePrefix+"_"+std::to_string(runID)+".bin";
        }
        
        /**********************************************************************/
        std::string getTxtFilename(const size_t& runID)
        {
            return outFolderName+outFolderSuffix+"/"+fileNamePrefix+"_"+std::to_string(runID)+".txt";
        }
        
        /**********************************************************************/
        std::string getTxtSegmentFilename(const size_t& runID)
        {
            return outFolderName+outFolderSuffix+"/"+"segments"+std::to_string(runID)+".txt";
        }
        
        /**********************************************************************/
        bool isBinGood(const size_t& frameID)
        {
            return std::ifstream(getBinFilename(frameID).c_str(), std::ios::in|std::ios::binary).good();

        }

        /**********************************************************************/
        bool isTxtGood(const size_t& frameID)
        {
            return std::ifstream(getTxtFilename(frameID).c_str(), std::ios::in).good();
            
        }

    };
    
}
#endif

