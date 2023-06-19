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
                 const std::string& outFolderSuffix_in="");
        
        std::string getBinFilename(const size_t& runID);
        std::string getTxtFilename(const size_t& runID);
        std::string getTxtSegmentFilename(const size_t& runID);
        bool isBinGood(const size_t& frameID);
        bool isTxtGood(const size_t& frameID);
        
        template<typename T>
        static void binWrite(std::ofstream& of,const T& o)
        {
            of.write((char *) &o, (sizeof o));
        }
                
        template<typename T>
        static void binRead(std::ifstream& file,
                            T*& memblock,
                            const size_t& arraySize)
        {
            memblock = new T [arraySize];
            file.read (reinterpret_cast<char*>(memblock), arraySize*sizeof (T));
        }
        
    };
    
}
#endif

