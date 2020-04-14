/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#include <string>
#include <experimental/filesystem>
#include <DDconfigIO.h>

using namespace model;

namespace fs = std::experimental::filesystem;


int main(int argc, char** argv)
{
    const std::string path = "./evl";

    const std::string subS1="evl_";
    const std::string subS2=".bin";
    
    const bool writeSegments = argc > 1 ? std::atoi(argv[1]) : false;


    for (const auto & p : fs::directory_iterator(path))
    {
        
        const size_t found1=p.path().string().find(subS1);
        const size_t found2=p.path().string().find(subS2);
        
        if (   found1 != std::string::npos
            && found2 != std::string::npos)
        {
            const std::string stringID=p.path().string().substr(found1+subS1.length(),found2-(found1+subS1.length()));
            const size_t runID=std::stol(stringID);
            DDconfigIO<3> converter;
            converter.bin2txt(runID,writeSegments);
        }
    }

    return 0;
    
}


