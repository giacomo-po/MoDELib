/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DDtraitsIO_H_
#define model_DDtraitsIO_H_

#include <string>

namespace model
{
    
    
    struct DDtraitsIO
    {
        
        const std::string simulationFolder;
        const std::string inputFilesFolder;
        const std::string evlFolder;
        const std::string auxFolder;
        const std::string fFolder;
        const std::string ddFile;
        const std::string fFile;
        const std::string flabFile;
        const std::string polyFile;
        const std::string materialFile;
        const std::string microstructureFile;
        const std::string meshFile;

        DDtraitsIO(const std::string& folderName);
        
        
    };
    
}
#endif

