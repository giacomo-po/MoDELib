/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#include <MicrostructureGenerator.h>

using namespace model;

int main(int argc, char** argv)
{
    
    const std::string folderName(argc>1? std::string(argv[1]) : "./");
    
    DislocationDynamicsBase<3> ddBase( folderName);
    MicrostructureGenerator mg( ddBase);
    return 0;
}
