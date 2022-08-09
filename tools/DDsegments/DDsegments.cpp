/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#include <VTKsegments.h>

using namespace model;

int main(int argc, char** argv)
{
    
    const std::string folderName(argc>1? std::string(argv[1]) : "./");
    const std::string vtkFile(argc>2? std::string(argv[2]) : "dxa_ovito_0");

    VTKsegments vs(folderName);
    vs.readVTK(vtkFile);
    
    return 0;
}
