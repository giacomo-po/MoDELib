/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#include <model/Mesh/GmshReader.h>

int main (int argc, char * const argv[])
{
    model::GmshReader reader("block_structured1.msh");
    
    std::cout<<reader.nodes().size()<<" nodes"<<std::endl;
    std::cout<<reader.elements().size()<<" elements"<<std::endl;
    return 0;
}


