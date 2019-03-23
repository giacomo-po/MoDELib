/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2012 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#define _MODEL_NON_SINGULAR_DD_ 0 /* 1 = Cai's non-singular theory, 2 = Lazar's non-singular gradient theory */

#include <model/DislocationDynamics/DislocationNetwork.h>

using namespace model;

int main (int argc, char* argv[])
{
    // Create a DislocationNetwork object
        
    DislocationNetwork<3,1,CatmullRom,UniformOpen> DN(argc,argv);
    // Run the simulation
    DN.runSteps();
    
    return 0;
}

