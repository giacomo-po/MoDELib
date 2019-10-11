/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2012 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#define _MODEL_NON_SINGULAR_DD_ 1 // 0 classical theory, 1 Cai's regularization method, 2 Lazar's regularization method
#define userLoadController "./LoadController.h"

#include <DefectiveCrystal.h>

using namespace model;

int main (int argc, char* argv[])
{
    // Create the DislocationNetwork object
    //    DislocationNetwork<3,0,Hermite> DN(argc,argv);
    DefectiveCrystal<3,0,Hermite> DC(argc,argv);
    
    // Run time steps
    DC.runGlideSteps();
    
    return 0;
}

