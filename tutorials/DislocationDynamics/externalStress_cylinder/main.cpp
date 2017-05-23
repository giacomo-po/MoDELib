/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2012 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#define _MODEL_NON_SINGULAR_DD_ 1 /* 1 = Cai's non-singular theory, 2 = Lazar's non-singular gradient theory */

#include <model/DislocationDynamics/DislocationNetwork.h>

using namespace model;

int main (int argc, char* argv[])
{
    // Create a DislocationNetwork object with:
    // dim=3, splineContinuity=1 (tangents are continuous), CatmullRom splines,
    // 16 quadrature points per segment Uniformly distributed
    DislocationNetwork<3,1,CatmullRom,UniformOpen> DN(argc,argv);
    
    // alternatively use GaussLegendre quadrature
    //DislocationNetwork<3,1,CatmullRom,16,GaussLegendre> DN(argc,argv);

    // Run time steps
    //DN.runSteps();
    DN.runDDD();
    

	
    return 0;
}
