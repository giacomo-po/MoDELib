/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2012 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

/* Define the non-singluar method used for calculations
 * _MODEL_NON_SINGULAR_DD_ = 0 classical theory
 * _MODEL_NON_SINGULAR_DD_ = 1 Cai's regularization method
 * _MODEL_NON_SINGULAR_DD_ = 2 Lazar's regularization method
 */
#define _MODEL_NON_SINGULAR_DD_ 0
#include <model/DislocationDynamics/DislocationNetwork.h>

using namespace model;

int main (int argc, char* argv[])
{
    // Create the DislocationNetwork object
    typedef DislocationNetwork<3,1,CatmullRom,16,UniformOpen> DislocationNetworkType;
    DislocationNetworkType DN(argc,argv);
    
    
    double L=1000.0;    // max cell size
    int N=10;           // subdivisions of cell size
    if (argc>1)
    {
        L=atof(argv[1]);
        if (argc>2)
        {
            N=atoi(argv[2]);
        }
    }
    
    // Reach a stable number of nodes
    size_t nodeOrder=0;
    while(nodeOrder!=DN.nodeOrder())
    {
        nodeOrder=DN.nodeOrder();
        DN.remesh();
    }
    
    // Compute and output the PK force for different cell sizes
    // without multipole expansion
    DN.shared.use_StressMultipole=false;
    for(int n=1;n<=N;++n)
    {
        DislocationNetworkType::SpatialCellObserverType::setCellSize(n*L/N);
        DN.runSteps();
        DN.clearParticles(); // this also destroys all cells
    }

    // Compute and output the PK force for different cell sizes
    // with multipole expansion
    DN.shared.use_StressMultipole=true;
    for(int n=1;n<=N;++n)
    {
        DislocationNetworkType::SpatialCellObserverType::setCellSize(n*L/N);
        DN.runSteps();
        DN.clearParticles(); // this also destroys all cells
    }

    
    return 0;
}
