/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

// The following line is needed for shape functions or order 5 or more
//#define EIGEN_STACK_ALLOCATION_LIMIT 1000000


#include <iostream>
#include <model/Utilities/SequentialOutputFile.h>
#include <model/Mesh/SimplicialMesh.h>
#include <model/DislocationDynamics/Polycrystals/Polycrystal.h>


using namespace model;



int main(int argc, char** argv)
{
    
    // Take meshID as a user input
    int meshID(0);
    if (argc>1)
    {
        meshID=atoi(argv[1]);
    }
    
    // Create a 3d-SimplicialMesh object
    SimplicialMesh<3> mesh;
    // Read the mesh files ./T/T_meshID.txt and ./N/N_meshID.txt
    mesh.readMesh(meshID);
    
    Polycrystal<3> poly(mesh);
    
    std::cout<<poly.grain(1).region.regionID<<std::endl;
    
    
    std::cout<<poly.grainBoundary(1,2).regionBoundary.size()<<std::endl;
    
//    poly.grainBoundary(3,5)
    
    return 0;
}




