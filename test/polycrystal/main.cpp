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
#include <model/LatticeMath/LatticeVector.h>


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
    Eigen::Matrix<double,3,1> p0;
    p0<<500.,500.0,499.0;
    
    Eigen::Matrix<double,3,1> p1;
    p1<<500.,500.0,501.0;
    
//    LatticeVector<3> L(p,)
    
    LatticeVector<3> L0(poly.latticeVectorFromPosition(p0));
    LatticeVector<3> L1(poly.latticeVectorFromPosition(p1));
    LatticeVector<3> L2=L0;  // assignment operator
    LatticeVector<3> L3=L0+L1;  // assignment operator
    
//    std::cout<<poly.latticeVectorFromPosition(p)<<std::endl;
    
    return 0;
}




