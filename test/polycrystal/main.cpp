/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

// The following line is needed for shape functions or order 5 or more
//#define EIGEN_STACK_ALLOCATION_LIMIT 1000000


#include <model/Utilities/SequentialOutputFile.h>
#include <model/Mesh/SimplicialMesh.h>
#include <model/DislocationDynamics/Polycrystals/Polycrystal.h>
#include <model/LatticeMath/LatticeVector.h>


using namespace model;



int main(int argc, char** argv)
{
    
    // Take meshID as a user input
    int meshID(2);
    if (argc>1)
    {
        meshID=atoi(argv[1]);
    }


    // Create a 3d-SimplicialMesh object
    SimplicialMesh<3> mesh;
    // Read the mesh files ./T/T_meshID.txt and ./N/N_meshID.txt
    mesh.readMesh(meshID);
    
    Polycrystal<3> poly(mesh);
    poly.read("./DDinput.txt");

    poly.grain(1).selectMaterial(29);
    poly.grain(2).selectMaterial(13);

    
//    Eigen::Matrix<double,3,3> C2G1;
//    C2G1.row(1)<<0.0,1.0,0.0;
//    C2G1.row(2)<<(Eigen::Matrix<double,1,3>()<<3.0,0.0,-1.0).finished().normalized();
//    C2G1.row(0)<<C2G1.row(1).cross(C2G1.row(2));
//    std::cout<<std::setprecision(15)<<std::scientific<<C2G1<<std::endl;
//    poly.grain(1).rotate(C2G1);
//
//    Eigen::Matrix<double,3,3> C2G2;
//    C2G2.row(1)<<0.0,1.0,0.0;
//    C2G2.row(2)<<(Eigen::Matrix<double,1,3>()<<-3.0,0.0,-1.0).finished().normalized();
//    C2G2.row(0)<<C2G2.row(1).cross(C2G2.row(2));
//    std::cout<<std::setprecision(15)<<std::scientific<<C2G2<<std::endl;
//    poly.grain(2).rotate(C2G2);

    
    Eigen::Matrix<double,3,1> p;
    p<<25.78,11,110.6;
    
    std::cout<<poly.grainBoundary(1,2).latticePlane(1).snapToLattice(p)<<std::endl;
    std::cout<<poly.grainBoundary(1,2).latticePlane(2).snapToLattice(p)<<std::endl;

    

    
    return 0;
}




