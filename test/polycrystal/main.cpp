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
    int meshID(2);
    if (argc>1)
    {
        meshID=atoi(argv[1]);
    }
    
//    Eigen::Matrix<double,3,3> A;
//    A<< 0.0,1.0,1.0,
//    /**/1.0,0.0,1.0,
//    /**/1.0,1.0,0.0;
    


    // Create a 3d-SimplicialMesh object
    SimplicialMesh<3> mesh;
    // Read the mesh files ./T/T_meshID.txt and ./N/N_meshID.txt
    mesh.readMesh(meshID);
    
    Polycrystal<3> poly(mesh);

    poly.grain(1).selectMaterial(29);
    poly.grain(2).selectMaterial(13);

    
    Eigen::Matrix<double,3,3> C2G1;
    C2G1.row(1)<<0.0,1.0,0.0;
    C2G1.row(2)<<(Eigen::Matrix<double,1,3>()<<3.0,0.0,-1.0).finished().normalized();
    C2G1.row(0)<<C2G1.row(1).cross(C2G1.row(2));
    std::cout<<C2G1<<std::endl;
    poly.grain(1).rotate(C2G1);

    Eigen::Matrix<double,3,3> C2G2;
    C2G2.row(1)<<0.0,1.0,0.0;
    C2G2.row(2)<<(Eigen::Matrix<double,1,3>()<<-3.0,0.0,-1.0).finished().normalized();
    C2G2.row(0)<<C2G2.row(1).cross(C2G2.row(2));
        std::cout<<C2G2<<std::endl;
    poly.grain(2).rotate(C2G2);

    
    poly.grainBoundary(1,2).createLatticePlanes();
    
    Eigen::Matrix<double,3,1> p;
    p<<0,11,10;
    
    std::cout<<poly.grainBoundary(1,2).latticePlane(1).snapToLattice(p)<<std::endl;
    std::cout<<poly.grainBoundary(1,2).latticePlane(2).snapToLattice(p)<<std::endl;

    
//    
//    Eigen::Matrix<double,3,3> C2G2;
//    C2G2<<0,		1,		0,
//    0.9486832981,		0,		0.316227766,
//    0.316227766,	0,		-0.9486832981;
//    poly.grain(2).setLatticeBasis(C2G2*A);
//    
//    std::cout<<poly.grain(1).region.regionID<<std::endl;
//    
//    
//    std::cout<<poly.grainBoundary(1,2).regionBoundary.size()<<std::endl;
//    
////    poly.grainBoundary(3,5)
//    Eigen::Matrix<double,3,1> p0;
//    p0<<500.,500.0,499.0;
//    
//    Eigen::Matrix<double,3,1> p1;
//    p1<<500.,500.0,501.0;
//    
////    LatticeVector<3> L(p,)
//    
//    LatticeVector<3> L0(poly.latticeVectorFromPosition(p0));
//    LatticeVector<3> L1(poly.latticeVectorFromPosition(p1));
//    LatticeVector<3> L2=L0;  // assignment operator
//    LatticeVector<3> L3=L0+L1;  // assignment operator
//    
//    ReciprocalLatticeVector<3> R0(poly.reciprocalLatticeVectorFromPosition(p0));  // assignment operator
//    ReciprocalLatticeVector<3> R1(poly.reciprocalLatticeVectorFromPosition(p1));  // assignment operator
//    
//    LatticeDirection<3> LD0=R0.cross(R1);  // assignment operator
//
//    
//    ReciprocalLatticeDirection<3> RD0=L0.cross(L1);  // assignment operator
//
//    
//    std::cout<<LD0.cartesian()<<std::endl;
//    std::cout<<LD0<<std::endl;
////    
//    std::cout<<RD0.cartesian()<<std::endl;
//    std::cout<<RD0<<std::endl;
//    
//    
//    LatticeLine line0(L1,L1-L0);
//    line0.snapToLattice(p1);
//    std::cout<<"contains? "<<line0.contains(L0)<<std::endl;
//    
//    LatticeLine line1(L0,L0);
//
//    LineLineIntersection lli(line0,line1);
//    std::cout<<lli<<std::endl;
    
    return 0;
}




