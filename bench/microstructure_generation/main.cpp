/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <model/Mesh/SimplicialMesh.h>

using namespace model;



int main(int argc, char** argv)
{
    
    const int dim=3;
    typedef Eigen::Matrix<double,dim,1> VectorDim;
    
    
    // Take meshID as a user input
    int meshID(1);
    if (argc>1)
    {
        meshID=atoi(argv[1]);
    }
    
    // Create a 3d-SimplicialMesh object
    SimplicialMesh<dim> mesh;
    // Read the mesh files ./T/T_meshID.txt and ./N/N_meshID.txt
    mesh.readMesh(meshID);
    

    VectorDim P;
    P<<98.0,67.6,0.7;
    
    std::pair<bool,const Simplex<dim,dim>*> pair=mesh.search(P);
    
    if(pair.first)
    {
        std::cout<<"Point inside"<<std::endl;
        std::cout<<"Tet ID="<<pair.second->xID<<std::endl;
    }
    
    return 0;
}




