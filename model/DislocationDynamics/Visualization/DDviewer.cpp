/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#include <model/Mesh/SimplicialMesh.h> // defines mode::cout
#include <model/DislocationDynamics/Visualization/DDglut.h>


/************************************/
int main(int argc, char** argv)
{
//    model::SimplicialMesh<3> mesh(1);

	return model::DDglut(argc,argv);
    
}
