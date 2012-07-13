/* This file is part of finite element solution of BVP attached with mmdl "the Mechanics of Material Defects Library".
 *
 * Copyright (C) 2011 by Mamdouh Mohamed <mamdouh.s.mohamed@gmail.com>,
 * Copyright (C) 2011 by Giacomo Po <giacomopo@gmail.com>.
 *
 * mmdl is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef bvpfe_UpdateBoundaryConditions_H_
#define bvpfe_UpdateBoundaryConditions_H_

#include <Eigen/Dense>
#include <mmdl/Dislocations/DislocationSharedObjects.h>
#include <mmdl/Dislocations/Materials/Copper.h>

typedef std::pair< unsigned int, std::pair<unsigned int, double> > inputBCsType;

std::vector <inputBCsType> update_usr_BCs() const {

    /* write your own code to set the new boundary conditions
     * You can access the information for node i by shared.domain.nodeContainer[i].
     * The node position is "shared.domain.nodeContainer[i].P" , which is Eigen::Matrix<double,3,1>
     * The traction on the node is "shared.domain.nodeContainer[i].traction" , which is Eigen::Matrix<double,3,1>
     * "shared.domain.nodeContainer[i].isBoundaryNode" is true only if the node is on the boundary
     * You can set the displacement boundary condition for any node i by calling "shared.domain.nodeContainer[i].setBC(ix, value)"
     * where "ix" is an index for the DOF you want to fix (0 = x, 1=y, 2=z), and "value" is the BC value you want to set.
     * We can also remove a specific DOF for any node i by calling "shared.domain.nodeContainer[i].removeBC(ix)"
     
     */
    
    /*Always check if the node is a boundary node before setting its boundary condition, or removing it by "if(shared.domain.nodeContainer[i].isBoundaryNode)"*/

    
    
  //double zSurf = 16.0e03;
  
  std::vector <inputBCsType> usrBCsContainer;
  /*
  inputBCsType u_BC;

    static double prescribed_displ(0.0);
    
  for (unsigned int in = 0; in < shared.domain.nodeContainer.size(); in++ ) {
      if (!shared.domain.nodeContainer[in].isBoundaryNode) continue;
      if ( shared.domain.nodeContainer[in].P(2)!=zSurf) continue;
	  u_BC.first = in;
	  u_BC.second.first = 2;
	  u_BC.second.second = prescribed_displ;
	  usrBCsContainer.push_back(u_BC);
  }
    prescribed_displ += 0.1;
    */
 return usrBCsContainer;   
    /*--First, remove previous boundary conditions.
     *You can comment this line if you will just add new boundary conditions without changing the existing ones (you need to be very careful) --*/

    //shared.domain.removeBoundaryConditions();


}


#endif

