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
#include "mmdl/BVP/SphericalIndenter.h"
#include <mmdl/Dislocations/Materials/Copper.h>

typedef std::pair< unsigned int, std::pair<unsigned int, double> > inputBCsType;

std::vector <inputBCsType> update_usr_BCs() const {

  double zSurf = 16.0e03;
  
  std::vector <inputBCsType> usrBCsContainer;
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
 return usrBCsContainer;   
    /*--First, remove previous boundary conditions.
     *You can comment this line if you will just add new boundary conditions without changing the existing ones (you need to be very careful) --*/

    //shared.domain.removeBoundaryConditions();

    /*====================================================================*/
    /*write your own code to set the new boundary conditions*/
    /*====================================================================*/

    /* You can access the information for node i by shared.domain.nodeContainer[i].
     * The node position is "shared.domain.nodeContainer[i].P" , which is Eigen::Matrix<double,3,1>
     * The traction on the node is "shared.domain.nodeContainer[i].traction" , which is Eigen::Matrix<double,3,1>
     * "shared.domain.nodeContainer[i].isBoundaryNode" is true only if the node is on the boundary
     * You can set the displacement boundary condition for any node i by calling "shared.domain.nodeContainer[i].setBC(ix, value)"
     * where "ix" is an index for the DOF you want to fix (0 = x, 1=y, 2=z), and "value" is the BC value you want to set.
     * We can also remove a specific DOF for any node i by calling "shared.domain.nodeContainer[i].removeBC(ix)"

     */

    /*Always check if the node is a boundary node before setting its boundary condition, or removing it by "if(shared.domain.nodeContainer[i].isBoundaryNode)"*/

}
//==================================================================================================
// function to check if some nodes are over-constraint, so they need to be freed and resolve the BVP
//==================================================================================================

bool overConstraintNodes () const {
  //enum {dim=3};

  double zSurf = 20000.0;

  #include "mmdl/BVP/Tetrahedron.h"
  #include <map>
  #include <vector>
  
  std::map<unsigned int, double> constraintNodes;

  bvpfe::Tetrahedron* pTet;
  Eigen::Matrix<double,12,12> kTet;
  Eigen::Matrix<double,12,1> F , U;
  
  unsigned int id ;
  unsigned int nRelaxed = 0;
  

    bool resolve = false;
    bool onSurface;
    
    double t0=clock();
    std::cout<<"Checking overcontraint nodes.... ";
        
    for (unsigned int it=0; it < shared.domain.tetContainer.size(); it++ ){
      pTet = &shared.domain.tetContainer[it];
            
      onSurface = false;
      
      for (unsigned int in=0; in<4; in++){
	if(pTet->eleNodes[in]->P(2)==zSurf && pTet->eleNodes[in]->isBC(2)) onSurface = true;
      }
      
      if(!onSurface) continue;
      
      kTet = pTet->getElementStiffness();

      for (unsigned int ii = 0; ii < 4; ii++){
	for (unsigned int jj = 0; jj < 3; jj++){
	  id = (ii*3) + jj;
	  U(id) = pTet->eleNodes[ii]->u(jj);
	}
      }
      
      F = kTet * U;
      
      for (unsigned int in = 0; in < 4; in++){
	if(!pTet->eleNodes[in]->isBC(2)) continue;
	id = (in*3) + 2; 
	if(constraintNodes.find(pTet->eleNodes[in]->sID) == constraintNodes.end() ) 
	  constraintNodes.insert(std::pair<unsigned int,double>(pTet->eleNodes[in]->sID,F(id)));
	else
	  constraintNodes[pTet->eleNodes[in]->sID] += F(id);
      }
      
    }
    
    for (std::map<unsigned int, double>::iterator it=constraintNodes.begin() ; it != constraintNodes.end(); it++ ) {
      if(it->second < 0.0) continue; 
      
      shared.domain.nodeContainer[it->first].removeBC(2);
      nRelaxed++;
      resolve = true;
      //std::cout<< "Node freed " << shared.domain.nodeContainer[it->first].sID<< "    " << it->second <<std::endl;
    }
      
      std::cout<<"DONE: "<<nRelaxed << " nodes were relaxed ..." <<"["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]"<<std::endl;
      /*if(!isFinal){
	onBoarder = false;
	for(unsigned int inn=0; inn< shared.domain.nodeContainer[it->first].neighbor.size(); inn++){
	  if ((!shared.domain.nodeContainer[it->first].neighbor[inn]->isBC(2))&&(shared.domain.nodeContainer[it->first].neighbor[inn]->isBoundaryNode))  onBoarder = true; 
	}
	if(!onBoarder)  continue;
      }*/
	
      //nodes2beFreed.push_back(it->first); 
      //std::cout<< it->second << std::endl;   
      
    return resolve;

}

//=======================================================================================================
// Function to calculate the indentation load
//======================================================================================================

double indentationLoad() const  {
  #include "mmdl/BVP/Triangle.h"
  //enum {dim=3};
  double zSurf = 20000.0;

  
    bvpfe::Triangle* pTri;
    unsigned int iTet;
    double load  = 0.0;
    double factor;
    //Eigen::Matrix<double,dim,dim> stress = Eigen::Matrix<double,dim,dim>::Zero();
    Eigen::Matrix<double,dim,1> forceVec;
    
    for (unsigned int it=0; it < shared.domain.triContainer.size(); it++ ){
      pTri = shared.domain.triContainer[it];

      if (pTri->eleNodes[0]->P(2)!=zSurf || pTri->eleNodes[1]->P(2)!=zSurf || pTri->eleNodes[2]->P(2)!=zSurf) continue;
      
      // Do the load calculations only if any of the triangle nodes is underneath the indenter
      if (pTri->eleNodes[0]->isBC(2) || pTri->eleNodes[1]->isBC(2) || pTri->eleNodes[2]->isBC(2)) {
	
	// if the triangle element is completely underneath the indenter, its full load will be considered
	if (pTri->eleNodes[0]->isBC(2) && pTri->eleNodes[1]->isBC(2) && pTri->eleNodes[2]->isBC(2)) factor = 1.0;
	
	//  otherwise, only half of it will be considered
	else factor = 0.5;
	
	iTet = pTri->neighTetIndx;
	forceVec = pTri->deformed_area() * ( shared.domain.tetContainer[iTet].getStress() * pTri->triNormalDeformed() );
	load +=  factor*forceVec(2);
	
      }
      
      /*
      // Check if the element is underneath the indenter
      
      if (pTri->eleNodes[0]->isBC(2) && pTri->eleNodes[1]->isBC(2) && pTri->eleNodes[2]->isBC(2)) {
	iTet = pTri->neighTetIndx;
	forceVec = pTri->deformed_area() * ( shared.domain.tetContainer[iTet].getStress() * pTri->triNormalDeformed() );
	load +=  forceVec(2);
      }
      */
    }

    return -1*load*48.0e-6*0.2556*0.2556;  // load in milli Newton units 
}

#endif

