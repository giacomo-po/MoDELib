/* This file is part of finite element solution of BVP attached with model "the Mechanics of Material Defects Library".
 *
 * Copyright (C) 2011 by Mamdouh Mohamed <mamdouh.s.mohamed@gmail.com>, 
 * Copyright (C) 2011 by Giacomo Po <giacomopo@gmail.com>.
 * 
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef bvpfe_node_H_
#define bvpfe_node_H_

#include <algorithm>
#include "model/Utilities/StaticID.h"
#include "model/BVP/Point.h"
#include "model/BVP/Dof.h"
//#include "model/BVP/Triangle.h"

namespace bvpfe{


	//class Node;
	template<short unsigned int dim = 3>
	class Node :	public Point<dim>, 
			public model::StaticID<Node<dim> >,
			public Dof<dim> {		

#include <model/BVP/commonTypeDefs.h>				
						
		public :
						
			bool isBoundaryNode;			

			//VectorDim uInf;
			std::vector<unsigned int> triIDs;   // IDs for surface triangle that share this node			
			VectorDim Pc;	                 // current node position         	
						
			std::vector<Node*> neighbor;    // set of neighbor nodes

			//Dof<dim> u;                          // dof of the node (displacement)
			
			VectorDim traction;                // value for the traction that is used to interpolate for the traction at any point of the surface triangle
			
			Node (VectorDim& P_in) : Point<dim>::Point(P_in){
			  
				neighbor.push_back(this);     // make the node neighbor for it self (just to make creating sparse system easier)
				//std::cout<< this->sID<< "  " << neighbor[0]->sID << std::endl;
				Pc=this->P;
				
				traction = VectorDim::Zero();
				isBoundaryNode=false;
			}

			//================================================
			// function to set nodes nieghbors
			//==============================================

			inline void setNeighbor (Node* pN)
			{
				bool exist;
				exist = 0;

				for (unsigned int i = 0; i<neighbor.size(); i++)
				{
					if (neighbor[i]->sID == pN->sID)  
					{exist = 1;   break; }
				}

				if (! exist)
				{
					neighbor.push_back(pN);	
					pN->neighbor.push_back(this);
				}
			}
						
			//================================
			// function to update the nodes position
			//================================
						
			void displaceNode()
			{
				this->Pc=this->P + (this->u + this->uInf*0);
			}

			//====================================
			// function used in sorting neighbors vector
			//====================================

			void sortNeighbors()
			{
				Node* temp;

				for (unsigned int i = 0; i<neighbor.size(); i++ )
				{
					for (unsigned int j = i+1 ; j<neighbor.size(); j++ )
					{
						if(neighbor[i]->sID > neighbor[j]->sID) 
						{
							temp = neighbor[i];
							neighbor[i] = neighbor[j];
							neighbor[j] = temp;
						}
					} 	
				} 
			
				
			}
	};
		
	

}  //  namespace bvpfe
#endif
