/* This file is part of finite element solution of BVP attached with model "the Mechanics of Defects Evolution Library".
 *
 * Copyright (C) 2011 by Mamdouh Mohamed <mamdouh.s.mohamed@gmail.com>, 
 * Copyright (C) 2011 by Giacomo Po <giacomopo@gmail.com>.
 * 
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef bvpfe_element_H_
#define bvpfe_element_H_

#include <vector>
//#include "model/Dislocations/Materials/Copper.h"
#include "model/BVP/Node.h"
//#include <model/Geometry/Simplex.h>
//#include <model/Geometry/Mesh/SimplexEnums.h>
#include <model/Mesh/Simplex.h>

#include "model/Quadrature/Quadrature.h"

namespace bvpfe{
	
	
	
	
	//	template<unsigned int n>
	template<short unsigned int dim>
	class Element {
		
	public: 
		
		//enum{Nnodes=model::Simplex<dim>::Nnodes};
//		enum{Nnodes=model::SimplexEnums<dim>::nVertices};
		enum{Nnodes=model::Simplex<dim,dim>::nVertices};
		
		//static model::Copper material;
		
		//		unsigned int Nnodes;             // number of nodes per element
		std::vector< Node<3>* > eleNodes;           // vector contains pointers to element's nodes
		
		public :
		
		Element ()
		{
			//			Nnodes = n; 
			eleNodes.reserve(Nnodes); 
		}
		
		// ------ function to insert the element's nodes ---------------
		
		void insertNode (Node<3>* nodePtr)
		{
			if (eleNodes.size() < Nnodes )
			{
				eleNodes.push_back(nodePtr);
			}
			else
			{
				std::cout<< "Error (Element.h: insertNode()): Excess number of Nodes" << std::endl;
				assert(0);
			}
		}
		//==============================================================
		// function to set the neighbor nodes for each node
		//==============================================================
		void setNodesNeighbors ()
		{
			for (int i = 0; i<Nnodes; i++ )
			{
				for (int j = i+1; j<Nnodes; j++ )
				{
					eleNodes[i]->setNeighbor(eleNodes[j]);	
				}
			}
		}
		
	};
	
	//template<short unsigned int n>
	//model::Copper Element<n>::material;
	
}  //  namespace bvpfe
#endif
