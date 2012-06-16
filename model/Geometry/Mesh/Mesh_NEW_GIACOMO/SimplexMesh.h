#ifndef  mmdl_SIMPLEXMESH_H_
#define  mmdl_SIMPLEXMESH_H_

#include <boost/ptr_container/ptr_map.hpp>
#include <Eigen/Dense>

//#include <mmdl/Mesh/Simplex.h>
//#include "mmdl/Mesh/Mesh_NEW_GIACOMO/Simplex.h"

//#include <mmdl/Utilities/AddressBook.h>


#include "mmdl/Geometry/Mesh/Mesh_NEW_GIACOMO/SimplexObserver.h"
#include "mmdl/Geometry/Mesh/Mesh_NEW_GIACOMO/SimplexID.h"
#include "mmdl/Geometry/Mesh/Mesh_NEW_GIACOMO/Simplex.h"


//#include <mmdl/Mesh/MeshVertex.h>
//#include <mmdl/Mesh/SimplexEnums.h>
//#include <mmdl/Mesh/ReferenceContainer.h>


namespace mmdl {

	
//	template<short unsigned int dim>
//	class SimplexObserver : public AddressBook<Simplex<dim>,0>,
//	/*                   */ public SimplexObserver<dim-1>{};
//	
//	template<>
//	class SimplexObserver<0> : public AddressBook<Simplex<0>,0>{};
	

	template<short unsigned int dim>
	class SimplexMesh : /*private boost::ptr_map<size_t,MeshVertex<dim> >,*/
	/*        */ private boost::ptr_map<size_t,Simplex<dim> >{


		enum {Nnodes=SimplexEnums<dim>::Nnodes};
		typedef Eigen::Matrix<double,dim,1>	VectorDimD;
		typedef Eigen::Matrix<int,Nnodes,1> VectorNnodesI;
		
		
		
//		template<unsigned int N>
//		ReferenceContainer<N+1,MeshVertex<dim> > referenceNodes(const VectorNnodesI& ){
//		}
		
		
	public:
//		void read(){
//		
//			// insert nodes
////			for (int n=0;n<Nnodes;++n){
////				VectorDimD P;
////				std::auto_ptr<MeshVertex> pV(new MeshVertex(P));
////				size_t vertexID = pV->sID;
////				this->insert(vertexID, pV);
////			}
//
//			// insert simplices
////			for (int n=0;n<Nsimplices;++n){
//			
//			switch (dim) {
//				case 1:
//					<#statements#>
//					break;
//				case 2:
//					<#statements#>
//					break;
//				case 3:
//					<#statements#>
//					break;
//				default:
//					assert(0 && "SIMPLEX MESH ONLY AVAILABLE IN 1d, 2d and 3d.")
//					break;
//			}
//			
////				ReferenceContainer<Nnodes,MeshVertex<dim> > RC(V0 & V1 & V2 & V3);
////			}
//			
//		}
		
		
		
	};
	

	
	
	
}
#endif