#ifndef  mmdl_MESHVERTEX_H_
#define mmdl_MESHVERTEX_H_

#include <iostream>
#include <mmdl/Mesh/Mesh_NEW_GIACOMO/ReferenceContainer.h>
#include <mmdl/Utilities/StaticID.h>

//#include <mmdl/Utilities/CRTP.h>

namespace mmdl {
	
	

	
	
	template <short unsigned int dim>
	class MeshVertex : public StaticID<MeshVertex<dim> >, 
				      public ReferenceContainee<MeshVertex<dim> > {
		
	public:	
		
	MeshVertex(){
	std::cout<<"Creating MeshVertex<"<<dim<<"> "<<this->sID<<std::endl;	
		
	}
	
//	ReferenceContainer<2,MeshVertex<dim> > operator&(const MeshVertex<dim>& other){
//	return ReferenceContainer<2,MeshVertex<dim> >(*this,other);
//	} 		
		
	};
	
	
	
	
}
#endif