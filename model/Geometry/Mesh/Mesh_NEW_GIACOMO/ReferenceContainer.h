#ifndef  mmdl_REFERENCECONTAINER_H_
#define mmdl_REFERENCECONTAINER_H_

#include <iostream>
#include <assert.h>

namespace mmdl {
	
	

	
	
	template <unsigned int N, typename T>
	struct ReferenceContainer {
		
		
		const T& ref;
		ReferenceContainer<N-1,T> LowerRefContainer;
		
		// Constructor
		ReferenceContainer(const ReferenceContainer<N-1,T>& LowerRefContainerIN , const T& refIN) : LowerRefContainer(LowerRefContainerIN), ref(refIN){
		std::cout<<"Creting ReferenceContainer<"<<N<<">"<<std::endl;	
		}
		
		// Copy Constructor
		ReferenceContainer(const ReferenceContainer<N,T>& other) : LowerRefContainer(other.LowerRefContainer), ref(other.ref){
		std::cout<<"Creting ReferenceContainer<"<<N<<">"<<std::endl;	
		}
		
		// operator &
		ReferenceContainer<N+1,T> operator&(const T& refIN){
		return  ReferenceContainer<N+1,T>(*this, refIN);
		}
		
		// operator ()
		const T& operator()(const unsigned int& k){
			assert(k<N);
			return (k==N-1)? ref : LowerRefContainer(k);
		}
		
		
	};
	
	
	
	template <typename T>
	struct ReferenceContainer<2,T> {
		
		enum {N=2};
		const T& ref0;
		const T& ref1;
		
		// Constructor
		ReferenceContainer(const T& refIN0, const T& refIN1) : ref0(refIN0), ref1(refIN1){
		std::cout<<"Creting ReferenceContainer<"<<2<<">"<<std::endl;	
		}
			
		// Copy Constructor
		ReferenceContainer(const ReferenceContainer<N,T>& other) : ref0(other.ref0), ref1(other.ref1){
		std::cout<<"Creting ReferenceContainer<"<<N<<">"<<std::endl;	
		}
		
		// operator &
		ReferenceContainer<N+1,T> operator&(const T& refIN){
		return  ReferenceContainer<N+1,T>(*this, refIN);
		}
			
		// operator ()
		const T& operator()(const unsigned int& k){
			assert(k<2);
			return (k==0)? ref0 : ref1;
		}
		
	};
	
	
	
	template <typename Derived>
	struct ReferenceContainee {
		
	ReferenceContainer<2,Derived > operator&(const Derived& other){
	return ReferenceContainer<2,Derived>(*static_cast<Derived*>(this),other);
	} 	
		
	};
	
	
	
	
}
#endif