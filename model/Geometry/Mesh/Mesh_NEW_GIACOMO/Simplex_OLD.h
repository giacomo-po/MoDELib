/* This file is part of mmdl, the Mechanics of Material Defects Library.
 *
 * Copyright (C) 2011 by Giacomo Po <giacomopo@gmail.com>.
 *
 * mmdl is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef  mmdl_SIMPLEX_H_
#define  mmdl_SIMPLEX_H_

#include <boost/shared_ptr.hpp>

#include <mmdl/Utilities/StaticID.h>


//#include <mmdl/Mesh/ReferenceContainer.h>
//#include <mmdl/Mesh/ReferenceContainer.h>
//#include <mmdl/Mesh/MeshVertex.h>
#include "mmdl/Mesh/Mesh_NEW_GIACOMO/SimplexEnums.h"
//#include <mmdl/Mesh/SimplexEnums.h>

#include <Eigen/Dense>
#include <math.h>  // for std::abs
#include <utility> // for std::pair

namespace mmdl {
	
	
	
	
	
//	//	template<short unsigned int order>
//	//	class Simplex : public StaticID<Simplex<order> >,
//	//	/*           */ public ReferenceContainer<SimplexEnums<order>::Nnodes,MeshVertex<order> > {
//	
//	//		public:
//	
//	//		enum {Nnodes=SimplexEnums<order>::Nnodes};
//	//		enum {Nedges=SimplexEnums<order>::Nedges};
//	
//	//		Simplex(const ReferenceContainer<SimplexEnums<order>::Nnodes,MeshVertex<order> >& NodeRefsIN) : ReferenceContainer<SimplexEnums<order>::Nnodes,MeshVertex<order> >(NodeRefsIN){}
//	
//	
//	
//	// returns the k-the simplex of order lowerorder. Eg. in a 3d simplex get<2>(0) returns the fist face
//	// Eg. in a 3d simplex get<0>(0) returns the fist node
//	//	template <short unsigned int lowerorder>
//	//	const Simplex<lowerorder>& get(const short unsigned int k) const {
//	//	return ;
//	//	}
//	
//	//	};
//	
//	
//	template<typename T>
//	class Min{
//		 T min;
//
//	public:
//		template<typename...Tpacket>
//		Min(Tpacket... a) {}
//	
//	};
//	
//	
//	
//	
//		
//	
//	
////	template<short unsigned int N, typename T>
////	class SimplexID : private SimplexID<N-1,T>{
////		
////		const T stored_int;
////		
////		template <T...Tpacket>
//////		T keepMin(const T& i1, const T& i2, const Tpacket&... someT) const {
////			T keepMin(const T& i1, const T& i2, Tpacket) const {
////
////			assert(i1!=i2);
////		return (i1<i2)? keepMin(i1,someT...) : keepMin(i2,someT...); 
////		}
////
////		T keepMin(const T& i1, const T& i2) const {
////			assert(i1!=i2);
////			return (i1<i2)? i1 : i2;
////		}
////		
////	public:
////		template<T...Tpacket>
//////		SimplexID(const T& i_in, const Tpacket&...someT) : SimplexID<N-1,T>(someT...), stored_int(keepMin(i_in,someT...))
////		SimplexID(const T& i_in, Tpacket) : SimplexID<N-1,T>(Tpacket...), stored_int(keepMin(i_in,Tpacket...))
////		/*: i(keepMax(i_in,intPacket...)), SimplexID<N-1>(?????)*/ {} //here p
////		
////		
////		 T operator()(const int& k) const {
////			assert(k>=0 && k<N);
////			return (k==0)? stored_int : SimplexID<N-1,T>::operator()(k-1);
////		}
////		
////	};
//	
//	
//	
//	
////	sizeof...(Tpacket)
//
//	
////	template<typename T>
////	class SimplexID<1,T>{
////		T stored_int;
////		
////	public:
////		SimplexID(const T&i_in) : stored_int(i_in){}
////		
////		
////		const T& operator()(const int& k) const {
////			assert(k==0);
////			return stored_int;
////		}
////		
////	};
//	
//	
//	template<short unsigned int,short unsigned int>
//	//	class Simplex;
//	struct SimplexFaceExtractor;
//	
//	
//	
//	
//	// Each Simplex should create the simplices of order order-1 that it is composed of.
//	// Each Simplex should be 
//	
//	
//	template <unsigned int N, typename T>
//	class FixesSizeArray : FixesSizeArray<N-1,T> {
//		T elem;
//		
//	public:
//		T& operator()(const int& k) {
//			assert(k>=0 && k<N);
//			return (k==N-1)? elem : FixesSizeArray<N-1,T>(k);
//		}
//		
//		const T& operator()(const int& k) const {
//			assert(k>=0 && k<N);
//			return (k==N-1)? elem : FixesSizeArray<N-1,T>(k);
//		}
//	};
//	
//	template <typename T>
//	class FixesSizeArray<1,T> {
//		T elem;
//		
//		
//	public:
//		
//		T& operator()(const int& k) {
//			assert(k==0);
//			return elem;
//		}
//		
//		const T& operator()(const int& k) const {
//			assert(k==0);
//			return elem;
//		}
//		
//	};
	
	
	/******************************************************/
	/* Simplex<order>                                     */
	/******************************************************/
	template<short unsigned int order>
	struct SortSimplexID {
	
		typedef Simplex<dim,order> SimplexType;
		typedef SimplexType::SimplexIDType SimplexIDType;
		
		template<short unsigned int dim>
		static SimplexIDType sortIDs(const SimplexType& S){
		
		//return 
		}
		
	};
	
	/******************************************************/
	/* Simplex<order>                                       */
	/******************************************************/
	template<short unsigned int dim, short unsigned int order>
	class Simplex {
		
		static_assert( dim>0, "SIMPLEX DIM MUST BE > 0." ) ;
		static_assert( order<=dim, "SIMPLEX ORDER MUST BE <= DIM." ) ;
		
		typedef Simplex<order-1> Face;
		typedef Simplex<0> Vertex;
		
//		FixesSizeArray<order+1,boost::shared_ptr<Face> > faceContainer;
		// The pointer to the SubNetwork
		//boost::shared_ptr<SubNetworkType> psn;
		
		const Simplex<0>& last;
		const Face& lower;
		
	public:
		
		typedef Eigen::Matrix<unsigned int,1,order> SimplexIDType;
		
		typedef Eigen::Matrix<double,order,1> VectororderD;
		
		
		
		//enum {order=1};
		//enum {Nnodes=2};
		//enum {Nedges=1};
		enum {Nvertex=SimplexEnums<order>::Nnodes};
		
		/******************************************************/
		// Constructor
		Simplex (const Simplex<order-1>& lowerIN, const Simplex<0>& lastIN) : last(lastIN), lower(lowerIN){
			if(1>0){
				//faceContainer(0).reset;
				//this->psn.reset(new SubNetworkType(this->p_derived()));
			}
			else{
				//			psn=psnOther;		// redirect psn to the new Subnetwork
			}
			
			
			//	std::cout<<"Creating Simplex<"<<order<<">"<<std::endl;
		}
		
		/******************************************************/
		// Copy contructor
		Simplex(const Simplex<order>& other) : last(other.last), lower(other.lower){
			//	std::cout<<"Copying Simplex<"<<order<<">"<<std::endl;	
		}
		
		/******************************************************/
		// operator*
		Simplex<order+1> operator*(const Simplex<0>& other){
			return Simplex<order+1>(*this, other);
			//return boost::shared_ptr<Simplex<order+1> >
			//this->psn.reset(new SubNetworkType(this->p_derived()));
		}
		
		/******************************************************/
		// vertex
		const Simplex<0>& vertex(const int& k) const {
			/*! A const reference to the k-th vertex in the simplex. 
			 *  Vertices are of type Simplex<0>.
			 */
			assert(k>=0 && k<Nvertex);
			return (k==Nvertex-1)? last : lower.vertex(k);
		}
		
//		/******************************************************/
//		// face
//		Simplex<order-1> face(const int& k) const {
//			/*! A const reference to the k-th face in the simplex. 
//			 *  Faces are of type Simplex<order-1>. Faces are numbered such that face(k) does not include vertex(k).
//			 */
//			assert(k>=0 && k<Nvertex);
//			//		return SimplexFaceExtractor<order,order-1>(*this).face(k+1);
//			
//			return SimplexFaceExtractor<order,order-1>::face(*this,k+1);
//		}
		
		
		/******************************************************/
		// faceNormal
		VectororderD faceNormal(const int& k) const {
			
			// Algorithm is: 
			// make a Gram-Schmidt orthonormalization of the edges starting with the ones 
			// belonging to face k. Last element in the normalization is the normal.
			
		}
		
		//@
		//template <short unsigned int lowerorder>
		//	Simplex<order-1> face(const int& k) const {
		//	assert(k<=Nvertex);
		//	return Simplex<order-1>(this);
		//	}
		
	};
	
	/******************************************************/
	/* Simplex: template specialization Simplex<0>        */
	/******************************************************/	
	template<>
	class Simplex<0> : public StaticID<Simplex<0> > {
		
		
	public:
		//enum {order=1};
		enum {Nnodes=2};
		enum {Nedges=1};
		
		/* Constructor: this sould take a position ****************************/
		Simplex(){
			std::cout<<"Creating new Simplex<0> with ID "<<this->sID<<std::endl;
		}
		
		/* Destructor *********************************************************/
		~Simplex(){
			std::cout<<"Destroying Simplex<0> with ID "<<this->sID<<std::endl;
		}
		
		// operator *
		Simplex<1> operator*(const Simplex<0>& other){
//			if ()
			
			return Simplex<1>(*this,other);	
		}
		
		
		// operator ()
		const Simplex<0>& vertex(const int& k) const {
			assert(k==0);
			return *this;
		}
	};
	
//	template<short unsigned int simplexorder, short unsigned int currentorder=simplexorder-1>
//	struct SimplexFaceExtractor{
//		
//		enum {Nvertex=SimplexEnums<simplexorder>::Nnodes};
//		
//		static Simplex<currentorder> face(const Simplex<simplexorder> simplex, const int& k){
//			//		std::cout<<"k is"<<k<<std::endl;
//			//		std::cout<<"k%N is"<<k%Nvertex<<std::endl;
//			//				std::cout<<"simplex.vertex(k%Nvertex) is "<<simplex.vertex(k%Nvertex).sID<<std::endl;
//			
//			
//			//	std::cout<<"(k-1+simplexorder)%Nvertex is"<<(k+currentorder)%Nvertex<<std::endl;
//			
//			//	std::cout<<"simplex.vertex((k-1+simplexorder)%Nvertex) is "<<simplex.vertex((k+currentorder)%Nvertex).sID<<std::endl;
//			return Simplex<currentorder>( SimplexFaceExtractor<simplexorder,currentorder-1>::face(simplex,k%Nvertex),simplex.vertex((k+currentorder)%Nvertex) );		
//		}
//	};	
//	
//	template<short unsigned int simplexorder>
//	struct SimplexFaceExtractor<simplexorder,0>{
//		
//		enum {Nvertex=SimplexEnums<simplexorder>::Nnodes};
//		
//		static  const Simplex<0>& face(const Simplex<simplexorder> simplex, const int& k){
//			
//			return simplex.vertex(k%Nvertex);	
//		}
//	};
	
	
	
	
}
#endif
