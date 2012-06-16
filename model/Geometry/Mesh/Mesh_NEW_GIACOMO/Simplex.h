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
#include <mmdl/Utilities/CompareVectorsByComponent.h>

//#include <mmdl/Utilities/TypeTraits.h>


//#include <mmdl/Mesh/ReferenceContainer.h>
//#include <mmdl/Mesh/ReferenceContainer.h>
//#include <mmdl/Mesh/MeshVertex.h>
#include "mmdl/Mesh/Mesh_NEW_GIACOMO/SimplexEnums.h"
#include <mmdl/Mesh/Mesh_NEW_GIACOMO/ObservingMap.h>

//#include <mmdl/Mesh/SimplexEnums.h>

#include <Eigen/Dense>
#include <math.h>  // for std::abs
#include <utility> // for std::pair

namespace mmdl {
	
	
	// A SIMPLEX OF ORDER 0 IS A VERTEX
	// A SIMPLEX OF ORDER N IS OBTAINED FROM A SIMPLEX OF ORDER N-1 and A SIMPLEX OF ORDER N. THE CONNECTION CREATES NEW EDGES(1-SIMPLEX) AND FACES(N-1 SIMPEX)
	// HOWEVER THESE NEW EDGES AND FACES MAY BE SHARED AMONG SIMLICES IN A MESH, THEREFORE THEY SHOULD BE HANDLED BY SHARED POINTERS
	// ONLY THE UPPER LEVEL SIMPLICES SHOULD BE STORED. EVERYTHING ELSE SHOULD BE HANDLED BY SHARED PTRS SO THAT DESTRUCTION OD A TETRA ALSO DESTROYS EVERYTHIN ELSE
	
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
	template<short unsigned int dim, short unsigned int order>
	class Simplex : public ObservedByMap<Eigen::Matrix<unsigned int,order,1>,Simplex<dim,order>,CompareVectorsByComponent<unsigned int,order> > {
		
		static_assert( dim>0, "SIMPLEX DIM MUST BE > 0." ) ;
		static_assert( order<=dim, "SIMPLEX ORDER MUST BE <= DIM." ) ;
		
		typedef Simplex<dim,order-1> Face;
		typedef Simplex<dim,0> Vertex;
		
		//		FixesSizeArray<order+1,boost::shared_ptr<Face> > faceContainer;
		// The pointer to the SubNetwork
		//boost::shared_ptr<SubNetworkType> psn;
		
		//		const Vertex& last;
		//		const Face& lower;
		
	public:
		
		//		typedef Eigen::Matrix<unsigned int,1,order> SimplexIDType;
		//		
		//		typedef Eigen::Matrix<double,order,1> VectororderD;
		//		
		//		
		//		
		//		//enum {order=1};
		//		//enum {Nnodes=2};
		//		//enum {Nedges=1};
		//		enum {Nvertex=SimplexEnums<order>::Nnodes};
		//		
		//		/******************************************************/
		//		// Constructor
		//		Simplex (const Simplex<order-1>& lowerIN, const Simplex<0>& lastIN) : last(lastIN), lower(lowerIN){
		//			if(1>0){
		//				//faceContainer(0).reset;
		//				//this->psn.reset(new SubNetworkType(this->p_derived()));
		//			}
		//			else{
		//				//			psn=psnOther;		// redirect psn to the new Subnetwork
		//			}
		//			
		//			
		//			//	std::cout<<"Creating Simplex<"<<order<<">"<<std::endl;
		//		}
		//		
		//		/******************************************************/
		//		// Copy contructor
		//		Simplex(const Simplex<order>& other) : last(other.last), lower(other.lower){
		//			//	std::cout<<"Copying Simplex<"<<order<<">"<<std::endl;	
		//		}
		//		
		//		/******************************************************/
		//		// operator*
		//		Simplex<order+1> operator*(const Simplex<0>& other){
		//			return Simplex<order+1>(*this, other);
		//			//return boost::shared_ptr<Simplex<order+1> >
		//			//this->psn.reset(new SubNetworkType(this->p_derived()));
		//		}
		//		
		//		/******************************************************/
		//		// vertex
		//		const Simplex<0>& vertex(const int& k) const {
		//			/*! A const reference to the k-th vertex in the simplex. 
		//			 *  Vertices are of type Simplex<0>.
		//			 */
		//			assert(k>=0 && k<Nvertex);
		//			return (k==Nvertex-1)? last : lower.vertex(k);
		//		}
		
		
		
	};
	
	
	
	
	/******************************************************/
	/* Simplex<order>                                       */
	/******************************************************/
	template<short unsigned int cols>
	struct SortSimplexID {

//		typedef std::set<size_t> SetType;
//		typedef Simplex<dim,order> SimplexType;
//		typedef SimplexType::SimplexIDType SimplexIDType;
		
//		template<short unsigned int dim>
//		static Eigen::Matrix<unsigned int,1,1> getID(const Simplex<dim,0>& S0){
//			return (Eigen::Matrix<unsigned int,1,1>()<<S0.sID).finished();
//		}

		
		typedef Eigen::Matrix<unsigned int,1,cols> VectorType;
		
		static VectorType getID(const VectorType& V){
			std::set<size_t> set;
			for (int k=0;k<cols;++k){
				assert(set.insert(V(k)).second && "COULD NOT INSERT SIMPLEX ID IN SET FOR SORTING.");
			}
			
//			assert(set.insert(S0.sID).second && "COULD NOT INSERT SIMPLEX ID IN SET FOR SORTING.");
//			assert(set.insert(S1.sID).second && "COULD NOT INSERT SIMPLEX ID IN SET FOR SORTING.");
			Eigen::Matrix<unsigned int,1,cols> temp;
			int k(0);
			for (std::set<size_t>::const_iterator iter=set.begin();iter!=set.end();++iter){
				temp(k)=(*iter);
				++k;
			}
			return temp;
		}
		
//		template<short unsigned int dim>
//		static Eigen::Matrix<unsigned int,1,3> getID(const Simplex<dim,0>& S0,
//													 const Simplex<dim,0>& S1,
//													 const Simplex<dim,0>& S2){
//			std::set<size_t> set;
//			assert(set.insert(S0.sID).second && "COULD NOT INSERT SIMPLEX ID IN SET FOR SORTING.");
//			assert(set.insert(S1.sID).second && "COULD NOT INSERT SIMPLEX ID IN SET FOR SORTING.");
//			assert(set.insert(S2.sID).second && "COULD NOT INSERT SIMPLEX ID IN SET FOR SORTING.");
//			Eigen::Matrix<unsigned int,1,3> temp;
//			int k(0);
//			for (std::set<size_t>::const_iterator iter=set.begin();iter!=set.end();++iter){
//				temp(k)=(*iter);
//				++k;
//			}
//			return temp;
//		}
//		
//		template<short unsigned int dim>
//		static Eigen::Matrix<unsigned int,1,4> getID(const Simplex<dim,0>& S0,
//													 const Simplex<dim,0>& S1,
//													 const Simplex<dim,0>& S2,
//													 const Simplex<dim,0>& S3){
//			std::set<size_t> set;
//			assert(set.insert(S0.sID).second && "COULD NOT INSERT SIMPLEX ID IN SET FOR SORTING.");
//			assert(set.insert(S1.sID).second && "COULD NOT INSERT SIMPLEX ID IN SET FOR SORTING.");
//			assert(set.insert(S2.sID).second && "COULD NOT INSERT SIMPLEX ID IN SET FOR SORTING.");
//			assert(set.insert(S3.sID).second && "COULD NOT INSERT SIMPLEX ID IN SET FOR SORTING.");
//			Eigen::Matrix<unsigned int,1,4> temp;
//			int k(0);
//			for (std::set<size_t>::const_iterator iter=set.begin();iter!=set.end();++iter){
//				temp(k)=(*iter);
//				++k;
//			}
//			return temp;
//		}
		
		
	};
	
	
	
	
	
	

	
	
	
	template<short unsigned int dim>
	class Simplex<dim,1> : public ObservedByMap<Eigen::Matrix<unsigned int,1,SimplexEnums<1>::Nnodes>,
	/*                                       */ Simplex<dim,1>,
	/*                                       */ CompareVectorsByComponent<unsigned int,SimplexEnums<1>::Nnodes> >,
	/*                  */ public  ObservingMap<Eigen::Matrix<unsigned int,1,SimplexEnums<0>::Nnodes>,
	/*                                       */ Simplex<dim,0>,
	/*                                       */ CompareVectorsByComponent<unsigned int,SimplexEnums<0>::Nnodes> >{
	
	
		typedef ObservedByMap<Eigen::Matrix<unsigned int,1,SimplexEnums<1>::Nnodes>,
		/*                                       */ Simplex<dim,1>,
		/*                                       */ CompareVectorsByComponent<unsigned int,SimplexEnums<1>::Nnodes> > ObservedType;
		
		typedef ObservingMap<Eigen::Matrix<unsigned int,1,SimplexEnums<0>::Nnodes>,
		/*                                       */ Simplex<dim,0>,
		/*                                       */ CompareVectorsByComponent<unsigned int,SimplexEnums<0>::Nnodes> > ObservingMap0;
		
		// number of unordered 1-sets in a pool of 2
		boost::shared_ptr<Simplex<dim,0> > V0;
		boost::shared_ptr<Simplex<dim,0> > V1;
	
	
		boost::shared_ptr<Simplex<dim,0> > checkS0(const unsigned int k){
		// NOOOOOOOOO THIS CANNOT BE DONE BECAUSE SHERE POINTERS CAN ONLY BE ASSIGNED FROM SHARE POINTERS
		
		}
		
		
	public:
		Simplex(const Eigen::Matrix<unsigned int,1,2> IDs) : ObservedType::ObservedByMap(SortSimplexID<2>::getID(IDs)){
		
		}
	
	}; 
	
	
//	template<short unsigned int dim>
//	class Simplex<dim,2>{
//		
//		// number of unordered 1-sets in a pool of 3		
//		boost::shared_ptr<Simplex<dim,0> > V0;
//		boost::shared_ptr<Simplex<dim,0> > V1;
//		boost::shared_ptr<Simplex<dim,0> > V2;
//		
//		// number of unordered 2-sets in a pool of 3
//		boost::shared_ptr<Simplex<dim,1> > E1; // (V0*V1)
//		boost::shared_ptr<Simplex<dim,1> > E2; // (V1*V2)
//		boost::shared_ptr<Simplex<dim,1> > E3; // (V2*V3)
//		
//		
//	}; 
//	
//	
//	template<short unsigned int dim>
//	class Simplex<dim,3>{
//		
//		// number of unordered 1-sets in a pool of 4
//		boost::shared_ptr<Simplex<dim,0> > V0;
//		boost::shared_ptr<Simplex<dim,0> > V1;
//		boost::shared_ptr<Simplex<dim,0> > V2;
//		boost::shared_ptr<Simplex<dim,0> > V3;
//
//		// number of unordered 2-sets in a pool of 4
//		boost::shared_ptr<Simplex<dim,1> > E01; // (V0*V1)
//		boost::shared_ptr<Simplex<dim,1> > E02; // (V0*V2)
//		boost::shared_ptr<Simplex<dim,1> > E03; // (V0*V3)
//		boost::shared_ptr<Simplex<dim,1> > E12; // (V1*V2)
//		boost::shared_ptr<Simplex<dim,1> > E13; // (V1*V3)
//		boost::shared_ptr<Simplex<dim,1> > E23; // (V2*V3)		
//		
//		// number of unordered 3-sets in a pool of 4
//		
//	}; 
	
	
	
	
	
	/******************************************************/
	/* Simplex: template specialization Simplex<0>        */
	/******************************************************/	
	template<short unsigned int dim>
	class Simplex<dim,0> : public StaticID<Simplex<dim,0> >,
	/*                  */ public ObservedByMap<Eigen::Matrix<unsigned int,1,SimplexEnums<0>::Nnodes>,
	/*                                       */ Simplex<dim,0>,
	/*                                       */ CompareVectorsByComponent<unsigned int,SimplexEnums<0>::Nnodes> >,
	/*                  */ public ObservingMap<Eigen::Matrix<unsigned int,1,SimplexEnums<1>::Nnodes>,
	/*                                       */ Simplex<dim,1>,
	/*                                       */ CompareVectorsByComponent<unsigned int,SimplexEnums<1>::Nnodes> >{
		
		enum {order=0};
		
		typedef Eigen::Matrix<unsigned int,1,SimplexEnums<order>::Nnodes> SimplexIDType;
		typedef CompareVectorsByComponent<unsigned int,SimplexEnums<order>::Nnodes> CompareSimplexIDType;
		typedef ObservedByMap<SimplexIDType,Simplex<dim,order>,CompareSimplexIDType> BaseObservedType;
		
		typedef Eigen::Matrix<unsigned int,1,SimplexEnums<order+1>::Nnodes> NextSimplexIDType;
		typedef CompareVectorsByComponent<unsigned int,SimplexEnums<order+1>::Nnodes> CompareNextSimplexIDType;
		typedef ObservingMap<NextSimplexIDType,Simplex<dim,order+1>,CompareNextSimplexIDType> BaseObservingType;
		
	public:
		//enum {order=1};
//		enum {Nnodes=2};
//		enum {Nedges=1};
		
		/* Constructor: this sould take a position ****************************/
		Simplex() : BaseObservedType::ObservedByMap((SimplexIDType()<<this->sID).finished()){
			std::cout<<"Creating new Simplex<0> with ID "<<this->sID<<std::endl;
		}
		
		/* Destructor *********************************************************/
		~Simplex(){
			std::cout<<"Destroying Simplex<0> with ID "<<this->sID<<std::endl;
		}
		
//		/* operator * *********************************************************/
//		void operator*(const Simplex<dim,0>& other){
////			if ()
//			
////			return Simplex<1>(*this,other);	
//		}
//		
//		
//		// operator ()
//		const Simplex<0>& vertex(const int& k) const {
//			assert(k==0);
//			return *this;
//		}
	};
	

	
	
	
	
}
#endif
