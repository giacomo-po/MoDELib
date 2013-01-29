/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_NETWORKNODE_H_
#define model_NETWORKNODE_H_

#ifndef VERBOSELEVEL
#define VERBOSELEVEL 0
#endif

#include <iomanip>
#include <map>
#include <set>
#include <assert.h>
#include <algorithm>
#include <limits.h>

#include <boost/ptr_container/ptr_map.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/utility.hpp>

#include <model/Network/Operations/includeNetworkOperations.h>
#include <model/Utilities/StaticID.h>
#include <model/Utilities/CRTP.h>


namespace model {
	
	template <typename Derived>
	class NetworkNode : boost::noncopyable,
	/*               */ public CRTP<Derived>,		
	/*               */ public StaticID<Derived>{		
		
	public:
#include <model/Network/NetworkTypedefs.h>
	private:
//#include "model/Network/SubNetworkComponent.h"
		
		
		boost::shared_ptr<SubNetworkType> psn;

		
	protected:
		
		NeighborContainerType Neighborhood;			
		NeighborContainerType OutNeighborhood;
		NeighborContainerType InNeighborhood;
		
	public:
		
		/*****************************************************************************************/
		/* Costructor with node arguments ********************************************************/
		NetworkNode() : psn(new SubNetworkType(this->p_derived())){		
#if VERBOSELEVEL == 0
			std::cout<<"Creating Node "<<this->sID<<std::endl;
#endif
			// Insert this->p_derived() in the Neighborhood
			Neighborhood.insert(std::make_pair(this->sID, boost::tuples::make_tuple(this->p_derived(),(LinkType*) NULL,0) ));
			
		}
		
		/*****************************************************************************************/
		/* Costructor from EdgeExpansion *********************************************************/
		//		NetworkNode(const EdgeExpansion<LinkType>& ee) : state(0),
		NetworkNode(const ExpandingEdge<LinkType>& ee) : psn(ee.E.pSN()){		
#if VERBOSELEVEL == 0
			std::cout<<"Creating Node "<<this->sID<<" by EdgeExpansion."<<std::endl;
#endif
			// Insert this->p_derived() in the Neighborhood
			Neighborhood.insert(std::make_pair(this->sID, boost::tuples::make_tuple(this->p_derived(),(LinkType*) NULL,0) ));
			
			// Manage SubNetwork
			psn->add(this->p_derived());
		}
				
		/*****************************************************************************************/
		/* Destructor ****************************************************************************/
		~NetworkNode(){
#if VERBOSELEVEL == 0
			std::cout<<"Deleting Node "<<this->sID<<std::endl;
#endif
			//! 1- Remove this from Neighborhood	
			Neighborhood.erase(this->sID);
			
			//! 2- Remove this from the SubNetwork	
			this->psn->remove(this->p_derived());	// remove this in the new subnetwork
		}
		
		//////////////////////////////////////////////////////////////////////////////
		// formSubNetwork
		void formSubNetwork(const boost::shared_ptr<SubNetworkType> & psnOther){
			if (psn!=psnOther){
				psn->remove(this->p_derived());
				psn=psnOther;		// redirect psn to the new Subnetwork
				psn->add(this->p_derived());	// add this in the new subnetwork
			}
		}
		
		
		
		//////////////////////////////////////////////////////////////////////////////
		// sndID
		size_t snID() const {
			return psn->snID(this->p_derived());
		}
		
		
		
		//////////////////////////////////////////////////////////////////////////////
		//! Returns a const pointer to the parent SubNetwork
		const boost::shared_ptr<SubNetworkType> & pSN() const {
			return psn;
		}
		
		
		/*****************************************************************************************/
		/* addToNeighborhood *********************************************************************/
		void addToNeighborhood(LinkType* const pL){
			
			Derived* pN=NULL;
			size_t key=0;
			short int dir;
			
			if (pL->source==this->p_derived()){
				pN=pL->sink;
				key=pN->sID;
				dir=1;
				assert(OutNeighborhood.insert(std::make_pair(key, boost::tuples::make_tuple(pN,pL,dir) )).second && "CANNOT INSERT IN OUT_NEIGHBORHOOD");
			}
			
			if (pL->sink==this->p_derived()){
				pN=pL->source;
				key=pN->sID;
				dir=-1;
				assert(InNeighborhood.insert(std::make_pair(key, boost::tuples::make_tuple(pN,pL,dir) )).second && "CANNOT INSERT IN IN_NEIGHBORHOOD");
			}
			
			//			bool success=Neighborhood.insert(std::make_pair(key, boost::tuples::make_tuple(pN,pL,dir) )).second;
			//			assert(success);
			assert(Neighborhood.insert(std::make_pair(key, boost::tuples::make_tuple(pN,pL,dir) )).second && "CANNOT INSERT IN NEIGHBORHOOD.");
			
		}		
		
		/*****************************************************************************************/
		/* removeFromNeighborhood ****************************************************************/
		void removeFromNeighborhood(LinkType* const pL){
			
			Derived* pN=NULL;
			size_t key=0;
			
			if (pL->source==this->p_derived()){
				pN=pL->sink;
				key=pN->sID;
				OutNeighborhood.erase(key);
			}
			
			if (pL->sink==this->p_derived()){
				pN=pL->source;
				key=pN->sID;
				InNeighborhood.erase(key);
			}
			
			bool success=Neighborhood.erase(key);
			assert(success);
			
		}	
		
		/*****************************************************************************************/
		/* outFlow *******************************************************************************/
        FlowType outFlow() const {
//            FlowType Fout;
            FlowType Fout(FlowType::Zero()); // generalize
			Fout*=0.0;
			for (typename NeighborContainerType::const_iterator     NeighborIter=OutNeighborhood.begin();NeighborIter!=OutNeighborhood.end();++NeighborIter){
				Fout+=boost::tuples::get<1>(NeighborIter->second)->flow;
			}
            return Fout;
        }
        
		/*****************************************************************************************/
		/* inFlow ********************************************************************************/
        FlowType inFlow() const {
//            FlowType Fin;
            FlowType Fin(FlowType::Zero());
			Fin*=0.0;
			for (typename NeighborContainerType::const_iterator     NeighborIter=InNeighborhood.begin();NeighborIter!=InNeighborhood.end();++NeighborIter){
				Fin+=boost::tuples::get<1>(NeighborIter->second)->flow;
			}
            return Fin;
        }

		/*****************************************************************************************/
		/* depthFirstSearch **********************************************************************/
		bool depthFirstSearch (const size_t& ID, const size_t& N = ULONG_MAX) const {
			std::set<size_t> searchedNodes;
			return depthFirstSearch(searchedNodes,ID,N);
		}
		
		bool depthFirstSearch (std::set<size_t>& searchedNodes, const size_t& ID, const size_t& N = ULONG_MAX) const{
			bool reached(this->sID==ID);
			if (N!=0 && !reached){
				assert(searchedNodes.insert(this->sID).second && "CANNOT INSERT CURRENT NODE IN SEARCHED NODES"); // this node has been searched
				for (typename NeighborContainerType::const_iterator NeighborIter=Neighborhood.begin();NeighborIter!=Neighborhood.end();++NeighborIter){
					if (searchedNodes.find(boost::tuples::get<0>(NeighborIter->second)->sID)==searchedNodes.end()){  // neighbor not searched
						reached=boost::tuples::get<0>(NeighborIter->second)->depthFirstSearch(searchedNodes,ID,N-1); // ask if neighbor can reach
						if (reached){
							break;
						}
					}
				}	
			}
			return reached;
		}
		
		/*****************************************************************************************/
		/* depthFirstExecute *********************************************************************/
		template <typename T>
		void depthFirstExecute(void (Derived::*Nfptr)(const T&),void (LinkType::*Lfptr)(const T&), const T & input, const size_t& N = ULONG_MAX){
			std::set<size_t> searchedNodes;
			std::set<std::pair<size_t,size_t> > searchedLinks;
			depthFirstExecute(searchedNodes,searchedLinks,Nfptr,Lfptr,input, N);
		}
		
		template <typename T>
		void depthFirstExecute(std::set<size_t>& searchedNodes, std::set<std::pair<size_t,size_t> >& searchedLinks, void (Derived::*Nfptr)(const T&),void (LinkType::*Lfptr)(const T&), const T & input, const size_t& N = ULONG_MAX){
			(this->p_derived()->*Nfptr)(input); // execute Nfptr on this node
			if (N!=0){
				assert(searchedNodes.insert(this->sID).second && "CANNOT INSERT CURRENT NODE IN SEARCHED NODES"); // this node has been searched
				for (typename NeighborContainerType::iterator NeighborIter=Neighborhood.begin();NeighborIter!=Neighborhood.end();++NeighborIter){
					if (!boost::tuples::get<2>(NeighborIter->second)==0){
						if (searchedLinks.find(boost::tuples::get<1>(NeighborIter->second)->nodeIDPair)==searchedLinks.end()){  // neighbor not searched
							(boost::tuples::get<1>(NeighborIter->second)->*Lfptr)(input); // execute Lfptr on connecting link
							assert(searchedLinks.insert(boost::tuples::get<1>(NeighborIter->second)->nodeIDPair).second && "CANNOT INSERT CURRENT LINK IN SEARCHED LINKS"); // this node has been searched
						}
					}
					if (searchedNodes.find(boost::tuples::get<0>(NeighborIter->second)->sID)==searchedNodes.end()){  // neighbor not searched
						boost::tuples::get<0>(NeighborIter->second)->depthFirstExecute(searchedNodes,searchedLinks,Nfptr,Lfptr,input, N-1); // continue executing on neighbor
					}
				}	
			}
		}
		
		////////////////////////////////////////////////////////
		// transmit<T>
		template <typename T>
		void depthFirstNodeExecute(void (Derived::*Nfptr)(const T&), const T & input, const size_t& N = ULONG_MAX){
			std::set<size_t> searchedNodes;
			depthFirstNodeExecute(searchedNodes,Nfptr,input, N);
		}
		
		template <typename T>
		void depthFirstNodeExecute(std::set<size_t>& searchedNodes, void (Derived::*Nfptr)(const T&), const T & input, const size_t& N = ULONG_MAX){
			(this->p_derived()->*Nfptr)(input); // execute Nfptr on this node
			if (N!=0){
				assert(searchedNodes.insert(this->sID).second && "CANNOT INSERT CURRENT NODE IN SEARCHED NODES"); // this node has been searched
				for (typename NeighborContainerType::iterator NeighborIter=Neighborhood.begin();NeighborIter!=Neighborhood.end();++NeighborIter){
					if (searchedNodes.find(boost::tuples::get<0>(NeighborIter->second)->sID)==searchedNodes.end()){  // neighbor not searched
						boost::tuples::get<0>(NeighborIter->second)->depthFirstNodeExecute(searchedNodes,Nfptr,input, N-1); // continue executing on neighbor
					}
				}	
			}
		}
		
		void depthFirstNodeExecute(void (Derived::*Nfptr)(void), const size_t& N = ULONG_MAX){
			std::set<size_t> searchedNodes;
			depthFirstNodeExecute(searchedNodes,Nfptr, N);
		}
		
		void depthFirstNodeExecute(std::set<size_t>& searchedNodes, void (Derived::*Nfptr)(void), const size_t& N = ULONG_MAX){
			(this->p_derived()->*Nfptr)(); // execute Nfptr on this node
			if (N!=0){
				assert(searchedNodes.insert(this->sID).second && "CANNOT INSERT CURRENT NODE IN SEARCHED NODES"); // this node has been searched
				for (typename NeighborContainerType::iterator NeighborIter=Neighborhood.begin();NeighborIter!=Neighborhood.end();++NeighborIter){
					if (searchedNodes.find(boost::tuples::get<0>(NeighborIter->second)->sID)==searchedNodes.end()){  // neighbor not searched
						boost::tuples::get<0>(NeighborIter->second)->depthFirstNodeExecute(searchedNodes,Nfptr, N-1); // continue executing on neighbor
					}
				}	
			}
		}
		

		
		
		////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////
		void resetPSN(){
			//! 1- Removes this from the current SubNetwork
			this->psn->remove(this->p_derived());
			//! 2- Creates a new SubNetwork containing this
			this->psn.reset(new SubNetworkType(this->p_derived()));		
			//! 3- Transmits 'formSubNetwork' to the neighbors
			typedef void (Derived::*node_member_function_pointer_type)(const boost::shared_ptr<SubNetworkType>&); 
			node_member_function_pointer_type Nmfp(&Derived::formSubNetwork);
//			Nmfp=&Derived::formSubNetwork;
			typedef void (LinkType::*link_member_function_pointer_type)(const boost::shared_ptr<SubNetworkType>&); 
			link_member_function_pointer_type Lmfp(&LinkType::formSubNetwork);
//			Lmfp=&LinkType::formSubNetwork;
			depthFirstExecute(Nmfp,Lmfp,this->psn);
		}
		
		
		////////////////////////////////////////////////////////
		// neighborhood
		const NeighborContainerType& neighborhood() const {
			return Neighborhood;
		}
		
		////////////////////////////////////////////////////////
		// closedNeighborNode
		NodeType* closedNeighborNode(const unsigned int& k) const {
			assert(k<Neighborhood.size() && "INDEX EXCEEDS SIZE");
			typename NeighborContainerType::const_iterator nIter(Neighborhood.begin());
			for (unsigned int n=0;n<k;++n){
				nIter++;
			}
			return (boost::tuples::get<0>(nIter->second));
		}
		
		////////////////////////////////////////////////////////
		// closedNeighborNode
		NodeType* openNeighborNode(const unsigned int& k) const {
			assert(k<(Neighborhood.size()-1) && "INDEX EXCEEDS SIZE");
			typename NeighborContainerType::const_iterator nIter(Neighborhood.begin());
			if(nIter->first == this->sID){
				nIter++;
			}
			for (unsigned int n=0;n<k;++n){
				if(nIter->first == this->sID){
					nIter++;
				}
				nIter++;
			}
			if(nIter->first == this->sID){
				nIter++;
			}
			return (boost::tuples::get<0>(nIter->second));
		}
		
		////////////////////////////////////////////////////////
		// closedNeighborNode
		LinkType* closedNeighborLink(const unsigned int& k) const {
			assert(k<Neighborhood.size() && "INDEX EXCEEDS SIZE");
			typename NeighborContainerType::const_iterator nIter(Neighborhood.begin());
			for (unsigned int n=0;n<k;++n){
				nIter++;
			}
			return (boost::tuples::get<1>(nIter->second));
		}
		
		////////////////////////////////////////////////////////
		// closedNeighborNode
		LinkType* openNeighborLink(const unsigned int& k) const {
			assert(k<(Neighborhood.size()-1) && "INDEX EXCEEDS SIZE");
			typename NeighborContainerType::const_iterator nIter(Neighborhood.begin());
			if(nIter->first == this->sID){
				nIter++;
			}
			for (unsigned int n=0;n<k;++n){
				if(nIter->first == this->sID){
					nIter++;
				}
				nIter++;
			}
			if(nIter->first == this->sID){
				nIter++;
			}
			return (boost::tuples::get<1>(nIter->second));
		}
		
		////////////////////////////////////////////////////////
		// outNeighborhood
		const NeighborContainerType & outNeighborhood() const {
			return OutNeighborhood;
		}
		
		////////////////////////////////////////////////////////
		// inNeighborhood
		const NeighborContainerType & inNeighborhood() const {
			return InNeighborhood;
		}
		
		////////////////////////////////////////////////////////
		// neighborID
		size_t neighborID(const size_t & k) const {
			return std::distance(Neighborhood.begin(),Neighborhood.find(k));
		}
		
        ////////////////////////////////////////////////////////
		// outOrder
		size_t outOrder() const {
			return OutNeighborhood.size();
		}
		
		////////////////////////////////////////////////////////
		// inOrder
		size_t inOrder() const {
			return InNeighborhood.size();
		}
        
        ////////////////////////////////////////////////////////
		// order
		size_t openOrder() const {
			return outOrder()+inOrder();
		}
        
        ////////////////////////////////////////////////////////
		// order
		size_t closedOrder() const {
			return openOrder()+1;
		}
        
		//////////////////////////////////////////////////////////////////////////////
		// is_source() // THIS IS MEANINGLESS IF LINK DIRECTIONS CAN BE SWITCHED ARBITRARILY
		bool is_source() const {
			return inOrder()==0 && outOrder()>0;
		}
		
		//////////////////////////////////////////////////////////////////////////////
		// is_sink() // THIS IS MEANINGLESS IF LINK DIRECTIONS CAN BE SWITCHED ARBITRARILY
		bool is_sink() const {
			return outOrder()==0 && inOrder()>0;
		}
		
		//////////////////////////////////////////////////////////////////////////////
		// is_isolated
		bool is_isolated() const {
			return inOrder()==0 && outOrder()==0;
		}
		
		//////////////////////////////////////////////////////////////////////////////
		// is_central
		bool is_central() const {
			return inOrder()>0 && outOrder()>0;
		}
		
		//////////////////////////////////////////////////////////////////////////////
		// is_balanced
		//		bool is_balanced() const {
		//			return outFlow() == inFlow();
		//		}
		
		bool is_balanced() const {
			//			return FlowCompare<FlowType>(outFlow(),inFlow());
			return (outFlow() - inFlow()).norm()<FLT_EPSILON;
		}
		
		
		
		
		
		
		
		//////////////////////////////////////////////////////////////////////////////
		// is_simple
		bool is_simple() const {
			return (outOrder()+inOrder())==2 && is_balanced();
		}
		
		////////////////////////////////////////////////////////
		// friend T& operator <<
		template <class T, typename OtherDerived>
		friend T& operator << (T& os, const NetworkNode<OtherDerived> & NN){
			
			os << "Node sID=" << std::setw(3)<<NN.sID
			<< " dID=" << std::setw(3)<<NN.dID()
			<< " snID=" << std::setw(3)<<NN.snID()<<" ("
			<< std::setw(3)<<NN.dID() <<") ["<<NN.p_derived()<<"] : ["<<NN.address(NN.sID)<<"]"<<
			" state = "<<NN.get_state()<<std::endl;
			
			os << "		OutLinks: (outOrder= "<<NN.outOrder()<<")"<<std::endl;
			for(typename NeighborContainerType::const_iterator	 NeighborIter=NN.OutNeighborhood.begin();NeighborIter!=NN.OutNeighborhood.end();++NeighborIter){
				os<< *boost::tuples::get<1>(NeighborIter->second)<<std::endl;//
			}
			
			os << "		InLinks: (inOrder= "<<NN.inOrder()<<")"<<std::endl;
			for(typename NeighborContainerType::const_iterator	 NeighborIter=NN.InNeighborhood.begin();NeighborIter!=NN.InNeighborhood.end();++NeighborIter){
				os<< *boost::tuples::get<1>(NeighborIter->second)<<std::endl;//
			}
			
			return os;
		}
		
	};
		
	
	//////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////
} // namespace model
#endif





//		////////////////////////////////////////////////////////
//		////////////////////////////////////////////////////////
//		// nodeTransmit
//		template <typename T>
//		void nodeTransmit(void (Derived::*Nfptr)(const T&), const T & input, const size_t& N = ULONG_MAX){
//			if (!this->get_state() && N!=0){
//				this->set_state();
//				(this->p_derived()->*Nfptr)(input);
//				for (typename NeighborContainerType::iterator NeighborIter=Neighborhood.begin();NeighborIter!=Neighborhood.end();++NeighborIter){
//					if (!boost::tuples::get<2>(NeighborIter->second)==0){
//						boost::tuples::get<0>(NeighborIter->second)->nodeTransmit(Nfptr,input, N-1);//
//					}
//				}	
//				this->reset_state();				
//			}
//		}
//		
//		void nodeTransmit(void (Derived::*Nfptr)(void), const size_t& N = ULONG_MAX){
//			if (!this->get_state() && N!=0){
//				this->set_state();
//				(this->p_derived()->*Nfptr)();
//				for (typename NeighborContainerType::iterator NeighborIter=Neighborhood.begin();NeighborIter!=Neighborhood.end();++NeighborIter){
//					if (!boost::tuples::get<2>(NeighborIter->second)==0){
//						boost::tuples::get<0>(NeighborIter->second)->nodeTransmit(Nfptr, N-1);//
//					}
//				}	
//				this->reset_state();				
//			}
//		}
//		
//		
//		void linkTransmit(void (LinkType::*Lfptr)(void), const size_t& N = ULONG_MAX){
//			if (!this->get_state() && N!=0){
//				this->set_state();
//				//(this->p_derived()->*Nfptr)();
//				for (typename NeighborContainerType::iterator NeighborIter=Neighborhood.begin();NeighborIter!=Neighborhood.end();++NeighborIter){
//					if (!boost::tuples::get<2>(NeighborIter->second)==0){
//						(boost::tuples::get<1>(NeighborIter->second)->*Lfptr)();
//						boost::tuples::get<0>(NeighborIter->second)->linkTransmit(Lfptr, N-1);//
//					}
//				}	
//				this->reset_state();				
//			}
//		}



//		//////////////////////////////////////////////////////////////////////////////
//		// is_simple
//		bool is_simple() const {
//			return outOrder()==1 && inOrder()==1 && is_balanced();
//		}






//		////////////////////////////////////////////////////////
//		// canReach	
//		bool canReach(const size_t& ID, const size_t& N = ULONG_MAX){
//			bool reached(this->sID==ID);
//			if (!this->get_state() && N!=0 && !reached){
//				std::cout<<"Node "<<this->sID<< "can reach "<<ID;
//				this->set_state();
//				for (typename NeighborContainerType::const_iterator NeighborIter=Neighborhood.begin();NeighborIter!=Neighborhood.end();++NeighborIter){
//					if (!boost::tuples::get<2>(NeighborIter->second)==0){
//						reached=boost::tuples::get<0>(NeighborIter->second)->canReach(ID,N-1);//
//						if (reached){
//							std::cout<<"yes, reached"<<std::endl;
////							this->reset_state();				
//							break;
//						}
//						
//					}
//				}	
//				this->reset_state();	
//				std::cout<<"? "<<reached<<std::endl;;
//			}
//
//			return reached;
//		}





//		////////////////////////////////////////////////////////
//		// canOutReach	
//		bool canOutReach(const size_t& ID, const size_t& N = ULONG_MAX){
//			bool reached(this->sID==ID);
//			if (!this->get_state() && N!=0 && !reached){
//				this->set_state();
//				for (typename NeighborContainerType::const_iterator NeighborIter=OutNeighborhood.begin();NeighborIter!=OutNeighborhood.end();++NeighborIter){
//					reached=boost::tuples::get<0>(NeighborIter->second)->canReach(ID,N-1);//
//					if (reached){
//						break;
//					}
//				}	
//				this->reset_state();				
//			}
//			return reached;
//		}


//		////////////////////////////////////////////////////////
//		// pathTo
//		std::vector<size_t> pathTo(const size_t& ID){
//			std::vector<size_t> temp;
//			if (!this->get_state()){
//				this->set_state();
//				bool reached(this->sID==ID);
//				if (!reached) {
//					for (typename NeighborContainerType::const_iterator NeighborIter=OutNeighborhood.begin();NeighborIter!=OutNeighborhood.end();++NeighborIter){
//						temp=boost::tuples::get<0>(NeighborIter->second)->pathTo(ID);//
//						if (temp.size()){ // there is a path from the current neighbor
//							temp.push_back(this->sID); // add this to the path
//							break;
//						}
//					}
//				}
//				else{
//					temp.push_back(this->sID);
//				}
//				this->reset_state();				
//			}
//			return temp;
//		}




//		////////////////////////////////////////////////////////
//		// flow THIS OPERATION IS NOT WELL DEFINED IF THE NODE IS UNBALANCED
//		FlowType flow() {
//			FlowType F;
//			//F-=F; // this should make F "0" independently of the initial value
//			if ( is_central() ){
//				assert( is_balanced() );
//				F=outFlow();
//			}
//			else if (is_sink()){
//				F=inFlow();
//			}
//			else if (is_source() ){
//				F=outFlow();
//			}
//			else{
//				assert(0);
//			}
//			
//			return F;
//			
//			// = outFlow()*(is_source()+is_central())+ inFlow()*is_sink()
//		}




//		/*****************************************************************************************/
//		/* initData ******************************************************************************/
//		const InitializationData<Derived>& initData() const {
//			/*! The InitializationData<Derived> object used to construct this Vertex.
//			 */
//			return *static_cast<const InitializationData<Derived>*> (this);
//		}

////////////////////////////////////////////////////////
// flowBalance
//		void flowBalance() const {
//			return outFlow()-inFlow();
//		}

//		//////////////////////////////////////////////////
//		// topologyChangeActions
//		void topologyChangeActions(){
//			/*! This function is called by the NeworkLink destructor. It sould be overwritten by
//			 *  the Derived class if some action is to be taken when topological changes in connection take place.
//			 */
//		}




//		////////////////////////////////////////////////////////
//		// outFlow
//		FlowType outFlow() const {
//			FlowType Fout;
//			Fout-=Fout; // this should make Fout "0" independently of the initial value
//			for (typename NeighborContainerType::const_iterator	 NeighborIter=OutNeighborhood.begin();NeighborIter!=OutNeighborhood.end();++NeighborIter){
//				Fout+=boost::tuples::get<1>(NeighborIter->second)->flow;
//			}
//			return Fout;
//		}
//		
//		////////////////////////////////////////////////////////
//		// inFlow
//		FlowType inFlow() const {
//			FlowType Fin;
//			Fin-=Fin;	// this should make Fin "0" independently of the initial value
//			for (typename NeighborContainerType::const_iterator	 NeighborIter=InNeighborhood.begin();NeighborIter!=InNeighborhood.end();++NeighborIter){
//				Fin+=boost::tuples::get<1>(NeighborIter->second)->flow;
//			}
//			return Fin;
//		}
