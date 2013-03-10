/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_NETWORKNODE_H_
#define model_NETWORKNODE_H_

#include <iomanip>
#include <map>
#include <set>
#include <assert.h>
#include <algorithm>
#include <limits.h>

#include <boost/ptr_container/ptr_map.hpp>
//#include <boost/tuple/tuple.hpp>
#include <tuple> // std::tuple replaces boost::tuple in c++11
//#include <boost/shared_ptr.hpp>
//#include <boost/utility.hpp>
#include <memory> // std::shared_ptr


#include <model/Utilities/StaticID.h>
#include <model/Utilities/CRTP.h>
//#include "model/Network/NetworkLink.h"
#include <model/Network/Operations/includeNetworkOperations.h>


namespace model {
    
    template <typename Derived>
	class NetworkLink; // class predeclaration
	
	template <typename Derived>
	class NetworkNode : boost::noncopyable,
	/*               */ public CRTP<Derived>,		
	/*               */ public StaticID<Derived>{		
		
	public:
#include <model/Network/NetworkTypedefs.h>
        friend class NetworkLink<LinkType>; // allow NetworkLink to call private NetworkNode::formSubNetwork
        
        
	private:
//#include "model/Network/SubNetworkComponent.h"
		
		
		std::shared_ptr<SubNetworkType> psn;

   	protected:
     
        /**********************************************************************/
		void resetPSN(){
			//! 1- Removes this from the current SubNetwork
			this->psn->remove(this->p_derived());
			//! 2- Creates a new SubNetwork containing this
			this->psn.reset(new SubNetworkType(this->p_derived()));
			//! 3- Transmits 'formSubNetwork' to the neighbors
			typedef void (Derived::*node_member_function_pointer_type)(const std::shared_ptr<SubNetworkType>&);
			node_member_function_pointer_type Nmfp(&Derived::formSubNetwork);
            //			Nmfp=&Derived::formSubNetwork;
			typedef void (LinkType::*link_member_function_pointer_type)(const std::shared_ptr<SubNetworkType>&);
			link_member_function_pointer_type Lmfp(&LinkType::formSubNetwork);
            //			Lmfp=&LinkType::formSubNetwork;
			depthFirstExecute(Nmfp,Lmfp,this->psn);
		}
        
        /**********************************************************************/
        void formSubNetwork(const std::shared_ptr<SubNetworkType> & psnOther)
        {
			if (psn!=psnOther){
				psn->remove(this->p_derived());
				psn=psnOther;		// redirect psn to the new Subnetwork
				psn->add(this->p_derived());	// add this in the new subnetwork
			}
		}
        
		/* addToNeighborhood **************************************************/
		void addToNeighborhood(LinkType* const pL){
			
			Derived* pN=NULL;
			size_t key=0;
			short int dir;
			
			if (pL->source==this->p_derived()){
				pN=pL->sink;
				key=pN->sID;
				dir=1;
				assert(OutNeighborhood.insert(std::make_pair(key, std::make_tuple(pN,pL,dir) )).second && "CANNOT INSERT IN OUT_NEIGHBORHOOD");
			}
			
			if (pL->sink==this->p_derived()){
				pN=pL->source;
				key=pN->sID;
				dir=-1;
				assert(InNeighborhood.insert(std::make_pair(key, std::make_tuple(pN,pL,dir) )).second && "CANNOT INSERT IN IN_NEIGHBORHOOD");
			}
			
			assert(Neighborhood.insert(std::make_pair(key, std::make_tuple(pN,pL,dir) )).second && "CANNOT INSERT IN NEIGHBORHOOD.");
			
		}
		
		/* removeFromNeighborhood *********************************************/
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
		
		
		NeighborContainerType Neighborhood;			
		NeighborContainerType OutNeighborhood;
		NeighborContainerType InNeighborhood;
		
	public:
		
		/* Costructor with node arguments *************************************/
		NetworkNode() : psn(new SubNetworkType(this->p_derived())){		
			// Insert this->p_derived() in the Neighborhood
			Neighborhood.insert(std::make_pair(this->sID, std::make_tuple(this->p_derived(),(LinkType*) NULL,0) ));
			
		}
		
		/* Costructor from EdgeExpansion **************************************/
		//		NetworkNode(const EdgeExpansion<LinkType>& ee) : state(0),
		NetworkNode(const ExpandingEdge<LinkType>& ee) : psn(ee.E.pSN()){		
			// Insert this->p_derived() in the Neighborhood
			Neighborhood.insert(std::make_pair(this->sID, std::make_tuple(this->p_derived(),(LinkType*) NULL,0) ));
			
			// Manage SubNetwork
			psn->add(this->p_derived());
		}
				
		/* Destructor *********************************************************/
		~NetworkNode(){
			//! 1- Remove this from Neighborhood	
			Neighborhood.erase(this->sID);
			
			//! 2- Remove this from the SubNetwork	
			this->psn->remove(this->p_derived());	// remove this in the new subnetwork
		}
		

		
		
		
		//////////////////////////////////////////////////////////////////////////////
		// sndID
		size_t snID() const {
			return psn->snID(this->p_derived());
		}
		
		
		
		//////////////////////////////////////////////////////////////////////////////
		//! Returns a const pointer to the parent SubNetwork
		const std::shared_ptr<SubNetworkType> & pSN() const {
			return psn;
		}
		
		

		
		/*****************************************************************************************/
		/* outFlow *******************************************************************************/
        FlowType outFlow() const {
//            FlowType Fout;
            FlowType Fout(FlowType::Zero()); // generalize
			Fout*=0.0;
			for (typename NeighborContainerType::const_iterator     NeighborIter=OutNeighborhood.begin();NeighborIter!=OutNeighborhood.end();++NeighborIter){
				Fout+=std::get<1>(NeighborIter->second)->flow;
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
				Fin+=std::get<1>(NeighborIter->second)->flow;
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
					if (searchedNodes.find(std::get<0>(NeighborIter->second)->sID)==searchedNodes.end()){  // neighbor not searched
						reached=std::get<0>(NeighborIter->second)->depthFirstSearch(searchedNodes,ID,N-1); // ask if neighbor can reach
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
					if (!std::get<2>(NeighborIter->second)==0){
						if (searchedLinks.find(std::get<1>(NeighborIter->second)->nodeIDPair)==searchedLinks.end()){  // neighbor not searched
							(std::get<1>(NeighborIter->second)->*Lfptr)(input); // execute Lfptr on connecting link
							assert(searchedLinks.insert(std::get<1>(NeighborIter->second)->nodeIDPair).second && "CANNOT INSERT CURRENT LINK IN SEARCHED LINKS"); // this node has been searched
						}
					}
					if (searchedNodes.find(std::get<0>(NeighborIter->second)->sID)==searchedNodes.end()){  // neighbor not searched
						std::get<0>(NeighborIter->second)->depthFirstExecute(searchedNodes,searchedLinks,Nfptr,Lfptr,input, N-1); // continue executing on neighbor
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
					if (searchedNodes.find(std::get<0>(NeighborIter->second)->sID)==searchedNodes.end()){  // neighbor not searched
						std::get<0>(NeighborIter->second)->depthFirstNodeExecute(searchedNodes,Nfptr,input, N-1); // continue executing on neighbor
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
					if (searchedNodes.find(std::get<0>(NeighborIter->second)->sID)==searchedNodes.end()){  // neighbor not searched
						std::get<0>(NeighborIter->second)->depthFirstNodeExecute(searchedNodes,Nfptr, N-1); // continue executing on neighbor
					}
				}	
			}
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
			return (std::get<0>(nIter->second));
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
			return (std::get<0>(nIter->second));
		}
		
		////////////////////////////////////////////////////////
		// closedNeighborNode
		LinkType* closedNeighborLink(const unsigned int& k) const {
			assert(k<Neighborhood.size() && "INDEX EXCEEDS SIZE");
			typename NeighborContainerType::const_iterator nIter(Neighborhood.begin());
			for (unsigned int n=0;n<k;++n){
				nIter++;
			}
			return (std::get<1>(nIter->second));
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
			return (std::get<1>(nIter->second));
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
				os<< *std::get<1>(NeighborIter->second)<<std::endl;//
			}
			
			os << "		InLinks: (inOrder= "<<NN.inOrder()<<")"<<std::endl;
			for(typename NeighborContainerType::const_iterator	 NeighborIter=NN.InNeighborhood.begin();NeighborIter!=NN.InNeighborhood.end();++NeighborIter){
				os<< *std::get<1>(NeighborIter->second)<<std::endl;//
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
            //		void resetPSN(){
            //			//! 1- Removes this from the current SubNetwork
            //			this->psn->remove(this->p_derived());
            //			//! 2- Creates a new SubNetwork containing this
            //			this->psn.reset(new SubNetworkType(this->p_derived()));
            //			//! 3- Transmits 'formSubNetwork' to the neighbors
            //			typedef void (Derived::*node_member_function_pointer_type)(const std::shared_ptr<SubNetworkType>&);
            //			node_member_function_pointer_type Nmfp(&Derived::formSubNetwork);
            ////			Nmfp=&Derived::formSubNetwork;
            //			typedef void (LinkType::*link_member_function_pointer_type)(const std::shared_ptr<SubNetworkType>&);
            //			link_member_function_pointer_type Lmfp(&LinkType::formSubNetwork);
            ////			Lmfp=&LinkType::formSubNetwork;
            //			depthFirstExecute(Nmfp,Lmfp,this->psn);
            //		}
            
            //////////////////////////////////////////////////////////////////////////////
            // formSubNetwork
            //		void formSubNetwork(const std::shared_ptr<SubNetworkType> & psnOther){
            //			if (psn!=psnOther){
            //				psn->remove(this->p_derived());
            //				psn=psnOther;		// redirect psn to the new Subnetwork
            //				psn->add(this->p_derived());	// add this in the new subnetwork
            //			}
            //		}
            
            //////////////////////////////////////////////////////////////////////////////
            // is_balanced
            //		bool is_balanced() const {
            //			return outFlow() == inFlow();
            //		}
            
            
            //		/*****************************************************************************************/
            //		/* addToNeighborhood *********************************************************************/
            //		void addToNeighborhood(LinkType* const pL){
            //
            //			Derived* pN=NULL;
            //			size_t key=0;
            //			short int dir;
            //
            //			if (pL->source==this->p_derived()){
            //				pN=pL->sink;
            //				key=pN->sID;
            //				dir=1;
            //				assert(OutNeighborhood.insert(std::make_pair(key, std::make_tuple(pN,pL,dir) )).second && "CANNOT INSERT IN OUT_NEIGHBORHOOD");
            //			}
            //
            //			if (pL->sink==this->p_derived()){
            //				pN=pL->source;
            //				key=pN->sID;
            //				dir=-1;
            //				assert(InNeighborhood.insert(std::make_pair(key, std::make_tuple(pN,pL,dir) )).second && "CANNOT INSERT IN IN_NEIGHBORHOOD");
            //			}
            //
            //			//			bool success=Neighborhood.insert(std::make_pair(key, std::make_tuple(pN,pL,dir) )).second;
            //			//			assert(success);
            //			assert(Neighborhood.insert(std::make_pair(key, std::make_tuple(pN,pL,dir) )).second && "CANNOT INSERT IN NEIGHBORHOOD.");
            //
            //		}
            //
            //		/*****************************************************************************************/
            //		/* removeFromNeighborhood ****************************************************************/
            //		void removeFromNeighborhood(LinkType* const pL){
            //
            //			Derived* pN=NULL;
            //			size_t key=0;
            //			
            //			if (pL->source==this->p_derived()){
            //				pN=pL->sink;
            //				key=pN->sID;
            //				OutNeighborhood.erase(key);
            //			}
            //			
            //			if (pL->sink==this->p_derived()){
            //				pN=pL->source;
            //				key=pN->sID;
            //				InNeighborhood.erase(key);
            //			}
            //			
            //			bool success=Neighborhood.erase(key);
            //			assert(success);
            //			
            //		}

