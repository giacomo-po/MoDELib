/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_NETWORKLINK_H_
#define model_NETWORKLINK_H_

#ifndef VERBOSELEVEL
#define VERBOSELEVEL 0
#endif

#include <limits>
#include <iomanip>
#include <utility>
#include <algorithm>

#include <boost/ptr_container/ptr_map.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/utility.hpp>
#include <boost/shared_ptr.hpp>

#include <model/Network/Operations/includeNetworkOperations.h>

#include "model/Utilities/AddressBook.h"
#include "model/Network/SubNetwork.h"
#include <model/Utilities/TypeTraits.h>

namespace model {
	
	/**************************************************************************************/
	/**************************************************************************************/
	template <typename  Derived>
	class NetworkLink : boost::noncopyable,
	/*               */ public AddressBook<Derived>{
		
#include <model/Network/NetworkTypedefs.h>
		
		
		boost::shared_ptr<SubNetworkType> psn;

		
	public:
		
		// constant pointers to source and sink nodes
		NodeType* const source;		
		NodeType* const sink;		
		const std::pair<size_t,size_t> nodeIDPair;
		const FlowType flow;
		
	private:
//#include "model/Network/SubNetworkComponent.h"
		
		
		/******************************************************/
		void topologyChangeActions(){
		
			//! 1- Adds this to the source node neighborood and to the sink node neighborood
			source->addToNeighborhood(this->p_derived());
			sink  ->addToNeighborhood(this->p_derived());
			
			
			//! 2 - Joins source and sink SubNetworks			
			if (source->pSN()==sink->pSN()){ // source and sink are already in the same subnetwork
				psn=source->pSN();				// redirect psn to the source psn
				psn->add(this->p_derived());	// add this to the existing subnetwork
			}
			else{ // source and sink are in different subnetworks
				
				typedef void (NodeType::*node_member_function_pointer_type)(const boost::shared_ptr<SubNetworkType>&); 
				node_member_function_pointer_type Nmfp;
				Nmfp=&NodeType::formSubNetwork;
				typedef void (Derived::*link_member_function_pointer_type)(const boost::shared_ptr<SubNetworkType>&); 
				link_member_function_pointer_type Lmfp;
				Lmfp=&Derived::formSubNetwork;
				
				// find the size of the source and sink 
				size_t sourceSNsize(source->pSN()->nodeOrder());
				size_t   sinkSNsize(  sink->pSN()->nodeOrder());
				if (sourceSNsize>=sinkSNsize){
					psn=source->pSN();					   // redirect psn to the source psn
					psn->add(this->p_derived());		   // add this to the source subnetwork
					sink->depthFirstExecute(Nmfp,Lmfp,this->psn);   // Transmits 'formSubNetwork' on the sink side
				}
				else{
					psn=sink->pSN();				       // redirect psn to the sink psn
					psn->add(this->p_derived());	       // add this to the source subnetwork
					source->depthFirstExecute(Nmfp,Lmfp,this->psn); // Transmits 'formSubNetwork' on the source side
				}
			}

			//! 4- Calls 'topologyChangeActions' in source and sink
			source->topologyChangeActions();
			sink  ->topologyChangeActions();
			
		}
		
	public:
		
		/*****************************************************************************************/
		/* Constructor with Flow *****************************************************************/
		NetworkLink(const std::pair<NodeType* const,NodeType* const> & NodePair_in,
		/*       */ const FlowType & Flow_in) :
		/* initialization list               */ source(NodePair_in.first), 
		/* initialization list               */ sink(NodePair_in.second),
		/* initialization list               */ nodeIDPair(std::make_pair(source->sID,sink->sID)),
		/* initialization list               */ flow(Flow_in){ 
			topologyChangeActions();
		}
		
		/*****************************************************************************************/
		/* Constructor from ExpandingEdge ********************************************************/
		NetworkLink(const std::pair<NodeType* const,NodeType* const> & NodePair_in,
		/*       */ const ExpandingEdge<LinkType>& ee) :
		/* initialization list               */ source(NodePair_in.first), 
		/* initialization list               */ sink(NodePair_in.second),
		/* initialization list               */ nodeIDPair(std::make_pair(source->sID,sink->sID)),
		/* initialization list               */ flow(ee.E.flow){			
			topologyChangeActions();
		}
		
		/*****************************************************************************************/
		/* Destructor ****************************************************************************/
		~NetworkLink(){

			//! 1- Remove this from the Neighborood of Source and Sink nodes
			source->removeFromNeighborhood(this->p_derived());			
			sink  ->removeFromNeighborhood(this->p_derived());
			
			//! 2- Remove this from the SubNetwork	
			this->psn->remove(this->p_derived());

			//! 3- If Now Source and Sink are disconnected then reset the SubNetwork in the sink
			const bool sourceCanReachSink(source->depthFirstSearch(sink->sID));

			if (!sourceCanReachSink){
				sink -> resetPSN();
			}
			
			//! 4- Calls 'topologyChangeActions' in source and sink
			source->topologyChangeActions();
			sink  ->topologyChangeActions();
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
		/* isIncident ****************************************************************************/
		bool isIncident(const size_t & k) const{ // MARKED FOR REMOVE. CHANGE Network::remove first.
			return (source->sID==k || sink->sID==k);
		}
		
		
		////////////////////////////////////////////////////////
		// operator <<
		template <class T, typename OtherSegment>
		friend T& operator << (T& os, const NetworkLink<OtherSegment>& NL){
			
			os<<"			Link "<<std::setw(3)<<NL.sID<<" ("<<std::setw(3)<<NL.dID()<<") ["<<NL.p_derived()<<"]"<<
			" "<< std::setw(3)<<NL.source->sID<<
			" -> "<< std::setw(3)<<NL.sink->sID<<
			" flow = "<<NL.flow;
			
			return os;
		}
		
	};
	
	//////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////
} // namespace model
#endif