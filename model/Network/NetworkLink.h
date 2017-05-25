/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_NETWORKLINK_H_
#define model_NETWORKLINK_H_


#include <limits>
#include <iomanip>
#include <utility>
#include <algorithm>
#include <memory> // std::shared_ptr

//#include <boost/ptr_container/ptr_map.hpp>

#include <model/Network/NetworkComponent.h>
#include <model/Network/Operations/includeNetworkOperations.h>
#include <model/Utilities/CRTP.h>
#include <model/Utilities/StaticID.h>
#include <model/Network/NetworkNode.h>
#include <model/Utilities/TypeTraits.h>
#include <model/Utilities/NonCopyable.h>

namespace model
{
	
	/**************************************************************************/
	/**************************************************************************/
	template <typename  Derived>
	class NetworkLink : public NonCopyable,
    /*               */ public CRTP<Derived>,
	/*               */ public StaticID<Derived>
    {
#include <model/Network/NetworkTypedefs.h>
		
        friend class NetworkNode<NodeType>; // allow NetworkNode to call private NetworkLink::formNetworkComponent
		
		std::shared_ptr<NetworkComponentType> psn;

		
	public:
		
		// constant pointers to source and sink nodes
		NodeType* const source;		
		NodeType* const sink;		
		const std::pair<size_t,size_t> nodeIDPair;
		const FlowType flow;
		
	private:
		
        /**********************************************************************/
        void changeSN(const NetworkComponent<NodeType,LinkType>& SN)
        {
        
            const std::map<size_t,NodeType* const> tempNodeMap(SN);
            const std::map<std::pair<size_t,size_t>,LinkType* const> tempLinkMap(SN);

            for (typename std::map<size_t,NodeType* const>::const_iterator vIter=tempNodeMap.begin();vIter!=tempNodeMap.end();++vIter)
            {
                vIter->second->formNetworkComponent(psn);
            }
            
            for (typename std::map<std::pair<size_t,size_t>,LinkType* const>::const_iterator lIter=tempLinkMap.begin();lIter!=tempLinkMap.end();++lIter)
            {
                lIter->second->formNetworkComponent(psn);
            }
            
        }
        
		/**********************************************************************/
		void makeTopologyChange()
        {
			//! 1- Adds this to the source node neighborood and to the sink node neighborood
			source->addToNeighborhood(this->p_derived());
			sink  ->addToNeighborhood(this->p_derived());
			
			
			//! 2 - Joins source and sink NetworkComponents			
			if (source->pSN()==sink->pSN()) // source and sink are already in the same NetworkComponent
            { 
				psn=source->pSN();				// redirect psn to the source psn
				psn->add(this->p_derived());	// add this to the existing NetworkComponent
			}
			else // source and sink are in different NetworkComponents
            {
				// find the size of the source and sink 
				size_t sourceSNsize(source->pSN()->nodeOrder());
				size_t   sinkSNsize(  sink->pSN()->nodeOrder());
				if (sourceSNsize>=sinkSNsize)
                {
					psn=source->pSN();					   // redirect psn to the source psn
					psn->add(this->p_derived());		   // add this to the source NetworkComponent
                    changeSN(*(sink->psn.get()));
//					sink->depthFirstExecute(Nmfp,Lmfp,this->psn);   // Transmits 'formNetworkComponent' on the sink side
				}
				else
                {
					psn=sink->pSN();				       // redirect psn to the sink psn
					psn->add(this->p_derived());	       // add this to the source NetworkComponent
                    changeSN(*(source->psn.get()));
//                    source->depthFirstExecute(Nmfp,Lmfp,this->psn); // Transmits 'formNetworkComponent' on the source side
				}
			}
			
		}
        
        /* formNetworkComponent ***********************************************/
		void formNetworkComponent(const std::shared_ptr<NetworkComponentType> & psnOther)
        {
			if (psn!=psnOther)
            {
				psn->remove(this->p_derived());
				psn=psnOther;		// redirect psn to the new NetworkComponent
				psn->add(this->p_derived());    // add this in the new NetworkComponent
			}
		}
		
	public:
		
		/* Constructor with Flow **********************************************/
		NetworkLink(const std::pair<NodeType* const,NodeType* const> & NodePair_in,
		/*       */ const FlowType & Flow_in) :
		/* init list */ source(NodePair_in.first), 
		/* init list */ sink(NodePair_in.second),
		/* init list */ nodeIDPair(std::make_pair(source->sID,sink->sID)),
		/* init list */ flow(Flow_in)
        {/*! Constructor with pointers to source and sink, and flow
          *  @param[in] NodePair_in the pair of source and sink pointers
          *  @param[in] Flow_in the input flow
          */
//            std::cout<<"Creating NetworkLink "<<source->sID<<"->"<<sink->sID<<std::endl;
			makeTopologyChange();
        }
		
		/* Constructor from EdgeRef *************************************/
		NetworkLink(const std::pair<NodeType* const,NodeType* const> & NodePair_in,
		/*       */ const EdgeRef<LinkType>& ee) :
		/* init list */ source(NodePair_in.first), 
		/* init list */ sink(NodePair_in.second),
		/* init list */ nodeIDPair(std::make_pair(source->sID,sink->sID)),
		/* init list */ flow(ee.E.flow)
        {
			
//            if (source==ee.E.source) // first new link in expansion
//            {
//                assert(sink!=ee.E.sink && "FIRST LINK IN EXPANSION FINDS WRONG SINK");
//            }
//            else // second new link in expansion
//            {
//                assert(sink==ee.E.sink && "SECOND LINK IN EXPANSION FINDS WRONG SINK");
//                assert(source!=ee.E.source && "SECOND LINK IN EXPANSION FINDS WRONG SOURCE");
//            }
//            
//			makeTopologyChange(ee);
//                        std::cout<<"Creating (2) NetworkLink "<<source->sID<<"->"<<sink->sID<<std::endl;
            makeTopologyChange();
		}
		
		/* Destructor *********************************************************/
		~NetworkLink()
        {
            
                        //std::cout<<"Destroying NetworkLink "<<source->sID<<"->"<<sink->sID<<std::endl;
            
			//! 1- Remove this from the Neighborood of Source and Sink nodes
			source->removeFromNeighborhood(this->p_derived());			
			sink  ->removeFromNeighborhood(this->p_derived());
			
			//! 2- Remove this from the NetworkComponent	
			this->psn->remove(this->p_derived());

			//! 3- If Now Source and Sink are disconnected then reset the NetworkComponent in the sink
			const bool sourceCanReachSink(source->depthFirstSearch(sink->sID));

			if (!sourceCanReachSink)
            {
				sink -> resetPSN();
			}
        }
		
		/* snID ***************************************************************/
		size_t snID() const
        {/*! @return The StaticID of the NetworkComponent containing this
          */
			return psn->snID(this->p_derived());
		}
		
		/* pSN ****************************************************************/
		const std::shared_ptr<NetworkComponentType> & pSN() const
        {/*! @return The shared_ptr to the NetworkComponent containing this
          */
			return psn;
		}
		
		/**********************************************************************/
		bool isIncident(const size_t & k) const
        { // MARKED FOR REMOVE. CHANGE Network::remove first.
			return (source->sID==k || sink->sID==k);
		}
		
		/**********************************************************************/
		template <class T, typename OtherSegment>
		friend T& operator << (T& os, const NetworkLink<OtherSegment>& NL)
        {
			os<<NL.source->sID<<" "<<NL.sink->sID<<" "<<NL.source->sID<<NL.flow;
			return os;
		}
		
	};
	
} // namespace model
#endif
