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

#include <tuple> // std::tuple replaces boost::tuple in c++11
#include <memory> // std::shared_ptr

#include <model/Network/NetworkComponent.h>
#include <model/Utilities/StaticID.h>
#include <model/Utilities/CRTP.h>
#include <model/Utilities/NonCopyable.h>

//#include "model/Network/NetworkLink.h"
#include <model/Network/Operations/includeNetworkOperations.h>


namespace model {
    
    template <typename Derived>
    class NetworkLink; // class predeclaration
    
    template <typename Derived>
    class NetworkNode : public NonCopyable,
    /*               */ public CRTP<Derived>,
    /*               */ public StaticID<Derived>
    {
        
    public:
#include <model/Network/NetworkTypedefs.h>
        friend class NetworkLink<LinkType>; // allow NetworkLink to call private NetworkNode::formNetworkComponent
        
        
    private:
        
        std::shared_ptr<NetworkComponentType> psn;
        
        
   	protected:
        
        /**********************************************************************/
        void resetPSN()
        {
            //! 1- Removes this from the current NetworkComponent
            this->psn->remove(this->p_derived());
            //! 2- Creates a new NetworkComponent containing this
            this->psn.reset(new NetworkComponentType(this->p_derived()));
            //! 3- Transmits 'formNetworkComponent' to the neighbors
            typedef void (Derived::*node_member_function_pointer_type)(const std::shared_ptr<NetworkComponentType>&);
            node_member_function_pointer_type Nmfp(&Derived::formNetworkComponent);
            //			Nmfp=&Derived::formNetworkComponent;
            typedef void (LinkType::*link_member_function_pointer_type)(const std::shared_ptr<NetworkComponentType>&);
            link_member_function_pointer_type Lmfp(&LinkType::formNetworkComponent);
            //			Lmfp=&LinkType::formNetworkComponent;
            depthFirstExecute(Nmfp,Lmfp,this->psn);
        }
        
        /**********************************************************************/
        void formNetworkComponent(const std::shared_ptr<NetworkComponentType> & psnOther)
        {/*!@param[in] psnOther a shared_ptr to another NetworkComponent
          *
          * If the shared_ptr to the NetworkComponent of *this is different from
          * psnOther, the former is reset (this may destroy the NetworkComponent)
          * and reassigned to psnOther.
          */
            if (psn!=psnOther)
            {
                psn->remove(this->p_derived());
                psn=psnOther;		// redirect psn to the new NetworkComponent
                psn->add(this->p_derived());	// add this in the new NetworkComponent
            }
        }
        
        /**********************************************************************/
        void addToNeighborhood(LinkType* const pL)
        {/*!@param[in] pL a pointer to a LinkType edge
          */
            
            if (pL->source==this->p_derived())
            {// this vertex is the source of edge *pL
                const NeighborType temp(pL->sink,pL,1);
                assert(OutNeighborhood.insert( std::make_pair(pL->sink->sID,temp) ).second && "CANNOT INSERT IN OUT_NEIGHBORHOOD");
                assert(   Neighborhood.insert( std::make_pair(pL->sink->sID,temp) ).second && "CANNOT INSERT IN NEIGHBORHOOD.");
            }
            else if (pL->sink==this->p_derived())
            {// this vertex is the sink of edge *pL
                const NeighborType temp(pL->source,pL,-1);
                assert(InNeighborhood.insert( std::make_pair(pL->source->sID,temp) ).second && "CANNOT INSERT IN IN_NEIGHBORHOOD");
                assert(  Neighborhood.insert( std::make_pair(pL->source->sID,temp) ).second && "CANNOT INSERT IN NEIGHBORHOOD.");
            }
            else
            {
                assert(0 && "CANNOT INSERT NON-INCIDENT EDGE");
            }
            
        }
        
        /**********************************************************************/
        void removeFromNeighborhood(LinkType* const pL)
        {
            
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
        
        /**********************************************************************/
        NetworkNode() :
        /* init list */ psn(new NetworkComponentType(this->p_derived()))
        {/*! Costructor with node arguments
          */
            Neighborhood.insert(std::make_pair(this->sID, NeighborType(this->p_derived(),(LinkType*) NULL,0) ));
        }
        
        /**********************************************************************/
        NetworkNode(const EdgeRef<LinkType>& ee) :
        /* init list */ psn(ee.E.pSN())
        {/*! Costructor from EdgeExpansion
          */
            // Insert this->p_derived() in the Neighborhood
            Neighborhood.insert(std::make_pair(this->sID, NeighborType(this->p_derived(),(LinkType*) NULL,0) ));
            
            // Manage NetworkComponent
            psn->add(this->p_derived());
        }
        
        /**********************************************************************/
        ~NetworkNode()
        {/*! The NetworkNode destructor performs two actions:
          */
            //! -1 Removes this from Neighborhood
            Neighborhood.erase(this->sID);
            
            //! -2 Removes this from the NetworkComponent
            this->psn->remove(this->p_derived());
        }
        
        /**********************************************************************/
        size_t snID() const
        {/*!\returns The NetworkComponent::snID() of the component
          * containing this.
          */
            return psn->snID(this->p_derived());
        }
        
        /**********************************************************************/
        const std::shared_ptr<NetworkComponentType> & pSN() const
        {/*!\returns A const reference to the shared-pointer to the
          * NetworkComponent containing this.
          */
            return psn;
        }
        
        /**********************************************************************/
        FlowType outFlow() const
        {
            //            FlowType Fout;
            FlowType Fout(FlowType::Zero()); // generalize
            Fout*=0.0;
            for (typename NeighborContainerType::const_iterator     NeighborIter=OutNeighborhood.begin();NeighborIter!=OutNeighborhood.end();++NeighborIter){
                Fout+=std::get<1>(NeighborIter->second)->flow;
            }
            return Fout;
        }
        
        /**********************************************************************/
        FlowType inFlow() const
        {
            //            FlowType Fin;
            FlowType Fin(FlowType::Zero());
            Fin*=0.0;
            for (typename NeighborContainerType::const_iterator     NeighborIter=InNeighborhood.begin();NeighborIter!=InNeighborhood.end();++NeighborIter){
                Fin+=std::get<1>(NeighborIter->second)->flow;
            }
            return Fin;
        }
        
        /**********************************************************************/
        bool depthFirstSearch (const size_t& ID, const size_t& N = ULONG_MAX) const
        {
            std::set<size_t> searchedNodes;
            return depthFirstSearch(searchedNodes,ID,N);
        }
        
        /**********************************************************************/
        bool depthFirstSearch (std::set<size_t>& searchedNodes, const size_t& ID, const size_t& N = ULONG_MAX) const
        {
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
        
        /**********************************************************************/
        template <typename T>
        void depthFirstExecute(void (Derived::*Nfptr)(const T&),void (LinkType::*Lfptr)(const T&), const T & input, const size_t& N = ULONG_MAX)
        {
            std::set<size_t> searchedNodes;
            std::set<std::pair<size_t,size_t> > searchedLinks;
            depthFirstExecute(searchedNodes,searchedLinks,Nfptr,Lfptr,input, N);
        }
        
        /**********************************************************************/
        template <typename T>
        void depthFirstExecute(std::set<size_t>& searchedNodes,
                               std::set<std::pair<size_t,size_t> >& searchedLinks,
                               void (Derived::*Nfptr)(const T&),void (LinkType::*Lfptr)(const T&),
                               const T & input,
                               const size_t& N = ULONG_MAX)
        {
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
        
        /**********************************************************************/
        template <typename T>
        void depthFirstNodeExecute(void (Derived::*Nfptr)(const T&),
                                   const T & input,
                                   const size_t& N = ULONG_MAX)
        {
            std::set<size_t> searchedNodes;
            depthFirstNodeExecute(searchedNodes,Nfptr,input, N);
        }
        
        /**********************************************************************/
        template <typename T>
        void depthFirstNodeExecute(std::set<size_t>& searchedNodes,
                                   void (Derived::*Nfptr)(const T&),
                                   const T & input,
                                   const size_t& N = ULONG_MAX)
        {
            (this->p_derived()->*Nfptr)(input); // execute Nfptr on this node
            if (N!=0)
            {
                assert(searchedNodes.insert(this->sID).second && "CANNOT INSERT CURRENT NODE IN SEARCHED NODES"); // this node has been searched
                for (typename NeighborContainerType::iterator NeighborIter=Neighborhood.begin();NeighborIter!=Neighborhood.end();++NeighborIter)
                {
                    if (searchedNodes.find(std::get<0>(NeighborIter->second)->sID)==searchedNodes.end())
                    {  // neighbor not searched
                        std::get<0>(NeighborIter->second)->depthFirstNodeExecute(searchedNodes,Nfptr,input, N-1); // continue executing on neighbor
                    }
                }
            }
        }
        
        /**********************************************************************/
        void depthFirstNodeExecute(void (Derived::*Nfptr)(void),
                                   const size_t& N = ULONG_MAX)
        {
            std::set<size_t> searchedNodes;
            depthFirstNodeExecute(searchedNodes,Nfptr, N);
        }
        
        /**********************************************************************/
        void depthFirstNodeExecute(std::set<size_t>& searchedNodes,
                                   void (Derived::*Nfptr)(void),
                                   const size_t& N = ULONG_MAX)
        {
            (this->p_derived()->*Nfptr)(); // execute Nfptr on this node
            if (N!=0)
            {
                assert(searchedNodes.insert(this->sID).second && "CANNOT INSERT CURRENT NODE IN SEARCHED NODES"); // this node has been searched
                for (typename NeighborContainerType::iterator NeighborIter=Neighborhood.begin();NeighborIter!=Neighborhood.end();++NeighborIter)
                {
                    if (searchedNodes.find(std::get<0>(NeighborIter->second)->sID)==searchedNodes.end())
                    {  // neighbor not searched
                        std::get<0>(NeighborIter->second)->depthFirstNodeExecute(searchedNodes,Nfptr, N-1); // continue executing on neighbor
                    }
                }
            }
        }
        
        /**********************************************************************/
        const NeighborContainerType& neighborhood() const
        {
            return Neighborhood;
        }
        
        /**********************************************************************/
        NodeType* closedNeighborNode(const unsigned int& k) const
        {
            assert(k<Neighborhood.size() && "INDEX EXCEEDS SIZE");
            typename NeighborContainerType::const_iterator nIter(Neighborhood.begin());
            for (unsigned int n=0;n<k;++n)
            {
                nIter++;
            }
            return (std::get<0>(nIter->second));
        }
        
        /**********************************************************************/
        NodeType* openNeighborNode(const unsigned int& k) const
        {
            assert(k<(Neighborhood.size()-1) && "INDEX EXCEEDS SIZE");
            typename NeighborContainerType::const_iterator nIter(Neighborhood.begin());
            if(nIter->first == this->sID)
            {
                nIter++;
            }
            for (unsigned int n=0;n<k;++n)
            {
                if(nIter->first == this->sID)
                {
                    nIter++;
                }
                nIter++;
            }
            if(nIter->first == this->sID)
            {
                nIter++;
            }
            return (std::get<0>(nIter->second));
        }
        
        /**********************************************************************/
        LinkType* closedNeighborLink(const unsigned int& k) const
        {
            assert(k<Neighborhood.size() && "INDEX EXCEEDS SIZE");
            typename NeighborContainerType::const_iterator nIter(Neighborhood.begin());
            for (unsigned int n=0;n<k;++n)
            {
                nIter++;
            }
            return (std::get<1>(nIter->second));
        }
        
        /**********************************************************************/
        LinkType* openNeighborLink(const unsigned int& k) const
        {
            assert(k<(Neighborhood.size()-1) && "INDEX EXCEEDS SIZE");
            typename NeighborContainerType::const_iterator nIter(Neighborhood.begin());
            if(nIter->first == this->sID)
            {
                nIter++;
            }
            for (unsigned int n=0;n<k;++n)
            {
                if(nIter->first == this->sID)
                {
                    nIter++;
                }
                nIter++;
            }
            if(nIter->first == this->sID)
            {
                nIter++;
            }
            return (std::get<1>(nIter->second));
        }
        
        /**********************************************************************/
        const NeighborContainerType & outNeighborhood() const
        {
            return OutNeighborhood;
        }
        
        /**********************************************************************/
        const NeighborContainerType & inNeighborhood() const
        {
            return InNeighborhood;
        }
        
        /**********************************************************************/
        size_t neighborID(const size_t & k) const
        {
            return std::distance(Neighborhood.begin(),Neighborhood.find(k));
        }
        
        /**********************************************************************/
        size_t outOrder() const
        {
            return OutNeighborhood.size();
        }
        
        /**********************************************************************/
        size_t inOrder() const
        {
            return InNeighborhood.size();
        }
        
        /**********************************************************************/
        size_t openOrder() const
        {
            return outOrder()+inOrder();
        }
        
        /**********************************************************************/
        size_t closedOrder() const
        {
            return openOrder()+1;
        }
        
        /**********************************************************************/
        bool is_source() const // THIS IS MEANINGLESS IF LINK DIRECTIONS CAN BE SWITCHED ARBITRARILY
        {
            return inOrder()==0 && outOrder()>0;
        }
        
        /**********************************************************************/
        bool is_sink() const // THIS IS MEANINGLESS IF LINK DIRECTIONS CAN BE SWITCHED ARBITRARILY
        {
            return outOrder()==0 && inOrder()>0;
        }
        
        /**********************************************************************/
        bool is_isolated() const
        {
            return inOrder()==0 && outOrder()==0;
        }
        
        /**********************************************************************/
        bool is_central() const
        {
            return inOrder()>0 && outOrder()>0;
        }
        
        /**********************************************************************/
        bool is_balanced() const
        {
            return (outFlow() - inFlow()).squaredNorm()==0;
        }
        
        /**********************************************************************/
        bool is_simple() const
        {
            return (outOrder()+inOrder())==2 && is_balanced();
        }
        
        /**********************************************************************/
        template <class T, typename OtherDerived>
        friend T& operator << (T& os, const NetworkNode<OtherDerived> & NN)
        {
            
            os << "Node sID=" << std::setw(3)<<NN.sID
            << " dID=" << std::setw(3)<<NN.dID()
            << " snID=" << std::setw(3)<<NN.snID()<<" ("
            << std::setw(3)<<NN.dID() <<") ["<<NN.p_derived()<<"] : ["<<NN.address(NN.sID)<<"]"<<
            " state = "<<NN.get_state()<<std::endl;
            
            os << "		OutLinks: (outOrder= "<<NN.outOrder()<<")"<<std::endl;
            for(typename NeighborContainerType::const_iterator	 NeighborIter=NN.OutNeighborhood.begin();NeighborIter!=NN.OutNeighborhood.end();++NeighborIter)
            {
                os<< *std::get<1>(NeighborIter->second)<<std::endl;//
            }
            
            os << "		InLinks: (inOrder= "<<NN.inOrder()<<")"<<std::endl;
            for(typename NeighborContainerType::const_iterator	 NeighborIter=NN.InNeighborhood.begin();NeighborIter!=NN.InNeighborhood.end();++NeighborIter)
            {
                os<< *std::get<1>(NeighborIter->second)<<std::endl;//
            }
            
            return os;
        }
        
    };
    

} // namespace model
#endif
