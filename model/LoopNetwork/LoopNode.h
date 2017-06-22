/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LoopNode_H_
#define model_LoopNode_H_

#include <iostream>
//#include <list>
#include <set>
#include <memory>
#include <assert.h>
#include <tuple>

//#include <iterator>

#include <model/Utilities/StaticID.h>
#include <model/Utilities/CRTP.h>
//#include <iomanip>
//#include <vector>
//#include <Eigen/Dense>
//#include <Eigen/StdVector>
#include <model/MPI/MPIcout.h>
#include <model/LoopNetwork/NodeObserver.h>
#include <model/LoopNetwork/LoopLink.h>

#define VerboseLoopNode(N,x) if(verboseLevel>=N){model::cout<<x;}


namespace model
{
    template<typename Derived>
    class LoopNode : public StaticID<Derived>,
    /*            */ public CRTP<Derived>,
    /*            */ private std::set<LoopLink<typename TypeTraits<Derived>::LinkType>*>,
    /*            */ private std::map<size_t,std::tuple<Derived* const ,typename TypeTraits<Derived>::LinkType* const,short int>>
    {
        
    public:

        typedef typename TypeTraits<Derived>::LinkType LinkType;
        typedef LoopLink<LinkType> LoopLinkType;
        typedef std::set<LoopLinkType*> LoopLinkContainerType;
        typedef std::map<size_t,LoopLinkContainerType> LinkByLoopContainerType;
        typedef std::tuple<Derived* const ,LinkType* const,short int>				NeighborType;
        typedef std::map<size_t,NeighborType>						    	NeighborContainerType;

        
        static int verboseLevel;
        
        /**********************************************************************/

        LoopNode(const LoopNode&) =delete;
    
        /**********************************************************************/
        LoopNode()
        {
            VerboseLoopNode(1,"Contructing LoopNode "<<name()<<std::endl);

//            std::cout<<"Constructing LoopNode "<<this->sID<<std::endl;
            NodeObserver<Derived>::addNode(this->p_derived());
        }
        
        /**********************************************************************/
        ~LoopNode()
        {
            VerboseLoopNode(1,"Destroying LoopNode "<<name()<<std::endl);

            NodeObserver<Derived>::removeNode(this->p_derived());
            
            assert(loopLinks().empty());
        }
    
        /**********************************************************************/
        const LoopLinkContainerType& loopLinks() const
        {
            return *this;
        }
        
        /**********************************************************************/
        LoopLinkContainerType& loopLinks()
        {
            return *this;
        }
        
        
        /**********************************************************************/
        LinkByLoopContainerType linksByLoopID() const
        {
            LinkByLoopContainerType temp;
            for(const auto& link : loopLinks())
            {
                temp[link->loop()->sID].insert(link);
            }
            return temp;
        }
        
        /**********************************************************************/
        const NeighborContainerType& neighbors() const
        {
            return *this;
        }
        
        /**********************************************************************/
        NeighborContainerType& neighbors()
        {
            return *this;
        }
        
        /**********************************************************************/
        void addToNeighborhood(LinkType* const pL)
        {/*!@param[in] pL a pointer to a LinkType edge
          */
            
            if (pL->source->sID==this->sID)
            {// this vertex is the source of edge *pL
                const NeighborType temp(pL->sink.get(),pL,1);
                const bool success=neighbors().insert( std::make_pair(pL->sink->sID,temp) ).second;
                assert(success && "CANNOT INSERT IN NEIGHBORHOOD.");
            }
            else if (pL->sink->sID==this->sID)
            {// this vertex is the sink of edge *pL
                const NeighborType temp(pL->source.get(),pL,-1);
                const bool success=neighbors().insert( std::make_pair(pL->source->sID,temp) ).second;
                assert(success  && "CANNOT INSERT IN NEIGHBORHOOD.");
            }
            else
            {
                assert(0 && "CANNOT INSERT NON-INCIDENT EDGE");
            }
            
        }
        
        /**********************************************************************/
        void removeFromNeighborhood(LinkType* const pL)
        {
            if (pL->source->sID==this->sID)
            {
                const size_t key=pL->sink->sID;
                int success=neighbors().erase(key);
                assert(success==1);
            }
            else if (pL->sink->sID==this->sID)
            {
                const size_t key=pL->source->sID;
                int success=neighbors().erase(key);
                assert(success==1);
            }
            else
            {
                assert(0 && "CANNOT REMOVE FROM NEIGHBORS");
            }
        }
        
        /**********************************************************************/
        void addLoopLink(LoopLinkType* const pL)
        {
//            std::cout<<"LoopNode "<<this->sID<<" adding Link "<<pL->source->sID<<"->"<<pL->sink->sID<<std::endl;
            // Store the connecting link
            const bool inserted=loopLinks().insert(pL).second;
            assert(inserted);

//            NodeObserver<Derived>::addNode(pL->source());
//            NodeObserver<Derived>::addNode(pL->sink());

            
            // form the prev/next structure
            for(const auto& other : loopLinks())
            {


                if(pL->loop().get()==other->loop().get() && pL!=other) //two links in same loop
                {
//                std::cout<<"other is "<<other->source->sID<<"->"<<other->sink->sID<<std::endl;
                    if(pL->source().get()==other->sink().get())
                    {
                        assert(pL->prev==nullptr || pL->prev==other);
                        assert(other->next==nullptr || other->next==pL);
                        pL->prev=other;
                        other->next=pL;
                    }
                    if (pL->sink().get()==other->source().get())
                    {
//                        if(pL->next!=nullptr)
//                        {
//                        std::cout<<"pL->next="<<pL->next->source->sID<<"->"<<pL->next->sink->sID<<std::endl;
//                        }
//                        if(other->prev!=nullptr)
//                        {
//                            std::cout<<"other->prev="<<other->prev->source->sID<<"->"<<other->prev->sink->sID<<std::endl;
//                        }
                        assert(pL->next==nullptr || pL->next==other);
                        assert(other->prev==nullptr || other->prev==pL);
                        pL->next=other;
                        other->prev=pL;
                    }
                    if(pL->source().get()==other->source().get() || pL->sink().get()==other->sink().get())
                    {
                        assert(0 && "Not a Loop");
                    }
                }
//                else
//                {
//                    assert(0 && " two links in same loop ");
//                }

            }

            
        }
        
        /**********************************************************************/
        void removeLoopLink(LoopLinkType* const pL)
        {
            const int erased=loopLinks().erase(pL);
            assert(erased==1);
            
            if(pL->next!=nullptr)
            {
                pL->next->prev=nullptr;
                pL->next=nullptr;
            }

            if(pL->prev!=nullptr)
            {
                pL->prev->next=nullptr;
                pL->prev=nullptr;
            }

        }
        
        /**********************************************************************/
        std::string name() const
        {/*!\returns the string "i" where i is this->sID
          */
            return std::to_string(this->sID) ;
        }
        
    };
    
    template<typename Derived>
    int LoopNode<Derived>::verboseLevel=1;

    
    
}
#endif
