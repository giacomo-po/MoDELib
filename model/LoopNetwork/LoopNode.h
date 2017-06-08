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
    /*            */ private std::set<LoopLink<typename TypeTraits<Derived>::LinkType>*>
    {
        
        
        typedef typename TypeTraits<Derived>::LinkType LinkType;
        typedef LoopLink<LinkType> LoopLinkType;
        typedef std::set<LoopLinkType*> LoopLinkContainerType;
        
    public:
        
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
        void addLoopLink(LoopLinkType* const pL)
        {
//            std::cout<<"LoopNode "<<this->sID<<" adding Link "<<pL->source->sID<<"->"<<pL->sink->sID<<std::endl;
            // Store the connecting link
            const bool inserted=this->insert(pL).second;
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
            const int erased=this->erase(pL);
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
