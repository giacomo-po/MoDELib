/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PeriodicDislocationLoop_H_
#define model_PeriodicDislocationLoop_H_


#include <memory>
#include <string>
#include <list>



#include <StaticID.h>
//#include <PeriodicLoopObserver.h>
#include <PeriodicGlidePlane.h>
#include <PolygonClipper.h>

#ifndef NDEBUG
#define VerbosePeriodicDislocationBase(N,x) if(this->verbosePeriodicDislocationBase>=N){model::cout<<x;}
#else
#define VerbosePeriodicDislocationBase(N,x)
#endif

namespace model
{
    
    struct PeriodicDislocationBase
    {
        static int verbosePeriodicDislocationBase;
        static void initFromFile(const std::string& fileName)
        {
            verbosePeriodicDislocationBase=TextFileParser(fileName).readScalar<int>("verbosePeriodicDislocationBase",true);
        }

    };
    int PeriodicDislocationBase::verbosePeriodicDislocationBase=0;

    template<typename DislocationNetworkType>
    struct PeriodicDislocationNode;
    
    template<typename DislocationNetworkType>
    class PeriodicDislocationLoop;
    
    template<typename DislocationNetworkType>
    struct PeriodicLoopLink;
    
    template <typename DislocationNetworkType>
    struct PeriodicDislocationLoopFactory;
    
    template <typename DislocationNetworkType>
    struct TypeTraits<PeriodicDislocationLoopFactory<DislocationNetworkType>>
    {
        typedef PeriodicDislocationLoop<DislocationNetworkType> ValueType;
        typedef GlidePlaneKey<DislocationNetworkType::dim> KeyType;
        typedef std::less<KeyType> CompareType;
    };
    
    
    /**************************************************************************/
    template <typename DislocationNetworkType>
    struct PeriodicDislocationLoopFactory : public PeriodicDislocationBase
    /*                                   */,public KeyConstructableWeakPtrFactory<PeriodicDislocationLoopFactory<DislocationNetworkType>>
    {
        static constexpr int dim=DislocationNetworkType::dim;
//        typedef std::array<long int, dim + 3> KeyType;
        typedef PeriodicDislocationLoopFactory<DislocationNetworkType> PeriodicDislocationLoopFactoryType;
        typedef PeriodicDislocationLoop<DislocationNetworkType> PeriodicDislocationLoopType;
        typedef std::shared_ptr<PeriodicDislocationLoopType> PeriodicDislocationLoopSharedPtrType;
        typedef KeyConstructableWeakPtrFactory<PeriodicDislocationLoopFactoryType> BaseType;
        typedef GlidePlane<dim> GlidePlaneType;

        PeriodicGlidePlaneFactory<dim> periodicGlidePlaneFactory;
        
        PeriodicDislocationLoopFactory(const Polycrystal<dim>& poly_in,
                                       GlidePlaneFactory<dim>& glidePlaneFactory_in) :
        /* init */ periodicGlidePlaneFactory(poly_in,glidePlaneFactory_in)
        {
            std::cout<<greenBoldColor<<"Creating PeriodicDislocationLoopFactory"<<std::endl;
        }
        
        PeriodicDislocationLoopSharedPtrType get(const GlidePlaneType& plane)
        {
            VerbosePeriodicDislocationBase(2,"PeriodicDislocationLoopFactory::get"<<std::endl;);
            return BaseType::get(periodicGlidePlaneFactory.periodicPlaneKey(plane.key));
        }
        
        /**********************************************************************/
        BaseType& periodicDislocationLoops()
        {
            return *this;
        }
        
        /**********************************************************************/
        const BaseType& periodicDislocationLoops() const
        {
            return *this;
        }
    };
    
    /**************************************************************************/
    template <typename DislocationNetworkType>
    struct NeighborConnectivityPair : public std::pair<typename DislocationNetworkType::VectorDim,std::set<const PeriodicLoopLink<DislocationNetworkType>*>>
    {
        
        typedef std::pair<typename DislocationNetworkType::VectorDim,std::set<const PeriodicLoopLink<DislocationNetworkType>*>> BaseType;
        
        NeighborConnectivityPair() :
        /* init */ BaseType(DislocationNetworkType::VectorDim::Zero(),std::set<const PeriodicLoopLink<DislocationNetworkType>*>())
        {
            
        }
        
        NeighborConnectivityPair(const NeighborConnectivityPair&) = delete;
        const NeighborConnectivityPair& operator=(const NeighborConnectivityPair&) = delete;

    };
    
    
    template <typename DislocationNetworkType>
    class PeriodicNeighborConnectivity : public PeriodicDislocationBase
    /*                                */,public std::map<const PeriodicDislocationNode<DislocationNetworkType>*,NeighborConnectivityPair<DislocationNetworkType>>
    {
        
        
        typedef PeriodicDislocationNode<DislocationNetworkType> PeriodicDislocationNodeType;
        typedef PeriodicLoopLink<DislocationNetworkType> PeriodicLoopLinkType;
        typedef std::map<const PeriodicDislocationNode<DislocationNetworkType>*,NeighborConnectivityPair<DislocationNetworkType>> BaseType;
        typedef typename DislocationNetworkType::VectorDim VectorDim;
        
        BaseType& base()
        {
            return *this;
        }
        
        const PeriodicDislocationNodeType* const node;
        
    public:
        
        PeriodicNeighborConnectivity(const PeriodicDislocationNodeType* const node_in) :
        /* init */ node(node_in)
        {
            
        }
        
        void removeUntwinnedEdges(const std::set<const PeriodicLoopLinkType*>& links)
        {
            
            for(const auto& link : links)
            {
                node->periodicLoop.removeUntwinnedEdge(link);
            }
            
        }
        
        void addUntwinnedEdges(const std::set<const PeriodicLoopLinkType*>& links)
        {
            
            for(const auto& link : links)
            {
                node->periodicLoop.addUntwinnedEdge(link);
            }
            
        }
        
        void addPeriodicLoopLink(PeriodicLoopLinkType *const periodicLoopLink)
        {
            if (periodicLoopLink->source->sID == node->sID)
            {// link is an out-link
                const auto neighborNode(periodicLoopLink->sink.get());
                auto& neighborPair(base()[neighborNode]);
                assert(neighborPair.second.find(periodicLoopLink)==neighborPair.second.end() && "NEIGHBOR ALREADY PRESENT IN NEIGHBORS SET");
                neighborPair.first+=periodicLoopLink->loopLink->loop()->burgers();
                neighborPair.second.insert(periodicLoopLink);
                addUntwinnedEdges(neighborPair.second);
                if(neighborPair.first.squaredNorm()<FLT_EPSILON)
                {
                    removeUntwinnedEdges(neighborPair.second);
                }
            }
            else if (periodicLoopLink->sink->sID == node->sID)
            {// link is an in-link
                const auto neighborNode(periodicLoopLink->source.get());
                auto& neighborPair(base()[neighborNode]);
                assert(neighborPair.second.find(periodicLoopLink)==neighborPair.second.end() && "NEIGHBOR ALREADY PRESENT IN NEIGHBORS SET");
                neighborPair.first-=periodicLoopLink->loopLink->loop()->burgers();
                neighborPair.second.insert(periodicLoopLink);
                addUntwinnedEdges(neighborPair.second);
                if(neighborPair.first.squaredNorm()<FLT_EPSILON)
                {
                    removeUntwinnedEdges(neighborPair.second);
                }
            }
            else
            {
                assert(false && "node MUST BE EITHER SOURCE OR SINK");
            }
        }
        
        void removePeriodicLoopLink(PeriodicLoopLinkType *const periodicLoopLink)
        {
            if (periodicLoopLink->source->sID == node->sID)
            {// link is an out-link
                const auto neighborNode(periodicLoopLink->sink.get());
                auto neighborIter(base().find(neighborNode));
                assert(neighborIter!=base().end() && "Node not found in neighbors");
                assert(neighborIter->second.second.find(periodicLoopLink)!=neighborIter->second.second.end() && "link not found in neighbors");
                neighborIter->second.first-=periodicLoopLink->loopLink->loop()->burgers();
                addUntwinnedEdges(neighborIter->second.second);
                if(neighborIter->second.first.squaredNorm()<FLT_EPSILON)
                {
                    VerbosePeriodicDislocationBase(2,"PeriodicNeighborConnectivity "<<node->sID<<" case A"<<std::endl;);
                    removeUntwinnedEdges(neighborIter->second.second);
                }
                else
                {
                    VerbosePeriodicDislocationBase(2,"PeriodicNeighborConnectivity "<<node->sID<<" case B"<<std::endl;);
                    node->periodicLoop.removeUntwinnedEdge(periodicLoopLink);
                }
                neighborIter->second.second.erase(periodicLoopLink);
            }
            else if (periodicLoopLink->sink->sID == node->sID)
            {// link is an in-link
                const auto neighborNode(periodicLoopLink->source.get());
                auto neighborIter(base().find(neighborNode));
                assert(neighborIter!=base().end() && "Node not found in neighbors");
                assert(neighborIter->second.second.find(periodicLoopLink)!=neighborIter->second.second.end() && "link not found in neighbors");
                neighborIter->second.first+=periodicLoopLink->loopLink->loop()->burgers();
                addUntwinnedEdges(neighborIter->second.second);
                if(neighborIter->second.first.squaredNorm()<FLT_EPSILON)
                {
                    VerbosePeriodicDislocationBase(2,"PeriodicNeighborConnectivity "<<node->sID<<" case C"<<std::endl;);
                    removeUntwinnedEdges(neighborIter->second.second);
                }
                else
                {
                    VerbosePeriodicDislocationBase(2,"PeriodicNeighborConnectivity "<<node->sID<<" case D"<<std::endl;);
                    node->periodicLoop.removeUntwinnedEdge(periodicLoopLink);
                }
                neighborIter->second.second.erase(periodicLoopLink);
            }
            else
            {
                assert(false && "node MUST BE EITHER SOURCE OR SINK");
            }
        }
        
    };
    
    
    template <typename DislocationNetworkType>
    struct SuperNodalConnectivity
    {
        typedef PeriodicLoopLink<DislocationNetworkType> PeriodicLoopLinkType;
        PeriodicLoopLinkType *inEdge;
        PeriodicLoopLinkType *outEdge;
        
        SuperNodalConnectivity() :
        /* init */ inEdge(nullptr)
        /* init */,outEdge(nullptr)
        {
        }
    };
    
    template<typename DislocationNetworkType>
    struct PeriodicDislocationNode : public StaticID<PeriodicDislocationNode<DislocationNetworkType>>
//    /*                            */,public Eigen::Matrix<double, DislocationNetworkType::dim-1, 1>
    /*                            */,public PeriodicNeighborConnectivity<DislocationNetworkType>
    {
        
        typedef PeriodicDislocationNode<DislocationNetworkType> PeriodicDislocationNodeType;
        typedef typename DislocationNetworkType::NodeType NodeType;
        typedef typename DislocationNetworkType::LoopType LoopType;
        typedef typename LoopType::LoopLinkType LoopLinkType;
        typedef PeriodicLoopLink<DislocationNetworkType> PeriodicLoopLinkType;
        typedef typename DislocationNetworkType::VectorDim VectorDim;
        typedef Eigen::Matrix<double,DislocationNetworkType::dim - 1, 1> VectorLowerDim;
        typedef SuperNodalConnectivity<DislocationNetworkType> SuperNodalConnectivityType;
        typedef std::map<size_t, SuperNodalConnectivityType> SuperNodalConnectivityContainerType;
        typedef PeriodicDislocationLoop<DislocationNetworkType> PeriodicDislocationLoopType;
        typedef std::set<PeriodicLoopLinkType*> InOutEdgeContainerType;
        
        PeriodicDislocationLoopType& periodicLoop;
        
        typedef std::map<NodeType*,VectorDim> RVEnodesSetType; //Value represents the shifts
        
//        std::map<VectorDim,std::pair<const NodeType* const,VectorLowerDim>,CompareVectorsByComponent<double,DislocationNetworkType::dim,float>> nodesMap;
        RVEnodesSetType rveNodes;

    private:
        
        VectorLowerDim P;

        SuperNodalConnectivityContainerType _loopConnectivities;
        
    public:


        /**********************************************************************/
        PeriodicDislocationNode(PeriodicDislocationLoopType& periodicLoop_in) : // REMOVE THIS CONSTRUCTOR
//        /* init */ VectorLowerDim(pos)
        /* init */ PeriodicNeighborConnectivity<DislocationNetworkType>(this)
        /* init */,periodicLoop(periodicLoop_in)
        /* init */,P(VectorLowerDim::Zero())
        {
            VerbosePeriodicDislocationBase(2,"Creating PeriodicDislocationNode "<<this->sID<<std::endl;);
        }
        
        /**********************************************************************/
        ~PeriodicDislocationNode()
        {
            VerbosePeriodicDislocationBase(2,"Destroying PeriodicDislocationNode "<<this->sID<<std::endl;);
        }
        
        void set_P(const VectorLowerDim& newP)
        {
            P=newP;
        }
        
        
        void addNode( NodeType* node,const VectorDim& shift)
        {
            rveNodes.emplace(node,shift);
            P=periodicLoop.periodicGlidePlane->getLocalPosition(node->get_P(),shift);
        }
                
        void removeNode( NodeType* node)
        {
            const size_t erased(rveNodes.erase(node));
            if (erased!=1)
            {
                std::cout<<" Erasing RVE nodes for "<<this->sID<<" periodic Node..... Printing RVE nodes"<<std::endl;
                for (const auto& nodetemp : rveNodes)
                {
                    std::cout<<nodetemp.first->sID<<" ---> "<<nodetemp.second.transpose()<<std::endl;
                }
            }
            assert(erased==1 && "UNABLE TO REMOVE DislocationNode from nodesMap");
        }
        
        void updateShifts(NodeType* nodeIn, const VectorDim& oldPosition, const VectorDim& newPosition ,const VectorDim& Loopshift)
        {
            //loop over outer boundary of the patch and determine the shift (there should be only one outerboundary)
            //This should in turn update the shifts in the Planar DislocationNode container
            VectorDim shift(VectorDim::Zero());
            // const VectorDim oldPosition(nodeIn->get_P()-nodeIn->get_V()*dt_in);
            // const VectorDim newPosition(nodeIn->get_P());

            for (const auto& gpPatch : periodicLoop.periodicGlidePlane->patches())
            {
                if ((gpPatch.second->shift - Loopshift).squaredNorm() < FLT_EPSILON)
                {
                    for (const auto &pEdges : gpPatch.second->edges())
                    {
                        SegmentSegmentDistance<DislocationNetworkType::dim> ssd(oldPosition, newPosition, pEdges->meshIntersection->P0, pEdges->meshIntersection->P1);
                        if (ssd.dMin < FLT_EPSILON)
                        {
                            //Intersection exist
                            shift += periodicLoop.periodicGlidePlane->getShift(*pEdges);
                        }
                    }
                }
            }
            assert(rveNodes.find(nodeIn)!=rveNodes.end());
            if (shift.squaredNorm()>FLT_EPSILON)
            {
                VectorDim shiftRVE(rveNodes.find(nodeIn)->second);
                rveNodes.find(nodeIn)->second+=shift;
                //update the 2D node shift in the rveNODE
                const auto periodicLoopIter(nodeIn->periodicNodeMap.find(&periodicLoop));
                assert(periodicLoopIter!=nodeIn->periodicNodeMap.end());
                const auto shiftIter(periodicLoopIter->second.find(shiftRVE));
                assert(shiftIter!=periodicLoopIter->second.end());
                assert(shiftIter->second.get()==this);
                const auto PNodetemp(shiftIter->second);
                periodicLoopIter->second.erase(shiftIter);
                shiftRVE+=shift;
                periodicLoopIter->second.emplace(shiftRVE,PNodetemp);
                //update finished
            }
            //Can be improved for removing the nodes 
        }

        VectorLowerDim get_P() const
        {
            assert(rveNodes.size());
            for(const auto& node : rveNodes)
            {
                const VectorDim shift(node.first->get_P()-periodicLoop.periodicGlidePlane->getGlobalPosition(P));
                // std::cout<<"shift is "<<shift.transpose()<<std::endl;
                LatticeVector<DislocationNetworkType::dim> latticeShift(periodicLoop.periodicGlidePlane->periodicGlidePlaneFactory.latticeVector(shift));
            }
            return P;
        }
        
        std::shared_ptr<NodeType> getRVESharedNode(const VectorDim& shift)
        {
            std::set< std::shared_ptr<NodeType>> temp;
            for(const auto& node : rveNodes)
            {
                if((node.second-shift).squaredNorm()<FLT_EPSILON)
                {
                    assert(periodicLoop.temporaryRVEsharedNodes.find(node.first)!=periodicLoop.temporaryRVEsharedNodes.end());
                    temp.insert(periodicLoop.temporaryRVEsharedNodes.find(node.first)->second);
                }
            }
         
            // std::cout<<"Temp size is "<<temp.size()<<std::endl;
            assert((temp.size()==1 || temp.size()==0));
            if (temp.size()==1)
            {
                return *temp.begin();
            }
            else
            {
                return nullptr;
            }
            
        }
        
        /**********************************************************************/
        PeriodicDislocationNode(const PeriodicDislocationNode&) = delete;
        const PeriodicDislocationNode& operator=(const PeriodicDislocationNode&) = delete;
        
        
        PeriodicNeighborConnectivity<DislocationNetworkType>& neighbors()
        {
            return *this;
        }

        const PeriodicNeighborConnectivity<DislocationNetworkType>& neighbors() const
        {
            return *this;
        }

        
        const PeriodicLoopLinkType* next(const PeriodicLoopLinkType* const link) const
        {
            VectorDim localBurgers(VectorDim::Zero());
            if(link->source.get()==this)
            {// out link
                localBurgers=link->loopLink->loop()->burgers();
            }
            else if(link->sink.get()==this)
            {// in link
                localBurgers=-link->loopLink->loop()->burgers();
            }
            else
            {
                assert(false && "node neither source nor sink");
            }
            
            std::set<const PeriodicDislocationNodeType*> nextNodes;
            for(const auto& neighbor : neighbors())
            {
                if((neighbor.second.first+localBurgers).squaredNorm()<FLT_EPSILON && neighbor.second.second.size()==1)
                {
                    nextNodes.insert(neighbor.first);
                }
            }
            
            if(nextNodes.size()==1)
            {
                return *neighbors().at(*nextNodes.begin()).second.begin();
            }
            else
            {// follow loop of link
                return link->next;
            }
            
        }
        /**********************************************************************/
        //Added by Yash
        LoopLinkType* getNeighborLinkinDifferentLoop(const PeriodicLoopLinkType &currentEdge, const size_t loopID)
        {
            //Grab only if neighbor belongs to a differnet loop
            const PeriodicLoopLinkType *temp(this->next(&currentEdge));
            if (temp->loopLink->loop()->sID != loopID)
            {
                return (temp->loopLink->loop()->linkStartingAt(temp->loopLink->source()->sID).second);
            }
            return nullptr;
        }
        /**********************************************************************/
        void addPeriodicLoopLink(PeriodicLoopLinkType *const periodicLoopLink)
        {
            if (periodicLoopLink->sink->sID == this->sID)
            { // edge ends at this node, so link is an inLink
                VerbosePeriodicDislocationBase(5,"Node "<<this->sID<<" ("<<periodicLoopLink->loopLink->sink()->sID<<"), loop "<<periodicLoopLink->loopLink->loop()->sID<<" adding inLink "<<periodicLoopLink<<std::endl;);
                // Update _loopConnectivities
                SuperNodalConnectivityType& loopConnectivity(loopConnectivities()[periodicLoopLink->loopLink->loop()->sID]);
                // std::cout<<"loopConnectivity.inEdge tag "<<loopConnectivity.inEdge->tag()<<std::endl;
                if (loopConnectivity.inEdge)
                {
                    std::cout<<"source->sink being inserted "<<periodicLoopLink->loopLink->source()->get_P().transpose()<<"-->"<<periodicLoopLink->loopLink->sink()->get_P().transpose()<<std::endl;
                    std::cout<<"source->sink already present "<<loopConnectivity.inEdge->loopLink->source()->get_P().transpose()<<"-->"<<loopConnectivity.inEdge->loopLink->sink()->get_P().transpose()<<std::endl;
                    std::cout<<"loopConnectivity.inEdge RVE tag "<<loopConnectivity.inEdge->loopLink->tag()<<std::endl;
                }

                assert(loopConnectivity.inEdge == nullptr || loopConnectivity.inEdge == periodicLoopLink);
                loopConnectivity.inEdge = periodicLoopLink;
                
                if (loopConnectivity.outEdge)
                {
                    // and outEdge of the same loop is connected to this node
                    assert(loopConnectivity.outEdge->prev == nullptr || loopConnectivity.outEdge->prev==periodicLoopLink);
                    assert(periodicLoopLink->next == nullptr || periodicLoopLink->next==loopConnectivity.outEdge);
                    loopConnectivity.outEdge->prev = periodicLoopLink;
                    periodicLoopLink->next = loopConnectivity.outEdge;
                }
                
            }
            else if (periodicLoopLink->source->sID == this->sID)
            {// edge starts at this node, so link is an outLink
                VerbosePeriodicDislocationBase(5,"Node "<<this->sID<<" ("<<periodicLoopLink->loopLink->source()->sID<<"), loop "<<periodicLoopLink->loopLink->loop()->sID<<" adding outLink "<<periodicLoopLink<<std::endl;);
                // Update _loopConnectivities
                SuperNodalConnectivityType& loopConnectivity(loopConnectivities()[periodicLoopLink->loopLink->loop()->sID]);
                assert(loopConnectivity.outEdge == nullptr || loopConnectivity.outEdge == periodicLoopLink);
                loopConnectivity.outEdge = periodicLoopLink;
                
                if (loopConnectivity.inEdge)
                {
                    // and outEdge of the same loop is connected to this node
                    assert(loopConnectivity.inEdge->next == nullptr || loopConnectivity.inEdge->next==periodicLoopLink);
                    assert(periodicLoopLink->prev == nullptr || periodicLoopLink->prev==loopConnectivity.inEdge);
                    loopConnectivity.inEdge->next = periodicLoopLink;
                    periodicLoopLink->prev = loopConnectivity.inEdge;
                }
            }
            else
            {
                assert(false && "CONNECTING LINK TO NON_INCIDENT NODE");
            }
            neighbors().addPeriodicLoopLink(periodicLoopLink);
        }
        
        /**********************************************************************/
        void removePeriodicLoopLink(PeriodicLoopLinkType *const periodicLoopLink)
        {
            
            if (periodicLoopLink->sink->sID == this->sID)
            {// edge ends at this node, so link is an inLink
                VerbosePeriodicDislocationBase(5,"Node "<<this->sID<<" ("<<periodicLoopLink->loopLink->sink()->sID<<"), loop "<<periodicLoopLink->loopLink->loop()->sID<<" removing inLink "<<periodicLoopLink<<std::endl;);

                // Update _loopConnectivities
                auto loopIter(loopConnectivities().find(periodicLoopLink->loopLink->loop()->sID));
                assert(loopIter != loopConnectivities().end() && "LOOP NOT FOUND IN loopConnectivities");
                SuperNodalConnectivityType& loopConnectivity(loopIter->second);
                assert(loopConnectivity.inEdge == periodicLoopLink);
                loopConnectivity.inEdge = nullptr;
                if (loopConnectivity.inEdge == nullptr && loopConnectivity.outEdge == nullptr)
                {
                    loopConnectivities().erase(loopIter);
                }
            }
            else if (periodicLoopLink->source->sID == this->sID)
            {// edge starts at this node, so link is an outLink
                // Update _loopConnectivities
                VerbosePeriodicDislocationBase(5,"Node "<<this->sID<<" ("<<periodicLoopLink->loopLink->source()->sID<<"), loop "<<periodicLoopLink->loopLink->loop()->sID<<" removing outLink "<<periodicLoopLink<<std::endl;);

                auto loopIter(loopConnectivities().find(periodicLoopLink->loopLink->loop()->sID));
                assert(loopIter != loopConnectivities().end() && "LOOP NOT FOUND IN loopConnectivities");
                SuperNodalConnectivityType& loopConnectivity(loopIter->second);
                assert(loopConnectivity.outEdge == periodicLoopLink);
                loopConnectivity.outEdge = nullptr;
                if (loopConnectivity.inEdge == nullptr && loopConnectivity.outEdge == nullptr)
                {
                    loopConnectivities().erase(loopIter);
                }
            }
            else
            {
                assert(false && "DISCONNECTING LINK FROM NON-INCIDENT NODES");
            }
            neighbors().removePeriodicLoopLink(periodicLoopLink);

        }
        
        /**********************************************************************/
        const SuperNodalConnectivityContainerType &loopConnectivities() const
        {
            return _loopConnectivities;
        }
        
        /**********************************************************************/
        SuperNodalConnectivityContainerType &loopConnectivities()
        {
            return _loopConnectivities;
        }
        
    };
    
    /**********************************************************************/
    /**********************************************************************/
    template<typename DislocationNetworkType>
    struct PeriodicLoopLink : public PeriodicDislocationBase
    {
        typedef typename DislocationNetworkType::LoopType LoopType;
        typedef typename LoopType::LoopLinkType LoopLinkType;
        typedef PeriodicDislocationNode<DislocationNetworkType> PeriodicDislocationNodeType;
        typedef PeriodicDislocationLoop<DislocationNetworkType> PeriodicDislocationLoopType;
        
        PeriodicDislocationLoopType& periodicLoop;
        const LoopLinkType* const loopLink;
        std::shared_ptr<PeriodicDislocationNodeType> source;
        std::shared_ptr<PeriodicDislocationNodeType>   sink;
//        PeriodicDislocationNodeType* const source;
//        PeriodicDislocationNodeType* const   sink;

        PeriodicLoopLink<DislocationNetworkType> *next;
        PeriodicLoopLink<DislocationNetworkType> *prev;
        
        /**********************************************************************/
        PeriodicLoopLink(PeriodicDislocationLoopType& pLoop,
                         const LoopLinkType* const pLink) :
        /* init */ periodicLoop(pLoop)
        /* init */,loopLink(pLink)
        /* init */,source((loopLink->source()->getSharedPeriodicNode(periodicLoop,loopLink->loop()->periodicShift)))
        /* init */,  sink((loopLink->sink()->getSharedPeriodicNode(periodicLoop,loopLink->loop()->periodicShift)))
        /* init */, next(nullptr)
        /* init */, prev(nullptr)
        {
            VerbosePeriodicDislocationBase(2,"Creating PeriodicLoopLink "<<loopLink->tag()<<std::endl;);

            source->addPeriodicLoopLink(this);
            sink->addPeriodicLoopLink(this);
        }
        
        /**********************************************************************/
        PeriodicLoopLink(const PeriodicLoopLink&) = delete;
        const PeriodicLoopLink& operator=(const PeriodicLoopLink&) = delete;

        /**********************************************************************/
        ~PeriodicLoopLink()
        {
            VerbosePeriodicDislocationBase(2,"Destroying PeriodicLoopLink "<<loopLink->tag()<<std::endl;);
            source->removePeriodicLoopLink(this); // Problem is here: both source and sink will try to modify untwinned
            sink->removePeriodicLoopLink(this);
            if (next)
            {
                next->prev = nullptr;
            }
            if (prev)
            {
                prev->next = nullptr;
            }
        }
    };
    
    /**********************************************************************/
    /**********************************************************************/
    template<typename DislocationNetworkType>
    class PeriodicDislocationLoop : public PeriodicDislocationBase
    /*                           */,public StaticID<PeriodicDislocationLoop<DislocationNetworkType>>
    /*                           */,public std::map<typename DislocationNetworkType::LoopLinkType*,PeriodicLoopLink<DislocationNetworkType>>
    /*                           */,private std::map<Eigen::Matrix<double, DislocationNetworkType::dim-1, 1>, const std::weak_ptr<PeriodicDislocationNode<DislocationNetworkType>>, CompareVectorsByComponent<double,DislocationNetworkType::dim-1,float>>
    /*                           */,private std::set<const PeriodicLoopLink<DislocationNetworkType>*>
    {
        static constexpr int dim=DislocationNetworkType::dim;
        typedef Eigen::Matrix<double,DislocationNetworkType::dim, 1> VectorDim;
        typedef Eigen::Matrix<double,DislocationNetworkType::dim - 1, 1> VectorLowerDim;
        typedef typename DislocationNetworkType::LoopType LoopType;
        typedef typename LoopType::LoopLinkType LoopLinkType;
        typedef PeriodicDislocationLoop<DislocationNetworkType> PeriodicDislocationLoopType;
        typedef PeriodicDislocationLoopFactory<DislocationNetworkType> PeriodicDislocationLoopFactoryType;
        typedef PeriodicDislocationNode<DislocationNetworkType> PeriodicDislocationNodeType;
        typedef PeriodicLoopLink<DislocationNetworkType> PeriodicLoopLinkType;
        typedef std::map<typename DislocationNetworkType::LoopLinkType*,PeriodicLoopLinkType> PeriodicLoopLinkContainerType;
        typedef std::map<Eigen::Matrix<double, DislocationNetworkType::dim - 1, 1>, const std::weak_ptr<PeriodicDislocationNodeType>, CompareVectorsByComponent<double, DislocationNetworkType::dim - 1, float>> NodeContainerType;
        typedef std::set<const PeriodicLoopLinkType*> UntwinnedEdgeContainerType;
        typedef GlidePlane<dim> GlidePlaneType;
        typedef typename GlidePlaneType::KeyType GlidePlaneKeyType;
        typedef std::vector<const PeriodicDislocationNodeType*> BoundaryContainerType;
        typedef std::vector<std::pair<BoundaryContainerType,VectorDim>> BoundariesContainerType;
        typedef std::map<size_t,std::pair<std::shared_ptr<PeriodicDislocationNodeType>,std::shared_ptr<typename DislocationNetworkType::NodeType>>> RVEnodeMapType;
        typedef std::map<typename DislocationNetworkType::NodeType*,std::shared_ptr<typename DislocationNetworkType::NodeType>> sharedRVENodesContainer;
        BoundariesContainerType _outerBoundaries;
        
    public:
        PeriodicDislocationLoopFactoryType& periodicDislocationLoopFactory;
        std::shared_ptr<PeriodicGlidePlane<dim>> periodicGlidePlane;
        sharedRVENodesContainer temporaryRVEsharedNodes;
        /**********************************************************************/
        PeriodicDislocationLoop(PeriodicDislocationLoopFactoryType& pdlf,
                                const GlidePlaneKeyType& key_in) :
        /* init */ periodicDislocationLoopFactory(pdlf)
        /* init */,periodicGlidePlane(periodicDislocationLoopFactory.periodicGlidePlaneFactory.get(key_in))
        {
            VerbosePeriodicDislocationBase(1,"Creating PeriodicDislocationLoop "<<this->sID<<std::endl;);
        }
        
        /**********************************************************************/
        ~PeriodicDislocationLoop()
        {
            VerbosePeriodicDislocationBase(1,"Destroying PeriodicDislocationLoop "<<this->sID<<std::endl;);
        }
        
        PeriodicDislocationLoop(const PeriodicDislocationLoop&) = delete;
        const PeriodicDislocationLoop& operator=(const PeriodicDislocationLoop&) = delete;

        /**********************************************************************/
        NodeContainerType& sharedNodes()
        {//!\returns the PeriodicDislocationNode container in the current time step
            return *this;
        }
        
        /**********************************************************************/
        const NodeContainerType& sharedNodes() const
        {//!\returns the PeriodicDislocationNode container in the current time step
            return *this;
        }
        /**********************************************************************/
        BoundariesContainerType &outerBoundaries()
        {
            return _outerBoundaries;
        }
        
        const BoundariesContainerType &outerBoundaries() const
        {
            return _outerBoundaries;
        }
        
        /**********************************************************************/
        const UntwinnedEdgeContainerType &untwinnedEdges() const
        {
            return *this;
        }
        
        UntwinnedEdgeContainerType &untwinnedEdges()
        {
            return *this;
        }
        
        /**********************************************************************/
        const PeriodicLoopLinkContainerType& loopLinks() const
        {
            return *this;
        }
        PeriodicLoopLinkContainerType& loopLinks()
        {
            return *this;
            
        }
        /**********************************************************************/
        void addUntwinnedEdge(const PeriodicLoopLink<DislocationNetworkType> *link)
        {
            untwinnedEdges().insert(link);
        }
        /**********************************************************************/
        void removeUntwinnedEdge(const PeriodicLoopLink<DislocationNetworkType> *link)
        {
            const size_t erased(untwinnedEdges().erase(link));
            assert(erased == 1 && "COULD NOT ERASE LINK FROM BOUNDARYLINKS");
        }
        void addLoopLink(LoopLinkType* const link)
        {// this needs to be an attach, which returns a PeriodicLoopLink stored where???
            const auto success(loopLinks().emplace(std::piecewise_construct,std::forward_as_tuple(link),std::forward_as_tuple(*this,link)));
            assert(success.second && "could not insert LoopLink in PeriodicLoopLinkContainerType");
        }
        void removeLoopLink(LoopLinkType* const link)
        {//this needs to be a detach, which returs a nullptr stored where???
            const size_t erased(loopLinks().erase(link));
            assert(erased==1 && "could not erase LoopLink from PeriodicLoopLinkContainerType");
        }
        
        /**********************************************************************/
        std::shared_ptr<PeriodicDislocationNode<DislocationNetworkType>> getSharedNode(const VectorDim& point, const VectorDim& shift)
        {
            return getSharedNode(periodicGlidePlane->getLocalPosition(point,shift));
        }
        
        
        /**********************************************************************/
        std::shared_ptr<PeriodicDislocationNode<DislocationNetworkType>> getSharedNode(const VectorLowerDim& point)
        {
            const auto iter(sharedNodes().find(point));
            if (iter == sharedNodes().end())
            {   // point does not exist
                //                sharedNodes().emplace(point,newNode.get());
                return sharedNodes().emplace(point, std::shared_ptr<PeriodicDislocationNodeType>(new PeriodicDislocationNodeType(*this))).first->second.lock();
                //                return newNode;
            }
            else
            { // point exists
                if (iter->second.expired())
                { // node deleted elsewhere
                    sharedNodes().erase(iter);
                    return sharedNodes().emplace(point, std::shared_ptr<PeriodicDislocationNodeType>(new PeriodicDislocationNodeType(*this))).first->second.lock();
                }
                else
                {
                    return iter->second.lock();
                }
            }
        }
        
        void createNewBoundary(const VectorDim& refBurgers,const PeriodicLoopLinkType* currentEdge, UntwinnedEdgeContainerType &untwinnedCopy, BoundaryContainerType& temp)
        {
            const size_t erased(untwinnedCopy.erase(currentEdge));
            if (erased != 1)
            {
                VerbosePeriodicDislocationBase(5,"currentEdge= "<<std::flush<<currentEdge->loopLink->tag()<<std::endl;);
                VerbosePeriodicDislocationBase(5,"untwinnedCopy= "<<std::endl;);
                for(const auto& untwinned : untwinnedCopy)
                {
                    VerbosePeriodicDislocationBase(5,untwinned->loopLink->tag()<<std::endl;);
                }
                 assert(erased == 1 && "could not find link in untwinnedEdges 2");
            }
            
            if((currentEdge->loopLink->loop()->burgers()-refBurgers).squaredNorm()<FLT_EPSILON)
            {// same Burgers, use sink
                temp.push_back(currentEdge->sink.get());
                if (currentEdge->sink->next(currentEdge)->sink.get() != temp.front())
                {
                    createNewBoundary(refBurgers,currentEdge->sink->next(currentEdge), untwinnedCopy,temp);
                }
            }
            else if((currentEdge->loopLink->loop()->burgers()+refBurgers).squaredNorm()<FLT_EPSILON)
            {// opposite Burgers, use source
                temp.push_back(currentEdge->source.get());
                if (currentEdge->source->next(currentEdge)->sink.get() != temp.front())
                {
                    createNewBoundary(refBurgers,currentEdge->source->next(currentEdge), untwinnedCopy,temp);
                }
            }
            else
            {
                assert(false && "currentEdge->burgers is not +/- refBurgers");
            }
        }
        /**********************************************************************/
        void updateOuterBoundaries()
        {
            VerbosePeriodicDislocationBase(2,"PeriodicDislocationLoop "<<this->sID<<" updating outer boundaries"<<std::endl;);
            outerBoundaries().clear();
            UntwinnedEdgeContainerType untwinnedCopy(this->untwinnedEdges());
            while (untwinnedCopy.size())
            {
                BoundaryContainerType temp;
                temp.reserve(untwinnedCopy.size());
                const VectorDim refBurgers((*untwinnedCopy.begin())->loopLink->loop()->burgers());
                createNewBoundary(refBurgers,*untwinnedCopy.begin(), untwinnedCopy,temp);
                if(temp.size()>=3)
                {
                    outerBoundaries().emplace_back(temp,refBurgers);
                }
            }
            VerbosePeriodicDislocationBase(2,", outerBoundaries().size() "<<outerBoundaries().size()<<std::endl;);
        }
        /**********************************************************************/
        void insertRVEloop(DislocationNetworkType& DN,
                           const std::vector<VectorLowerDim>& points,
                           const VectorDim& Burgers,
                           const std::shared_ptr<GlidePlaneType>& glidePlane,
                           const VectorDim& shift,
                           std::map<Eigen::Matrix<double, DislocationNetworkType::dim,1>, const std::shared_ptr<typename DislocationNetworkType::NodeType>, CompareVectorsByComponent<double,DislocationNetworkType::dim,float>>& rveNodesMap)
        {
            
            model::cout<<"        Initiating insertRVE loop with "<<points.size()<<" points ....."<<std::flush;
            std::vector<std::shared_ptr<typename DislocationNetworkType::NodeType>> nodes;
            for(const auto& point : points)
            {
                const auto sharedNode(getSharedNode(point));
                if(sharedNode->rveNodes.size())
                {// 2D has at least one rve node
                    auto rveNode(sharedNode->getRVESharedNode(shift)); // get appropriate node for shift
                    if (rveNode)
                    {
                        //rve node found
                        rveNode->confinedObject().clear();
                        rveNode->meshFaces().clear();
                        for (const auto &neighbor : rveNode->neighbors())
                        { // node may be moving away from boundary. Clear mesh faces of connected links
                          //Node may be a junction node moving out of boundary. Clear the confined object of connected links.
                            // std::get<1>(neighbor.second)->meshFaces().clear();
                            std::get<1>(neighbor.second)->confinedObject().clear();  //Discuss this with Dr. po
                        }
                        const VectorDim globalP(periodicGlidePlane->getGlobalPosition(point) + shift);
                        const auto rveIter(rveNodesMap.find(globalP));
                        if (rveIter != rveNodesMap.end())
                        {
                            std::cout<<rveIter->second->sID<<"-->"<<rveNode->sID<<std::endl;
                            assert(rveIter->second == rveNode && "rveNodes at same position");
                        }
                        rveNode->set_P_WO_Confinement(globalP);
                        // std::cout<<"rveNode looplink size is "<<rveNode->loopLinks().size()<<std::endl;
                        rveNode->confinedObject().updateGeometry(rveNode->get_P());

                        rveNode->p_Simplex = rveNode->get_includingSimplex(rveNode->get_P(), rveNode->p_Simplex); // update including simplex
                        nodes.push_back(rveNode);
                    }
                    else
                    {
                        //need to create a new RVE node if a 3D RVE node is not already present at the location
                        const VectorDim globalP(periodicGlidePlane->getGlobalPosition(point)+shift);
                        const auto rveIter(rveNodesMap.find(globalP));
                        if (rveIter != rveNodesMap.end())
                        {
                            nodes.push_back(rveIter->second);
                        }
                        else
                        {
                            std::shared_ptr<typename DislocationNetworkType::NodeType> newNode(new typename DislocationNetworkType::NodeType(&DN, globalP, VectorDim::Zero(), 1.0));
                            nodes.push_back(newNode);
                            rveNodesMap.emplace(globalP, newNode);
                            temporaryRVEsharedNodes.emplace(newNode.get(), newNode);
                        }
                        // TODO: CHECK IF THIS NODE IS NEEDED TO BE ADDED TO THE CONTAINER
                    }
                }
                else
                {// no 2d node found, check if rve node exists
                    const VectorDim globalP(periodicGlidePlane->getGlobalPosition(point)+shift);
                    const auto rveIter(rveNodesMap.find(globalP));
                    if(rveIter!=rveNodesMap.end())
                    {
                        nodes.emplace_back(rveIter->second);
                    }
                    else
                    {
                        std::shared_ptr<typename DislocationNetworkType::NodeType> newNode(new typename DislocationNetworkType::NodeType(&DN,globalP,VectorDim::Zero(),1.0));
                        // temporaryRVEsharedNodes.emplace(newNode.get(),newNode); //see this
                        rveNodesMap.emplace(globalP,newNode);
                        nodes.push_back(newNode);
                    }

                }
            }
            model::cout<<"        Reinserting loop with "<<points.size()<<" points ("<<nodes.size()<<" nodes)"<<std::flush;
            DN.insertLoop(nodes,Burgers,glidePlane,shift);
            
        }
        
        void updateSharedNodeswithConstraints() //This function assembles shared node container with constraints
        {
            VerbosePeriodicDislocationBase(5,"Updating shared nodes container for periodic dislocation loop with Constraints "<<this->sID<<std::flush;);
            
            sharedNodes().clear(); // since we search 2D nodes by position, we clear old positions from previous time step
            for(const auto& pair : loopLinks())
            {
                if (!(pair.first->source()->isBoundaryNode() && pair.second.source->loopConnectivities().size()>1))
                {
                    const auto temp(sharedNodes().emplace(pair.second.source->get_P(), pair.second.source));
                    if (!temp.second)
                    {
                        std::cout << "WARNING: could not add node to sharedNodes" << std::endl;
                    }
                }
            }
            VerbosePeriodicDislocationBase(5,"DONE [ "<<sharedNodes().size()<<" ] shared nodes "<<std::endl;);
        }

        void updateSharedNodeswithoutConstraints()
        {
            VerbosePeriodicDislocationBase(5,"Updating shared nodes container for periodic dislocation loop without constraints "<<this->sID<<std::flush;);
            
            sharedNodes().clear(); // since we search 2D nodes by position, we clear old positions from previous time step
            for(const auto& pair : loopLinks())
            {
                    const auto temp(sharedNodes().emplace(pair.second.source->get_P(), pair.second.source));
                    if (!temp.second)
                    {
                        std::cout << "WARNING: could not add node to sharedNodes" << std::endl;
                    }
            }
            VerbosePeriodicDislocationBase(5,"DONE [ "<<sharedNodes().size()<<" ] shared nodes "<<std::endl;);
        }

        void updateSharedRVENodes()
        {
            temporaryRVEsharedNodes.clear(); 
            for(const auto& pair : loopLinks())
            {
                temporaryRVEsharedNodes.emplace(pair.first->source().get(),pair.first->source());
            }
        }


        
        /**********************************************************************/
        void updateRVEloops(DislocationNetworkType& DN,
                            const double & dt_in,
                            std::map<Eigen::Matrix<double, DislocationNetworkType::dim,1>, const std::shared_ptr<typename DislocationNetworkType::NodeType>, CompareVectorsByComponent<double,DislocationNetworkType::dim,float>>& rveNodesMap
                            ,std::set<typename DislocationNetworkType::NodeType*>& nodesToBeMoved)
        {
            model::cout<<"        Moving DislocationNodes by glide (dt="<<dt_in<< ")"<<std::flush;
            // periodicGlidePlane->patches().clear();
            updateSharedNodeswithConstraints();
            updateSharedRVENodes();
//            RVEnodeMapType nodeMap(getNodeMap());
            const auto t6= std::chrono::system_clock::now();
            // for(const auto& nodePair : sharedNodes())
            // {
            //     if (!nodePair.second.expired())
            //     {
            //         for (const auto &rveNode : nodePair.second.lock()->rveNodes)
            //         {
            //             for (const auto &neighbor : rveNode.first->neighbors())
            //             { // node may be moving away from boundary. Clear mesh faces of connected links
            //                 std::get<1>(neighbor.second)->meshFaces().clear();
            //             }
            //             //Store old position and new position then update the position based on the updated shift
            //             // std::cout<<"Checking the size of the periodic patches "<<periodicGlidePlane->patches().size()<<std::endl;
            //             for (const auto& looptemp : rveNode.first->loops())
            //             {
            //                 //Update the shifts in the node
            //                 nodePair.second.lock()->updateShifts(rveNode.first,dt_in,looptemp->periodicShift);                            
            //             }
            //             rveNode.first->set_P_WO_Confinement(rveNode.first->get_P() + rveNode.first->get_V() * dt_in); // also changes 2D positions
            //         }
            //     }
            // }
            //Just print some data
            for (const auto &pair : loopLinks())
            {
                // if (!(pair.first->source()->isBoundaryNode() && pair.second.source->loopConnectivities().size() > 1))
                {
                    std::cout<<magentaColor<<"3D node ID and 2D position is "<<pair.first->source()->sID<<"--->"<<pair.second.source->sID<<"-->"
                    <<pair.second.source->get_P().transpose()<<"["<<pair.second.source->loopConnectivities().size()<<"]"
                    <<"[ "<<(pair.first->source()->isBoundaryNode()==0) <<"--->"<< (pair.first->source()->periodicNodeMap.find(this)->second.size()>1)<<" ]"
                    <<std::endl<<defaultColor;
                }
            }
            //After fixing twice motion of the nodes
            // for (const auto& nodetemp : DN.nodes())
            // {
            //      bool moveNode=true;
            //      for (const auto &neighbor : nodetemp.second->neighbors())
            //      { // node may be moving away from boundary. Clear mesh faces of connected links
            //          std::get<1>(neighbor.second)->meshFaces().clear();
            //      }
            //      //Store old position and new position then update the position based on the updated shift
            //      // std::cout<<"Checking the size of the periodic patches "<<periodicGlidePlane->patches().size()<<std::endl;
            //      for (const auto &looptemp : nodetemp.second->loops())
            //      {
            //          //Update the shifts in the node
            //          const auto& nodeShared(nodetemp.second->getSharedPeriodicNode(*this,looptemp->periodicShift));
            //          if (moveNode && sharedNodes().find(nodeShared->get_P()) != sharedNodes().end())
            //          {
            //              std::cout<<"Node shared loopConnectivity "<<nodeShared->loopConnectivities().size()<<std::endl;
            //              nodetemp.second->set_P_WO_Confinement(nodetemp.second->get_P() + nodetemp.second->get_V() * dt_in); // also changes 2D positions
            //              moveNode = false;
            //              updateSharedNodes(); //This will make things super slow
            //          }
            //          if (sharedNodes().find(nodeShared->get_P())!=sharedNodes().end())
            //          {
            //              std::cout<<"Node "<<nodetemp.second->sID<<" updating shift "<<std::endl; 
            //              nodeShared->updateShifts(nodetemp.second, dt_in, looptemp->periodicShift);
            //          }
            //      }
            //Fixing the motion for the third tim
            // for (const auto &nodetemp : DN.nodes())
            // {
            //     std::cout<<"Moving dislocation node "<<nodetemp.second->sID<<std::endl;
            //     bool moveNode = true;
            //     for (const auto &neighbor : nodetemp.second->neighbors())
            //     { // node may be moving away from boundary. Clear mesh faces of connected links
            //         std::get<1>(neighbor.second)->meshFaces().clear();
            //     }
            //     //Store old position and new position then update the position based on the updated shift
            //     // std::cout<<"Checking the size of the periodic patches "<<periodicGlidePlane->patches().size()<<std::endl;
            //     for (const auto &looptemp : nodetemp.second->loops())
            //     {
            //         std::cout<<"Updating for loop while moving "<<looptemp->sID<<"[ "<<looptemp->loopType<<" ]"<<std::endl;

            //         //Update the shifts in the node
            //         //Only possible if the loop is glissile
            //         if (looptemp->glidePlane!=nullptr)
            //         {
            //             std::cout<<"Updating for periodic loopID "<<this->sID<<std::endl;
            //             const auto &nodeShared(nodetemp.second->getSharedPeriodicNode(*this, looptemp->periodicShift));
            //             if (moveNode && sharedNodes().find(nodeShared->get_P()) != sharedNodes().end() &&
            //                 !(nodetemp.second->isBoundaryNode() && nodetemp.second->periodicNodeMap.find(this)->second.size() > 1)) //Condition important for self annihilatoin junction formation at bounday
            //             {
            //                 std::cout << nodetemp.second->sID << "Node shared loopConnectivity " << nodeShared->loopConnectivities().size() << std::endl;
            //                 nodetemp.second->set_P_WO_Confinement(nodetemp.second->get_P() + nodetemp.second->get_V() * dt_in); // also changes 2D positions
            //                 moveNode = false;
            //             }
            //         }
            //     }
                
            //     updateSharedNodeswithConstraints(); //This will make things super slow

            //     for (const auto &looptemp : nodetemp.second->loops())
            //     {
            //         std::cout<<"Updating for loop "<<looptemp->sID<<"[ "<<looptemp->loopType<<" ]"<<std::endl;
            //         if (looptemp->glidePlane!=nullptr)
            //         {
            //             const auto &nodeShared(nodetemp.second->getSharedPeriodicNode(*this, looptemp->periodicShift));

            //             if (sharedNodes().find(nodeShared->get_P()) != sharedNodes().end() &&
            //                 !(nodetemp.second->isBoundaryNode() && nodetemp.second->periodicNodeMap.find(this)->second.size() > 1))
            //             {
            //                 std::cout << "Node " << nodetemp.second->sID << " updating shift " << std::endl;
            //                 nodeShared->updateShifts(nodetemp.second, dt_in, looptemp->periodicShift);
            //             }
            //         }
            //         std::cout<<"Updated shift "<<std::endl;
            //     }
            // }
            //See what is the best location to achieve this
            //Fixing the motion for the fourth time
            for (const auto &nodetemp : temporaryRVEsharedNodes)
            {
                std::cout << "Moving dislocation node " << nodetemp.first->sID << std::endl;
                bool moveNode = nodesToBeMoved.find(nodetemp.first) != nodesToBeMoved.end();
                if (moveNode)
                {
                    const VectorDim oldPosition(nodetemp.first->snapToGlidePlanes(nodetemp.first->get_P()));
                    const VectorDim newPosition(nodetemp.first->snapToGlidePlanes(nodetemp.first->get_P()+nodetemp.first->get_V()*dt_in));

                    for (const auto &neighbor : nodetemp.first->neighbors())
                    { // node may be moving away from boundary. Clear mesh faces of connected links
                        std::get<1>(neighbor.second)->meshFaces().clear();
                    }
                    //Store old position and new position then update the position based on the updated shift
                    // std::cout<<"Checking the size of the periodic patches "<<periodicGlidePlane->patches().size()<<std::endl;
                    for (const auto &looptemp : nodetemp.first->loops())
                    {
                        std::cout << "Updating for loop while moving " << looptemp->sID << "[ " << looptemp->loopType << " ]" << std::endl;

                        //Update the shifts in the node
                        //Only possible if the loop is glissile
                        if (looptemp->glidePlane != nullptr)
                        {
                            std::cout << "Updating for periodic loopID " << this->sID << std::endl;
                            const auto &nodeShared(nodetemp.first->getSharedPeriodicNode(*looptemp->periodicLoop.get(), looptemp->periodicShift));
                            if (moveNode && sharedNodes().find(nodeShared->get_P()) != sharedNodes().end() &&
                                !(nodetemp.second->isBoundaryNode() && nodetemp.first->periodicNodeMap.find(looptemp->periodicLoop.get())->second.size() > 1)) //Condition important for self annihilatoin junction formation at bounday
                            {
                                std::cout << nodetemp.first->sID << "Node shared loopConnectivity " << nodeShared->loopConnectivities().size() << std::endl;
                                // nodetemp.first->set_P_WO_Confinement(nodetemp.first->get_P() + nodetemp.first->get_V() * dt_in); // also changes 2D positions
                                std::cout << nodetemp.first->sID << "Glide Plane size " << nodetemp.first->glidePlanes().size() << std::endl;
                                nodetemp.first->set_P_WO_Confinement(newPosition); // also changes 2D positions
                                const auto removedItems(nodesToBeMoved.erase(nodetemp.first));
                                assert(removedItems==1 && "Unable to update the containers of nodes to be moved"); 
                                moveNode = false; //Not to update multiple times
                            }
                        }
                    }

                    updateSharedNodeswithConstraints(); //This will make things super slow

                    for (const auto &looptemp : nodetemp.first->loops())
                    {
                        std::cout << "Updating for loop " << looptemp->sID << "[ " << looptemp->loopType << " ]" << std::endl;
                        if (looptemp->glidePlane != nullptr)
                        {
                            const auto &nodeShared(nodetemp.first->getSharedPeriodicNode(*looptemp->periodicLoop.get(), looptemp->periodicShift));

                            if (sharedNodes().find(nodeShared->get_P()) != sharedNodes().end() &&
                                !(nodetemp.first->isBoundaryNode() && nodetemp.first->periodicNodeMap.find(looptemp->periodicLoop.get())->second.size() > 1))
                            {
                                std::cout << "Node " << nodetemp.first->sID << " updating shift " << std::endl;
                                nodeShared->updateShifts(nodetemp.first, oldPosition, newPosition, looptemp->periodicShift);
                            }
                        }
                        std::cout << "Updated shift " << std::endl;
                    }
                }
            }
            model::cout << magentaColor << std::setprecision(3) << std::scientific << " [" << (std::chrono::duration<double>(std::chrono::system_clock::now() - t6)).count() << " sec]." << defaultColor << std::endl;

            updateSharedNodeswithoutConstraints();
            // DN.io().output(50);


            model::cout<<"        Updating outer boundaries "<<std::flush;
            const auto t7= std::chrono::system_clock::now();
            updateOuterBoundaries();
            model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<outerBoundaries().size()<<" boundaries ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t7)).count()<<" sec]."<<defaultColor<<std::endl;
            std::set<size_t> removeLoops;
            for(const auto& pair : loopLinks())
            {
                for (const auto &loopIter : pair.second.source->loopConnectivities())
                {
                    // std::cout << "Inserting " << loopIter.first << std::endl;
                    removeLoops.insert(loopIter.first);
                }
            }
            const auto t1= std::chrono::system_clock::now();
            model::cout<<"        Clipping outer boundaries"<<std::flush;
            typedef std::tuple<std::vector<VectorLowerDim>,VectorDim,std::shared_ptr<GlidePlaneType>,VectorDim> ReinsertTupleType;
            std::vector<ReinsertTupleType> reinsertVector;
            for (const auto& pair : outerBoundaries())
            {
//                std::vector<VectorDim> pgpPoints3D;
                std::vector<VectorLowerDim> pgpPoints2D;
                for (const auto& node : pair.first)
                {
//                     if(node->loopConnectivities().size()==1) // THERE IS A BUG HERE. IF MULTIPLE LOOPS IN SAME PATCH WE CAN HAVE loopConnectivities>1 and node NOT ON PATCH BOUNDARIES
//                     {// DO INSTEAD: collect all loops shifts of that node, and go in if statement if there is ony one
// //                        pgpPoints3D.emplace_back(periodicGlidePlane->getGlobalPosition(node->get_P()));
//                         pgpPoints2D.emplace_back(node->get_P());
//                     }
                    if(node->rveNodes.size()==1 &&
                    !(node->rveNodes.begin()->first->isBoundaryNode() && node->rveNodes.begin()->first->periodicNodeMap.find(this)->second.size() > 1)) // THERE WAS A BUG HERE. IF MULTIPLE LOOPS IN SAME PATCH WE CAN HAVE loopConnectivities>1 and node NOT ON PATCH BOUNDARIES
                    {// DO INSTEAD: collect all loops shifts of that node, and go in if statement if there is ony one
//                        pgpPoints3D.emplace_back(periodicGlidePlane->getGlobalPosition(node->get_P()));
                        pgpPoints2D.emplace_back(node->get_P());
                    }
                }

                periodicGlidePlane->addPatchesContainingPolygon(pgpPoints2D);
                
//                for(const auto& patchPair : periodicGlidePlane->patches())
//                {
//                    // CLIPPING
//                    std::vector<VectorLowerDim> patch2DPositions;
//                    for(const auto& patchLink : patchPair.second->edges())
//                    {
//                        patch2DPositions.push_back(*patchLink->source);
//                    }
//
//                    LoopPathClipper lpc(pgpPoints2D, patch2DPositions);
//                    lpc.makePaths();
//                    std::vector<std::vector<VectorLowerDim>> result = lpc.getClippedPolygons();
//
//                    for(const auto& points :  lpc.getClippedPolygons())
//                    {
//                        reinsertVector.emplace_back(points,pair.second,patchPair.second->glidePlane,patchPair.second->shift);
//                    }
//                }
                if (periodicGlidePlane->patches().size() == 1)
                {// only one patch, there cannot be patch intersections. Just return the outer2DNodes
                    const auto &iter(periodicGlidePlane->patches().begin());
                    reinsertVector.emplace_back(pgpPoints2D, pair.second, iter->second->glidePlane, iter->second->shift);
                }
                else
                {
                    for (const auto &patchPair : periodicGlidePlane->patches())
                    {
                        // CLIPPING
                        std::vector<VectorLowerDim> patch2DPositions;
                        for (const auto &patchLink : patchPair.second->edges())
                        {
                            patch2DPositions.push_back(*patchLink->source);
                        }
                        
                        LoopPathClipper lpc(pgpPoints2D, patch2DPositions);
                        lpc.makePaths();
                        for (const auto &points : lpc.getClippedPolygons())
                        {
                            if(points.size()>2)
                            {
                                reinsertVector.emplace_back(points, pair.second, patchPair.second->glidePlane, patchPair.second->shift);
                            }
                        }
                    }
                }
            }
            model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t1)).count()<<" sec]."<<defaultColor<<std::endl;

            updateSharedRVENodes();
            const auto t2= std::chrono::system_clock::now();
            model::cout<<"        Removing "<<removeLoops.size()<<" loops"<<std::flush;
            for (const auto& loopsID : removeLoops)
            {// Delete existing RVE loops
                DN.deleteLoop(loopsID);
            }
            model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t2)).count()<<" sec]."<<defaultColor<<std::endl;

            // RE-INSERT IN RVE
            const auto t3= std::chrono::system_clock::now();
            model::cout<<"        Reinserting " <<reinsertVector.size()<<" loops."<<std::flush;
            for(const auto tup : reinsertVector)
            {// Re-insert in RVE
                insertRVEloop(DN,std::get<0>(tup),std::get<1>(tup),std::get<2>(tup),std::get<3>(tup),rveNodesMap);
            }
            model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t3)).count()<<" sec]."<<defaultColor<<std::endl;


            
//            periodicGlidePlane->patches().clear();
            DN.updateGeometry(0.0);
            temporaryRVEsharedNodes.clear();

        }
        
    };
    
    
}
#endif
