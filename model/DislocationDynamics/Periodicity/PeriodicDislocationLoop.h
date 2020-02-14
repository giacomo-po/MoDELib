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
        typedef std::array<long int, dim + 3> KeyType;
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
    struct PeriodicDislocationNode : public PeriodicDislocationBase
    /*                            */,public StaticID<PeriodicDislocationNode<DislocationNetworkType>>
    /*                            */,public Eigen::Matrix<double, DislocationNetworkType::dim-1, 1>
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
        const VectorLowerDim& P;
        
    private:
        
        SuperNodalConnectivityContainerType _loopConnectivities;
        
    public:
        

        /**********************************************************************/
        PeriodicDislocationNode(PeriodicDislocationLoopType& periodicLoop_in,
                                const VectorLowerDim &pos) :
        /* init */ VectorLowerDim(pos)
        /* init */,PeriodicNeighborConnectivity<DislocationNetworkType>(this)
        /* init */,periodicLoop(periodicLoop_in)
        /* init */,P(*this)
        {
            VerbosePeriodicDislocationBase(2,"Creating PeriodicDislocationNode "<<this->sID<<std::endl;);
        }
        
        /**********************************************************************/
        ~PeriodicDislocationNode()
        {
            VerbosePeriodicDislocationBase(2,"Destroying PeriodicDislocationNode "<<this->sID<<std::endl;);
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
        void addPeriodicLoopLink(PeriodicLoopLinkType *const periodicLoopLink)
        {
            if (periodicLoopLink->sink->sID == this->sID)
            { // edge ends at this node, so link is an inLink
                VerbosePeriodicDislocationBase(5,"Node "<<this->sID<<" ("<<periodicLoopLink->loopLink->sink()->sID<<"), loop "<<periodicLoopLink->loopLink->loop()->sID<<" adding inLink "<<periodicLoopLink<<std::endl;);
                // Update _loopConnectivities
                SuperNodalConnectivityType& loopConnectivity(loopConnectivities()[periodicLoopLink->loopLink->loop()->sID]);
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
        PeriodicLoopLink<DislocationNetworkType> *next;
        PeriodicLoopLink<DislocationNetworkType> *prev;
        
        /**********************************************************************/
        PeriodicLoopLink(PeriodicDislocationLoopType& pLoop,
                         const LoopLinkType* const pLink) :
        /* init */ periodicLoop(pLoop)
        /* init */,loopLink(pLink)
        /* init */,source(periodicLoop.getSharedNode(loopLink->source()->get_P(),loopLink->loop()->periodicShift))
        /* init */,  sink(periodicLoop.getSharedNode(loopLink->  sink()->get_P(),loopLink->loop()->periodicShift))
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
        
        /**********************************************************************/
        void updateSourceSink()
        {
            const auto tempSource(periodicLoop.getSharedNode(loopLink->source()->get_P(),loopLink->loop()->periodicShift));
            const auto   tempSink(periodicLoop.getSharedNode(loopLink->  sink()->get_P(),loopLink->loop()->periodicShift));
            if(tempSource->sID!=source->sID || tempSink->sID!=sink->sID)
            {
                source->removePeriodicLoopLink(this);
                sink->removePeriodicLoopLink(this);
                source=tempSource;
                sink=tempSink;
                source->addPeriodicLoopLink(this);
                sink->addPeriodicLoopLink(this);
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
//        typedef std::map<size_t,LoopType*> LoopContainerType;
        typedef PeriodicDislocationLoop<DislocationNetworkType> PeriodicDislocationLoopType;
//        typedef PeriodicLoopObserver<PeriodicDislocationLoopType> PeriodicDislocationLoopObserverType;
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

        
        BoundariesContainerType _outerBoundaries;

        
    public:
        
        PeriodicDislocationLoopFactoryType& periodicDislocationLoopFactory;
        std::shared_ptr<PeriodicGlidePlane<dim>> periodicGlidePlane;
        
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
        NodeContainerType &nodes()
        {
            return *this;
        }
        
        /**********************************************************************/
        const NodeContainerType &nodes() const
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
            return getSharedNode(periodicGlidePlane->getLocalPosition(point-shift));
        }
        
        
        /**********************************************************************/
        std::shared_ptr<PeriodicDislocationNode<DislocationNetworkType>> getSharedNode(const VectorLowerDim& point)
        {
            const auto iter(nodes().find(point));
            if (iter == nodes().end())
            {   // point does not exist
                //                nodes().emplace(point,newNode.get());
                return nodes().emplace(point, std::shared_ptr<PeriodicDislocationNodeType>(new PeriodicDislocationNodeType(*this,point))).first->second.lock();
                //                return newNode;
            }
            else
            { // point exists
                if (iter->second.expired())
                { // node deleted elsewhere
                    nodes().erase(iter);
                    return nodes().emplace(point, std::shared_ptr<PeriodicDislocationNodeType>(new PeriodicDislocationNodeType(*this,point))).first->second.lock();
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
                           const RVEnodeMapType& nodeMap,
                           const std::vector<VectorLowerDim>& points,
                           const VectorDim& Burgers,
                           const std::shared_ptr<GlidePlaneType>& glidePlane,
                           const VectorDim& shift)
        {
            
            std::vector<std::shared_ptr<typename DislocationNetworkType::NodeType>> nodes;
            for(const auto& point : points)
            {
                const auto sharedNode(getSharedNode(point));
                const auto iter(nodeMap.find(sharedNode->sID));
                if(iter!=nodeMap.end())
                {// 2D node found, grab corresponding RVE node
                    iter->second.second->confinedObject().clear();
                    iter->second.second->meshFaces().clear();
                    static_cast<typename DislocationNetworkType::NodeType::SplineNodeType* const>(iter->second.second.get())->set_P(periodicGlidePlane->getGlobalPosition(point)+shift);
                    iter->second.second->confinedObject().updateGeometry(iter->second.second->get_P());
                    
                    iter->second.second->p_Simplex=iter->second.second->get_includingSimplex(iter->second.second->get_P(),iter->second.second->p_Simplex); // update including simplex

                    
                    nodes.push_back(iter->second.second);
                }
                else
                {
                    nodes.emplace_back(new typename DislocationNetworkType::NodeType(&DN,periodicGlidePlane->getGlobalPosition(point)+shift,VectorDim::Zero(),1.0));
                }
            }
            model::cout<<"        Reinserting loop with "<<points.size()<<" points ("<<nodes.size()<<" nodes)"<<std::flush;
            DN.insertLoop(nodes,Burgers,glidePlane,shift);
            
        }
        
        
        RVEnodeMapType getNodeMap() const
        {
            RVEnodeMapType nodeMap;
            const auto t0= std::chrono::system_clock::now();
            model::cout<<"        Constructing perdiodic nodes map"<<std::flush;
            for(const auto& node : nodes())
            {
                if(!node.second.expired())
                {
                    const auto periodicNode(node.second.lock());
                    if(periodicNode->loopConnectivities().size()==1)
                    {// not a boundary node
                        const auto rveNodeSource(periodicNode->loopConnectivities().begin()->second.outEdge->loopLink->source());
                        const auto   rveNodeSink(periodicNode->loopConnectivities().begin()->second.inEdge->loopLink->  sink());
                        assert(rveNodeSource==rveNodeSink);
                        nodeMap.emplace(periodicNode->sID,std::make_pair(periodicNode,rveNodeSource));
                    }
                }
            }
            model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<defaultColor<<std::endl;

            return nodeMap;
        }
        
        /**********************************************************************/
        void updateRVEloops(DislocationNetworkType& DN,const double & dt_in)
        {
            
            periodicGlidePlane->patches().clear();

            
            RVEnodeMapType nodeMap(getNodeMap());

            model::cout<<"        Moving DislocationNodes by glide (dt="<<dt_in<< ")"<<std::flush;
            const auto t6= std::chrono::system_clock::now();
            for(const auto& nodePair : nodeMap)
            {
                for(const auto& neighbor : nodePair.second.second->neighbors())
                {// node may be moving away from boundary. Clear mesh faces of connected links
                    std::get<1>(neighbor.second)->meshFaces().clear();
                }
                static_cast<typename DislocationNetworkType::NodeType::SplineNodeType* const>(nodePair.second.second.get())->set_P(nodePair.second.second->get_P()+nodePair.second.second->get_V()*dt_in);
                nodePair.second.second->updatePeriodicLoopLinks();
            }
            model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t6)).count()<<" sec]."<<defaultColor<<std::endl;

            nodeMap=getNodeMap(); // update nodeMap to store new PeriodicDislocationNodes created during move

            model::cout<<"        Updating outer boundaries "<<std::flush;
            const auto t7= std::chrono::system_clock::now();
            updateOuterBoundaries();
            model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<outerBoundaries().size()<<" boundaries ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t7)).count()<<" sec]."<<defaultColor<<std::endl;

            std::set<size_t> removeLoops;
            for(const auto& node : nodes())
            {
                if(!node.second.expired())
                {
                    const auto periodicNode(node.second.lock());
                    for(const auto& loopIter : periodicNode->loopConnectivities())
                    {
                        removeLoops.insert(loopIter.first);
                    }
                }
            }

            
            const auto t1= std::chrono::system_clock::now();
            model::cout<<"        Clipping outer boundaries"<<std::flush;
            typedef std::tuple<std::vector<VectorLowerDim>,VectorDim,std::shared_ptr<GlidePlaneType>,VectorDim> ReinsertTupleType;
            std::vector<ReinsertTupleType> reinsertVector;
            for (const auto& pair : outerBoundaries())
            {
                std::vector<VectorDim> pgpPoints3D;
                std::vector<VectorLowerDim> pgpPoints2D;
                for (const auto& node : pair.first)
                {
                    if(node->loopConnectivities().size()==1)
                    {
                        pgpPoints3D.emplace_back(periodicGlidePlane->getGlobalPosition(node->P));
                        pgpPoints2D.emplace_back(node->P);
                    }
                }

                periodicGlidePlane->addPatchesContainingPolygon(pgpPoints3D);
                
                for(const auto& patchPair : periodicGlidePlane->patches())
                {
                    // CLIPPING
                    std::vector<VectorLowerDim> patch2DPositions;
                    for(const auto& patchLink : patchPair.second->edges())
                    {
                        patch2DPositions.push_back(*patchLink->source);
                    }

                    LoopPathClipper lpc(pgpPoints2D, patch2DPositions);
                    lpc.makePaths();
                    std::vector<std::vector<VectorLowerDim>> result = lpc.getClippedPolygons();
                    
                    for(const auto& points :  lpc.getClippedPolygons())
                    {
                        reinsertVector.emplace_back(points,pair.second,patchPair.second->glidePlane,patchPair.second->shift);
                    }
                }
            }
            model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t1)).count()<<" sec]."<<defaultColor<<std::endl;

            
            const auto t2= std::chrono::system_clock::now();
            model::cout<<"        Removing loops"<<std::flush;
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
                insertRVEloop(DN,nodeMap,std::get<0>(tup),std::get<1>(tup),std::get<2>(tup),std::get<3>(tup));
            }
            model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t3)).count()<<" sec]."<<defaultColor<<std::endl;


            
//            periodicGlidePlane->patches().clear();
            DN.updateGeometry(0.0);

        }
        
    };
    
    
}
#endif
