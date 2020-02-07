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
    class PeriodicDislocationLoop;
    
    template<typename DislocationNetworkType>
    struct PeriodicLoopLink;
    
    template <typename DislocationNetworkType>
    struct PeriodicDislocationLoopFactory;
    
    template <typename DislocationNetworkType>
    struct TypeTraits<PeriodicDislocationLoopFactory<DislocationNetworkType>>
    {
        typedef PeriodicDislocationLoop<DislocationNetworkType> ValueType;
        typedef std::array<long int, DislocationNetworkType::dim + 3> KeyType;
        typedef std::less<std::array<long int, DislocationNetworkType::dim + 3>> CompareType;
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
    {
        
        typedef typename DislocationNetworkType::NodeType NodeType;
        typedef typename DislocationNetworkType::LoopType LoopType;
        typedef typename LoopType::LoopLinkType LoopLinkType;
        typedef PeriodicLoopLink<DislocationNetworkType> PeriodicLoopLinkType;
        typedef Eigen::Matrix<double,DislocationNetworkType::dim - 1, 1> VectorLowerDim;
        typedef SuperNodalConnectivity<DislocationNetworkType> SuperNodalConnectivityType;
        typedef std::map<size_t, SuperNodalConnectivityType> SuperNodalConnectivityContainerType;
        typedef PeriodicDislocationLoop<DislocationNetworkType> PeriodicDislocationLoopType;

        PeriodicDislocationLoopType& periodicLoop;
        const VectorLowerDim& P;
        
    private:
        
        SuperNodalConnectivityContainerType _loopConnectivities;
        SuperNodalConnectivityContainerType _neighborConnectivities;
        
    public:
        

        
        PeriodicDislocationNode(PeriodicDislocationLoopType& periodicLoop_in,
                                const VectorLowerDim &pos) :
        /* init */ VectorLowerDim(pos)
        /* init */,periodicLoop(periodicLoop_in)
        /* init */,P(*this)
        {
            VerbosePeriodicDislocationBase(2,"Creating PeriodicDislocationNode "<<this->sID<<std::endl;);
        }
        
        ~PeriodicDislocationNode()
        {
            VerbosePeriodicDislocationBase(2,"Destroying PeriodicDislocationNode "<<this->sID<<std::endl;);
        }
        
        PeriodicDislocationNode(const PeriodicDislocationNode&) = delete;
        const PeriodicDislocationNode& operator=(const PeriodicDislocationNode&) = delete;
        
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
                
                //update _neighborConnectivities
                SuperNodalConnectivityType& neighborConnectivity(neighborConnectivities()[periodicLoopLink->source->sID]);
                assert(neighborConnectivity.inEdge == nullptr || neighborConnectivity.inEdge == periodicLoopLink);
                neighborConnectivity.inEdge = periodicLoopLink;
                if (neighborConnectivity.outEdge)
                {
                    // and outEdge of the same loop is connected to this node
                    assert(neighborConnectivity.outEdge->twin == nullptr || neighborConnectivity.outEdge->twin == periodicLoopLink);
                    assert(periodicLoopLink->twin == nullptr || periodicLoopLink->twin == neighborConnectivity.outEdge);
                    neighborConnectivity.outEdge->twin = periodicLoopLink;
                    periodicLoopLink->twin = neighborConnectivity.outEdge;
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
                
                //update _neighborConnectivities
                SuperNodalConnectivityType& neighborConnectivity(neighborConnectivities()[periodicLoopLink->sink->sID]);
                assert(neighborConnectivity.outEdge == nullptr || neighborConnectivity.outEdge == periodicLoopLink);
                neighborConnectivity.outEdge = periodicLoopLink;
                if (neighborConnectivity.inEdge)
                {
                    // and outEdge of the same loop is connected to this node
                    assert(neighborConnectivity.inEdge->twin == nullptr || neighborConnectivity.inEdge->twin == periodicLoopLink);
                    assert(periodicLoopLink->twin == nullptr || periodicLoopLink->twin == neighborConnectivity.inEdge);
                    neighborConnectivity.inEdge->twin = periodicLoopLink;
                    periodicLoopLink->twin = neighborConnectivity.inEdge;
                }
            }
            else
            {
                assert(false && "CONNECTING LINK TO NON_INCIDENT NODE");
            }
        }
        
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
                
                // Update _neighborConnectivities
                auto neighborIter(neighborConnectivities().find(periodicLoopLink->source->sID));
                assert(neighborIter != neighborConnectivities().end() && "NEIGHBOR NOT FOUND IN neighborConnectivity");
                SuperNodalConnectivityType& neighborConnectivity(neighborIter->second);
                assert(neighborConnectivity.inEdge == periodicLoopLink);
                neighborConnectivity.inEdge = nullptr;
                if (neighborConnectivity.inEdge == nullptr && neighborConnectivity.outEdge == nullptr)
                {
                    neighborConnectivities().erase(neighborIter);
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
                // Update _neighborConnectivities
                auto neighborIter(neighborConnectivities().find(periodicLoopLink->sink->sID));
                assert(neighborIter != neighborConnectivities().end() && "NEIGHBOR NOT FOUND IN neighborConnectivity");
                SuperNodalConnectivityType& neighborConnectivity(neighborIter->second);
                assert(neighborConnectivity.outEdge == periodicLoopLink);
                neighborConnectivity.outEdge = nullptr;
                if (neighborConnectivity.inEdge == nullptr && neighborConnectivity.outEdge == nullptr)
                {
                    neighborConnectivities().erase(neighborIter);
                }
            }
            else
            {
                assert(false && "DISCONNECTING LINK FROM NON-INCIDENT NODES");
            }
        }
        
        //*********************************************************************
        const SuperNodalConnectivityContainerType &loopConnectivities() const
        {
            return _loopConnectivities;
        }
        
        /**********************************************************************/
        SuperNodalConnectivityContainerType &loopConnectivities()
        {
            return _loopConnectivities;
        }
        
        /**********************************************************************/
        const SuperNodalConnectivityContainerType &neighborConnectivities() const
        {
            return _neighborConnectivities;
        }
        
        /**********************************************************************/
        SuperNodalConnectivityContainerType &neighborConnectivities()
        {
            return _neighborConnectivities;
        }
    };
    
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
        PeriodicLoopLink<DislocationNetworkType> *twin;
        
        
        PeriodicLoopLink(PeriodicDislocationLoopType& pLoop,
                         const LoopLinkType* const pLink) :
        /* init */ periodicLoop(pLoop)
        /* init */,loopLink(pLink)
        /* init */,source(periodicLoop.getSharedNode(loopLink->source()->get_P(),loopLink->loop()->periodicShift))
        /* init */,  sink(periodicLoop.getSharedNode(loopLink->  sink()->get_P(),loopLink->loop()->periodicShift))
        /* init */, next(nullptr)
        /* init */, prev(nullptr)
        /* init */, twin(nullptr)
        {
            VerbosePeriodicDislocationBase(2,"Creating PeriodicLoopLink "<<loopLink->tag()<<std::endl;);

            source->addPeriodicLoopLink(this);
            sink->addPeriodicLoopLink(this);
            if (twin)
            {
                periodicLoop.removeUntwinnedEdge(twin);
            }
            else
            {
                periodicLoop.addUntwinnedEdge(this);
            }
        }
        
        PeriodicLoopLink(const PeriodicLoopLink&) = delete;
        const PeriodicLoopLink& operator=(const PeriodicLoopLink&) = delete;

        
        ~PeriodicLoopLink()
        {
            VerbosePeriodicDislocationBase(2,"Destroying PeriodicLoopLink "<<loopLink->tag()<<std::endl;);
            source->removePeriodicLoopLink(this);
            sink->removePeriodicLoopLink(this);
            if (next)
            {
                next->prev = nullptr;
            }
            if (prev)
            {
                prev->next = nullptr;
            }
            if (twin)
            {
                twin->twin = nullptr;
                periodicLoop.addUntwinnedEdge(twin);
            }
            else
            {
                periodicLoop.removeUntwinnedEdge(this);
            }
        }
        
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
    
    template<typename DislocationNetworkType>
    class PeriodicDislocationLoop : public PeriodicDislocationBase
    /*                           */,public StaticID<PeriodicDislocationLoop<DislocationNetworkType>>
    /*                           */,public std::map<size_t,typename DislocationNetworkType::LoopType*>
    /*                           */,public std::map<typename DislocationNetworkType::LoopLinkType*,PeriodicLoopLink<DislocationNetworkType>>
    /*                           */,private std::map<Eigen::Matrix<double, DislocationNetworkType::dim-1, 1>, const std::weak_ptr<PeriodicDislocationNode<DislocationNetworkType>>, CompareVectorsByComponent<double,DislocationNetworkType::dim-1,float>>
    /*                           */,private std::set<const PeriodicLoopLink<DislocationNetworkType>*>
    {
        static constexpr int dim=DislocationNetworkType::dim;
        typedef Eigen::Matrix<double,DislocationNetworkType::dim, 1> VectorDim;
        typedef Eigen::Matrix<double,DislocationNetworkType::dim - 1, 1> VectorLowerDim;
        typedef typename DislocationNetworkType::LoopType LoopType;
        typedef typename LoopType::LoopLinkType LoopLinkType;
        typedef std::map<size_t,LoopType*> LoopContainerType;
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
                    typedef std::vector<const PeriodicLoopLinkType*> BoundariesContainerType;

//        PeriodicLoopObserver<PeriodicDislocationLoopType>* const observer;
        
        BoundariesContainerType _outerBoundaries;

        
    public:
        
        PeriodicDislocationLoopFactoryType& periodicDislocationLoopFactory;
        std::shared_ptr<PeriodicGlidePlane<dim>> periodicGlidePlane;
        
        
        PeriodicDislocationLoop(PeriodicDislocationLoopFactoryType& pdlf,
                                const GlidePlaneKeyType& key_in) :
        /* init */ periodicDislocationLoopFactory(pdlf)
        /* init */,periodicGlidePlane(periodicDislocationLoopFactory.periodicGlidePlaneFactory.get(key_in))
        {
            VerbosePeriodicDislocationBase(1,"Creating PeriodicDislocationLoop "<<this->sID<<std::endl;);
        }
        
        ~PeriodicDislocationLoop()
        {
            VerbosePeriodicDislocationBase(1,"Destroying PeriodicDislocationLoop "<<this->sID<<std::endl;);
        }
        
        PeriodicDislocationLoop(const PeriodicDislocationLoop&) = delete;
        const PeriodicDislocationLoop& operator=(const PeriodicDislocationLoop&) = delete;

        BoundariesContainerType &outerBoundaries()
        {
            return _outerBoundaries;
        }
        
        const BoundariesContainerType &outerBoundaries() const
        {
            return _outerBoundaries;
        }
        
        const UntwinnedEdgeContainerType &untwinnedEdges() const
        {
            return *this;
        }
        
        UntwinnedEdgeContainerType &untwinnedEdges()
        {
            return *this;
        }
        
        const PeriodicLoopLinkContainerType& loopLinks() const
        {
            return *this;
            
        }
        
        PeriodicLoopLinkContainerType& loopLinks()
        {
            return *this;
            
        }
        
        const LoopContainerType& loops() const
        {
            return *this;
        }
        
        
        /**********************************************************************/
        NodeContainerType &nodes()
        {
            return *this;
        }
        
        const NodeContainerType &nodes() const
        {
            return *this;
        }
        
        //        LoopContainerType& loops()
        //        {
        //            return *this;
        //        }
        
        void addUntwinnedEdge(const PeriodicLoopLink<DislocationNetworkType> *link)
        {
            untwinnedEdges().insert(link);
        }
        
        void removeUntwinnedEdge(const PeriodicLoopLink<DislocationNetworkType> *link)
        {
            const size_t erased(untwinnedEdges().erase(link));
            assert(erased == 1 && "COULD NOT ERASE LINK FROM BOUNDARYLINKS");
        }
        
        void addLoop(const LoopType* )
        {
            // add to LoopContainerType
            // update outLoop reconstrcuction
            //            periodicPlane.addPatch(loop->periodicShift);
        }
        
        void removeLoop(const LoopType* )
        {
            // femove from LoopContainerType
            // update outLoop reconstrcuction
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
        
        /**********************************************************************/
        void createNewBoundary(const PeriodicLoopLinkType* currentEdge, UntwinnedEdgeContainerType &untwinnedCopy)
        {
            //            std::cout<<"createNewBoundary"<<std::endl;

            BoundariesContainerType temp;
            temp.reserve(untwinnedCopy.size());
            while (true)
            {
                //                std::cout<<"currentEdge "<<currentEdge->tag()<<std::endl;
                temp.push_back(currentEdge);
                const size_t erased(untwinnedCopy.erase(currentEdge));
                if (erased != 1)
                {
                    
                    std::cout << "Trying to erase " << currentEdge->tag() << std::endl;
                    std::cout << "untwinnedCopy is" << std::endl;
                    for (const auto &edgePtr : untwinnedCopy)
                    {
                        std::cout << "    " << edgePtr->tag() << std::endl;
                    }
                    
                    assert(erased == 1 && "could not find link in untwinnedEdges 2");
                }
                if (temp.back()->sink->sID == temp.front()->source->sID)
                {
                    break;
                }
                currentEdge = currentEdge->next;
                while (currentEdge->twin)
                {
                    currentEdge = currentEdge->twin->next;
                }
            }
            if (temp.size())
            {
                this->_outerBoundaries.push_back(temp);
            }
        }
        
        /**********************************************************************/
        void updateBoundaries()
        {
            std::cout<<"Updating OuterBoundary"<<std::endl;
            outerBoundaries().clear();
            UntwinnedEdgeContainerType untwinnedCopy(this->untwinnedEdges());
            while (untwinnedCopy.size())
            {
                for (const auto &edge : untwinnedCopy)
                {
                    if (edge->source->outEdges().size() == 1 && edge->source->inEdges().size() == 1)
                    { // must start from node with only one edge in and one edge out
                        createNewBoundary(edge, untwinnedCopy);
                        break;
                    }
                }
            }
            
        }
    };
    
    
}
#endif
