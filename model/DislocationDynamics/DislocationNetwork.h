/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po             <gpo@ucla.edu>.
 * Copyright (C) 2011 by Benjamin Ramirez       <ramirezbrf@gmail.com>
 * Copyright (C) 2011 by Mamdouh Mohamed        <msm07d@fsu.edu>
 * Copyright (C) 2011 by Tamer Crsoby           <tamercrosby@gmail.com>
 * Copyright (C) 2011 by Can Erel               <canerel55@gmail.com>
 * Copyright (C) 2011 by Yinan Cui              <cuiyinan@ucla.edu>
 * Copyright (C) 2017 by Sabyasachi Chatterjee  <sabyasac@andrew.cmu.edu>
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

// valgrind --leak-check=full --show-leak-kinds=all ./DDomp

#ifndef model_DISLOCATIONNETWORK_H_
#define model_DISLOCATIONNETWORK_H_

#ifdef _MODEL_MPI_
#define _MODEL_DD_MPI_
#endif


#ifdef _OPENMP
#include <omp.h>
#define EIGEN_DONT_PARALLELIZE // disable Eigen Internal openmp Parallelization
#endif

#include <vector>
#include <chrono>
#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <Eigen/StdDeque>


#include <LoopNetwork.h>
#include <TerminalColors.h>
#include <DislocationNetworkTraits.h>
#include <DislocationNetworkComponent.h>
#include <DislocationNode.h>
#include <DislocationSegment.h>
#include <DislocationLoop.h>
#include <GlidePlaneObserver.h>
#include <DislocationNetworkRemesh.h>
#include <DislocationJunctionFormation.h>
#include <DislocationCrossSlip.h>
//#include <Material.h>
#include <DislocationNetworkIO.h>
#include <DislocationParticle.h>
#include <DislocationStress.h>
//#include <ParticleSystem.h>
#include <MPIcout.h>
//#include <SingleFieldPoint.h>
#include <DDtimeIntegrator.h>
#include <EqualIteratorRange.h>
//#include <BoundingLineSegments.h>
#include <GrainBoundaryTransmission.h>
//#include <GrainBoundaryDissociation.h>
#include <BVPsolver.h>
#include <Polycrystal.h>
#include <DislocationNodeContraction.h>
#include <EshelbyInclusion.h>
#include <TextFileParser.h>
//#include <DisplacementPoint.h>
#include <DefectiveCrystalParameters.h>
#include <ExternalLoadControllerBase.h>
//#include <ExternalLoadController.h>
#include <DislocationInjector.h>

#ifdef _MODEL_GREATWHITE_
#include <MooseSolution.h>
#endif

namespace model
{
    
    
    
    template <int _dim, short unsigned int _corder, typename InterpolationType>
    class DislocationNetwork : public GlidePlaneObserver<_dim> //public DislocationNetworkBase<_dim>, // must be first in inheritance tree
    /* base                 */,public LoopNetwork<DislocationNetwork<_dim,_corder,InterpolationType> >
    //    /* base                 */ public ParticleSystem<DislocationParticle<_dim> >,
    /* base                 */,public std::map<size_t,EshelbyInclusion<_dim>>
    {
        
        
    public:
        
        
        static constexpr int dim=_dim; // make dim available outside class
        static constexpr int corder=_corder; // make dim available outside class
        typedef DislocationNetwork<dim,corder,InterpolationType> DislocationNetworkType;
        typedef LoopNetwork<DislocationNetworkType> LoopNetworkType;
        typedef DislocationNode<dim,corder,InterpolationType> NodeType;
        typedef DislocationSegment<dim,corder,InterpolationType> LinkType;
        typedef DislocationNetworkComponent<NodeType,LinkType> DislocationNetworkComponentType;
        typedef NetworkComponent<NodeType,LinkType> NetworkComponentType;
        typedef Eigen::Matrix<double,dim,dim>	MatrixDimD;
        typedef Eigen::Matrix<double,dim,1>		VectorDim;
        typedef typename TypeTraits<DislocationNetworkType>::LoopType LoopType;
        typedef GlidePlaneObserver<dim> GlidePlaneObserverType;
        //        typedef DislocationParticle<_dim> DislocationParticleType;
        //        typedef typename DislocationParticleType::StressField StressField;
        //        typedef typename DislocationParticleType::DisplacementField DisplacementField;
        //        typedef ParticleSystem<DislocationParticleType> ParticleSystemType;
        //        typedef typename ParticleSystemType::SpatialCellType SpatialCellType;
        //        typedef SpatialCellObserver<DislocationParticleType,_dim> SpatialCellObserverType;
        typedef BVPsolver<dim,2> BvpSolverType;
        typedef typename BvpSolverType::FiniteElementType FiniteElementType;
        typedef typename FiniteElementType::ElementType ElementType;
        typedef typename LoopNetworkType::IsNodeType IsNodeType;
        typedef DislocationNetworkIO<DislocationNetworkType> DislocationNetworkIOType;
        typedef Polycrystal<dim> PolycrystalType;
        //        typedef ExternalLoadControllerBase<dim> ExternalLoadControllerType;
        typedef std::map<size_t,EshelbyInclusion<_dim>> EshelbyInclusionContainerType;
        typedef NetworkLinkObserver<LinkType> NetworkLinkObserverType;
        typedef typename NetworkLinkObserverType::LinkContainerType NetworkLinkContainerType;
        
#ifdef DislocationNucleationFile
#include DislocationNucleationFile
#endif
        
#ifdef _MODEL_GREATWHITE_
#include <DislocationNetworkGreatWhite.h>
#endif
        
    private:
        
        
        /**********************************************************************/
        void updateImageLoops()
        {
            if(   simulationParameters.simulationType==DefectiveCrystalParameters::PERIODIC_IMAGES
               || simulationParameters.simulationType==DefectiveCrystalParameters::PERIODIC_FEM)
            {
                const auto t0= std::chrono::system_clock::now();
                model::cout<<"        Updating image  loops "<<std::flush;
                
                assert(mesh.regions().size()==1);
                const auto& region(*mesh.regions().begin()->second);
                
                
                std::vector<std::tuple<std::vector<std::shared_ptr<NodeType>>,VectorDim,size_t,int>> imageLoopVector;

                
                for(const auto& link : this->links())
                {
                    if(    link.second->isBoundarySegment()
                       && !link.second->hasZeroBurgers()
                       && !link.second->hasPeriodicLoop()
                       //&& !link.second->sink->isVirtualBoundaryNode && !link.second->source->isVirtualBoundaryNode
                       )
                    {// First make the check whether the nodes connected to the sinks are not virtualNodes and they lie on the boundary
                    
                        
                        
                        //model::cout<<"Creating DislocationNode "<<nodeIDinFile<<" ("<<kk<<" of "<<evl.nodes().size()<<")"<<std::endl;
//                        const size_t nodeID=this->insertDanglingNode(node.P,node.V,node.velocityReduction).first->first;

                        
                        
//                        const VectorDim imageSourceP(NodeType::imagePosition(link.second->source.get(),link.second->meshFaces()));
//                        const VectorDim imageSinkP(NodeType::imagePosition(link.second->sink.get(),link.second->meshFaces()));
//
//                        std::shared_ptr<NodeType> sourceImage;
//                        
//                        // search for imageSourceP among images of either source or sink
//                        
//                        if()
//                        {// position imageSourceP found
//                            // reset sourceImage to that image node
//                        }
                        
                        std::shared_ptr<NodeType> sourceImage(new NodeType(this,link.second->source.get(),link.second));
                        std::shared_ptr<NodeType>   sinkImage(new NodeType(this,link.second->  sink.get(),link.second));

                        
                        link.second->source->imageNodeContainer.insert(sourceImage->sID,sourceImage);
                        link.second->sink->imageNodeContainer.insert(sinkImage->sID,sinkImage);

                        std::vector<std::shared_ptr<NodeType>> imageNodes;
                        imageNodes.push_back(sourceImage);
                        imageNodes.push_back(sinkImage);
                        imageNodes.push_back(link.second->  sink);
                        imageNodes.push_back(link.second->source);
                        
                        imageLoopVector.emplace_back(imageNodes,
                                                       VectorDim::Zero(),
                                                       region.regionID,
                                                       DislocationLoopIO<dim>::PERIODICLOOP);
                    }
                }
                
                size_t nLoops(0);
                for(const auto& tup : imageLoopVector)
                {// Insert the new virtual loops
                    this->insertLoop(std::get<0>(tup),std::get<1>(tup),std::get<2>(tup),std::get<3>(tup));
                    nLoops++;
                }
                
                model::cout<<"("<<nLoops<<" periodic loops)"<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<defaultColor<<std::endl;
                
            }
        }
        
        
        void createImageLoops()
        {
            if(   simulationParameters.simulationType==DefectiveCrystalParameters::PERIODIC_IMAGES
               || simulationParameters.simulationType==DefectiveCrystalParameters::PERIODIC_FEM)
            {
                
                const auto t0= std::chrono::system_clock::now();
                model::cout<<"        Creating image  loops "<<std::flush;

                assert(mesh.regions().size()==1);
                const auto& region(*mesh.regions().begin()->second);
                
                std::vector<std::tuple<std::vector<std::shared_ptr<NodeType>>,VectorDim,VectorDim>> imageLoopVector;

                
                for(const auto& link : this->links())
                {
                    if(    link.second->isBoundarySegment()
                       && !link.second->hasZeroBurgers()
                       && link.second->isGlissile()
                       )
                    {
                        if(!link.second->source->masterNode && !link.second->sink->masterNode)
                        {
                            std::set<size_t> key;
                            for(const auto& face : link.second->meshFaces())
                            {
                                key.insert(face->sID);
                            }
                            
                            const auto sourceIter(link.second->source->imageNodeContainer.find(key));
                            assert(sourceIter!=link.second->source->imageNodeContainer.end());
                            
                            const auto sinkIter(link.second->sink->imageNodeContainer.find(key));
                            assert(sinkIter!=link.second->sink->imageNodeContainer.end());
                            
                            const VectorDim bndNormal(link.second->bndNormal());
                            const VectorDim glideNormal((*link.second->glidePlanes().begin())->unitNormal);

                            
                            const VectorDim glideDir(bndNormal-bndNormal.dot(glideNormal)*glideNormal);
                            
                            const VectorDim centerNodePos(0.5*(sourceIter->second->get_P()+sinkIter->second->get_P())+glideDir.normalized()*10.0);
                            
                            
                            
                            std::shared_ptr<NodeType> centerNode(new NodeType(this,centerNodePos,VectorDim::Zero(),1.0));
                            
                            std::vector<std::shared_ptr<NodeType>> imageNodes;
                            imageNodes.push_back(sinkIter->second);
                            imageNodes.push_back(sourceIter->second);
                            imageNodes.push_back(centerNode);
                            
                            imageLoopVector.emplace_back(imageNodes,
                                                         link.second->burgers(),
                                                         glideNormal);
                        }
                        else
                        {
                            assert(false && "FINISH HERE");
                        }

                    }
                    
                    

                }
                
                size_t nLoops(0);
                for(const auto& tup : imageLoopVector)
                {// Insert the new virtual loops
                    this->insertLoop(std::get<0>(tup),std::get<1>(tup),std::get<2>(tup),std::get<0>(tup)[0]->get_P(),region.regionID);
                    nLoops++;
                }
                model::cout<<"("<<nLoops<<" image loops)"<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<defaultColor<<std::endl;

            }
        }
        
        /**********************************************************************/
        void updateVirtualBoundaryLoops()
        {
            
            
            if(   simulationParameters.simulationType==DefectiveCrystalParameters::FINITE_FEM
               || simulationParameters.simulationType==DefectiveCrystalParameters::PERIODIC_FEM
               || simulationParameters.simulationType==DefectiveCrystalParameters::PERIODIC_IMAGES)
            {
                //                    if(useVirtualExternalLoops)
                //                    {
                
                const auto t0= std::chrono::system_clock::now();
                model::cout<<"        Updating virtual boundary loops "<<std::flush;
                
                
                // First clean up outdated boundary loops
                std::set<size_t> removeLoops;
                for(const auto& loop : this->loops())
                {
                    if((loop.second->isVirtualBoundaryLoop() && loop.second->links().size()!=4) || loop.second->isPureVirtualBoundaryLoop())
                    {// clean up left over loops from topological operations
                        removeLoops.insert(loop.second->sID);
                    }
                }
                
                for(const size_t& loopID : removeLoops)
                {// Remove the virtual loops with ID in removeLoops
                    this->deleteLoop(loopID);
                }
                
                
                
                // Now reconstruct virtual boundary loops
                std::vector<std::tuple<std::vector<std::shared_ptr<NodeType>>,VectorDim,size_t,int>> virtualLoopVector;
                for(const auto& link : this->links())
                {
                    if(link.second->isBoundarySegment() && !link.second->hasZeroBurgers())
                    {
                        virtualLoopVector.emplace_back(std::vector<std::shared_ptr<NodeType>>{link.second->sink,link.second->source,link.second->source->virtualBoundaryNode(),link.second->sink->virtualBoundaryNode()},
                                                       link.second->burgers(),
                                                       (*link.second->grains().begin())->grainID,
                                                       DislocationLoopIO<dim>::VIRTUALLOOP);
                    }
                }
                
                for(const auto& tup : virtualLoopVector)
                {// Insert the new virtual loops
                    this->insertLoop(std::get<0>(tup),std::get<1>(tup),std::get<2>(tup),std::get<3>(tup));
                }
                
                //                    }
                //                break;
                model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<defaultColor<<std::endl;
                
            }
            
            
            //            switch (simulationParameters.simulationType)
            //            {
            //                case DefectiveCrystalParameters::FINITE_FEM:
            //                {
            ////                    if(useVirtualExternalLoops)
            ////                    {
            //
            //                        // First clean up outdated boundary loops
            //                        std::set<size_t> removeLoops;
            //                        for(const auto& loop : this->loops())
            //                        {
            //                            if((loop.second->isVirtualBoundaryLoop() && loop.second->links().size()!=4) || loop.second->isPureVirtualBoundaryLoop())
            //                            {// clean up left over loops from topological operations
            //                                removeLoops.insert(loop.second->sID);
            //                            }
            //                        }
            //
            //
            //
            //                        for(const size_t& loopID : removeLoops)
            //                        {// Remove the virtual loops with ID in removeLoops
            //                            this->deleteLoop(loopID);
            //                        }
            //
            //
            //                        // Now reconstruct virtual boundary loops
            //                        std::vector<std::tuple<std::vector<std::shared_ptr<NodeType>>,VectorDim,size_t,int>> virtualLoopVector;
            //                        for(const auto& link : this->links())
            //                        {
            //                            if(link.second->isBoundarySegment() && !link.second->hasZeroBurgers())
            //                            {
            //                                virtualLoopVector.emplace_back(std::vector<std::shared_ptr<NodeType>>{link.second->sink,link.second->source,link.second->source->virtualBoundaryNode(),link.second->sink->virtualBoundaryNode()},
            //                                                               link.second->burgers(),
            //                                                               (*link.second->grains().begin())->grainID,
            //                                                               DislocationLoopIO<dim>::VIRTUALLOOP);
            //                            }
            //                        }
            //
            //                        for(const auto& tup : virtualLoopVector)
            //                        {// Insert the new virtual loops
            //                            this->insertLoop(std::get<0>(tup),std::get<1>(tup),std::get<2>(tup),std::get<3>(tup));
            //                        }
            //
            ////                    }
            //                    break;
            //                }
            //
            //                case DefectiveCrystalParameters::PERIODIC_IMAGES:
            //                {
            //
            //                    std::vector<std::tuple<std::vector<std::shared_ptr<NodeType>>,VectorDim,VectorDim,VectorDim,size_t>> DislocationLoopVector;
            //                    for(const auto& link : this->links())
            //                    {
            //
            //                        if(link.second->isBoundarySegment() && !link.second->hasZeroBurgers() && !link.second->sink->isVirtualBoundaryNode && !link.second->source->isVirtualBoundaryNode)
            //                        {// First make the check whether the nodes connected to the sinks are not virtualNodes and they lie on the boundary
            //
            //
            //                            const auto virtualLink(this->link(std::min(link.second->sink->virtualBoundaryNode()->sID,link.second->source->virtualBoundaryNode()->sID),std::max(link.second->sink->virtualBoundaryNode()->sID,link.second->source->virtualBoundaryNode()->sID)));
            //
            //                            if(!virtualLink.first)
            //                            {// Checking whether the virtual nodes of the original link are already connected or not
            //                                //                        std::cout<<"Size of virtual boundary loop from source  ID"<<link.second->source->virtualBoundaryNode()->sID<<"is"<<link.second->source->virtualBoundaryNode()->loops().size()<<std::endl;
            //                                //                        std::cout<<"Size of loopLinks from source  ID"<<link.second->source->virtualBoundaryNode()->sID<<"is"<<link.second->source->virtualBoundaryNode()->loopLinks().size()<<std::endl;
            //                                //                        std::cout<<"Size of neighbors from source  ID"<<link.second->source->virtualBoundaryNode()->sID<<"is"<<link.second->source->virtualBoundaryNode()->neighbors().size()<<std::endl;
            //                                //                        std::cout<<"Size of virtual boundary loop from sink ID"<<link.second->sink->virtualBoundaryNode()->sID<<"is"<<link.second->sink->virtualBoundaryNode()->loops().size()<<std::endl;
            //                                //                        std::cout<<"Size of loopLinks from sink  ID"<<link.second->sink->virtualBoundaryNode()->sID<<"is"<<link.second->sink->virtualBoundaryNode()->loopLinks().size()<<std::endl;
            //                                //                        std::cout<<"Size of neighbors from sink  ID"<<link.second->sink->virtualBoundaryNode()->sID<<"is"<<link.second->sink->virtualBoundaryNode()->neighbors().size()<<std::endl;
            //
            //                                if(   link.second->  sink->virtualBoundaryNode()->loops().size()==0
            //                                   && link.second->source->virtualBoundaryNode()->loops().size()==0)
            //                                {// If the virtual nodes are not connected to any loops, create a midnode based on the two virtual nodes and form a loop between them
            //
            //                                    VectorDim N = (*link.second->loopLinks().begin())->loop()->glidePlane->unitNormal;
            //                                    // Direction along which the new node should be placed
            //                                    VectorDim X=(-link.second->source->virtualBoundaryNode()->get_P()+link.second->sink->virtualBoundaryNode()->get_P()).cross(N).normalized();
            //                                    // Make sure the node is placed inside the simulation domain
            //                                    if (!(X.dot(link.second->source->bndNormal().normalized()) > 0)) X = -X;
            //                                    X=X*(10);
            //                                    // Create a new node based on two virtual nodes
            //                                    std::shared_ptr<NodeType> midNode(new NodeType(this,X+0.5*(link.second->source->virtualBoundaryNode()->get_P()+link.second->sink->virtualBoundaryNode()->get_P())
            //                                                                                   ,0.5*(link.second->source->virtualBoundaryNode()->get_V()+link.second->sink->virtualBoundaryNode()->get_V()),1));
            //                                    //                                    std::cout<<"Creating the mid node"<<midNode->sID<<std::endl;
            //                                    // Create a container to store the nodes and the direction of the loop to be inserted
            //                                    DislocationLoopVector.emplace_back(std::vector<std::shared_ptr<NodeType>>{link.second->source->virtualBoundaryNode(),midNode,link.second->sink->virtualBoundaryNode()},
            //                                                                       link.second->burgers(),N,link.second->source->virtualBoundaryNode()->get_P(),
            //                                                                       (*link.second->grains().begin())->grainID);
            //                                    //                                    std::cout<<"Creating the dislocationloop"<<std::endl;
            //                                }
            //                                else if(   link.second->  sink->virtualBoundaryNode()->loops().size()!=0
            //                                        && link.second->source->virtualBoundaryNode()->loops().size()==0)
            //                                {
            //                                    // Loop over the links associated with the virtual node that has a loop
            //                                    std::cout<<"The node to which expansion is taking place is"<<link.second->source->virtualBoundaryNode()->sID<<std::endl;
            //                                    for(const auto& pair: link.second->sink->virtualBoundaryNode()->linksByLoopID())
            //                                    {
            //                                        VectorDim N1 = (*pair.second.begin())->loop()->glidePlane->unitNormal;
            //                                        VectorDim B1 = (*pair.second.begin())->loop()->burgers();
            //                                        VectorDim N = (*link.second->loopLinks().begin())->loop()->glidePlane->unitNormal;
            //                                        VectorDim B = (*link.second->loopLinks().begin())->loop()->burgers();
            //                                        if ((N-N1).norm()==FLT_EPSILON && (B-B1).norm()==FLT_EPSILON )
            //                                        {
            //                                            std::cout<<"New loop needs to be merged on source"<<std::endl; //TODO not checked for the slip system for node expansion for the periodic loops
            //                                        }
            //                                        // Find the link which is a boundary segment
            //                                        if ((*pair.second.begin())->pLink->isBoundarySegment())
            //                                        {
            //                                            // Change the link of the other segment by using expand
            //                                            this->expand((*pair.second.rbegin())->pLink->source->sID,(*pair.second.rbegin())->pLink->sink->sID,link.second->source->virtualBoundaryNode());
            //                                            break;
            //                                        }
            //                                        else
            //                                        {
            //                                            this->expand((*pair.second.begin())->pLink->source->sID,(*pair.second.begin())->pLink->sink->sID,link.second->source->virtualBoundaryNode());
            //                                            break;
            //                                        }
            //                                    }
            //
            //                                    // Second case - Expand the existing loop on the virtual node instead of creating a new loop
            //                                    //                                    // Check is made to see whether one of the node has any loop associated with it
            //                                    //                                    if(link.second->source->virtualBoundaryNode()->loops().size()==0)
            //                                    //                                    {
            //                                    //
            //                                    //                                    }
            //                                    //
            //                                    //                                    else
            //                                    //                                    {
            //                                    //
            //                                    //                                    }
            //
            //                                }
            //                                else if(   link.second->  sink->virtualBoundaryNode()->loops().size()==0
            //                                        && link.second->source->virtualBoundaryNode()->loops().size()!=0)
            //                                {
            //                                    // Same procedure done for the other case
            //                                    std::cout<<"The node to which expansion is taking place is"<<link.second->sink->virtualBoundaryNode()->sID<<std::endl;
            //                                    for(const auto& pair: link.second->source->virtualBoundaryNode()->linksByLoopID())
            //                                    {
            ////                                        VectorDim N1 = (*pair.second.begin())->loop()->rightHandedUnitNormal();
            ////                                        VectorDim B1 = (*pair.second.begin())->loop()->burgers();
            ////                                        VectorDim N = (*link.second->loopLinks().begin())->loop()->rightHandedUnitNormal();
            ////                                        VectorDim B = (*link.second->loopLinks().begin())->loop()->burgers();
            //                                        //                                    if ((N-N1).norm()==FLT_EPSILON && (B-B1).norm()==FLT_EPSILON )
            //                                        //                                    {
            //                                        //                                        std::cout<<"New loop needs to be merged on sink"<<std::endl;
            //                                        //                                    }
            //                                        if ((*pair.second.begin())->pLink->isBoundarySegment())
            //                                        {
            //                                            this->expand((*pair.second.rbegin())->pLink->source->sID,(*pair.second.rbegin())->pLink->sink->sID,link.second->sink->virtualBoundaryNode());
            //                                            break;
            //                                        }
            //                                        else
            //                                        {
            //                                            this->expand((*pair.second.begin())->pLink->source->sID,(*pair.second.begin())->pLink->sink->sID,link.second->sink->virtualBoundaryNode());
            //                                            break;
            //                                        }
            //                                    }
            //                                }
            //                                else
            //                                {
            //
            //                                    assert(0 && "CASE NOT CONSIDERED SO FAR");
            //
            //                                }
            //                            }
            //                        }
            //
            //
            //                    }
            //                    for(const auto& tup : DislocationLoopVector)
            //                    {// Insert the new loops
            //                        this->insertLoop(std::get<0>(tup),std::get<1>(tup),std::get<2>(tup),std::get<3>(tup),std::get<4>(tup));
            //                    }
            //
            //
            //                    break;
            //                }
            //
            //                default:
            //                {
            //                    break;
            //                }
            //            }
            
            
            
        }
        
        
        
    public:
        
        const DefectiveCrystalParameters& simulationParameters;
        const SimplicialMesh<dim>& mesh;
        const Polycrystal<dim>& poly;
        const std::unique_ptr<BVPsolver<dim,2>>& bvpSolver;
        const std::unique_ptr<ExternalLoadControllerBase<dim>>& externalLoadController;
        const std::vector<VectorDim,Eigen::aligned_allocator<VectorDim>>& periodicShifts;
        DislocationNetworkRemesh<DislocationNetworkType> networkRemesher;
        DislocationJunctionFormation<DislocationNetworkType> junctionsMaker;
        DislocationNodeContraction<DislocationNetworkType> nodeContractor;
        GrainBoundaryTransmission<DislocationNetworkType> gbTransmission;
        
        //        int timeIntegrationMethod;
        //        int maxJunctionIterations;
        //        long int runID;
        //        double totalTime;
        //        double dt;
        //        double vMax;
        //        size_t Nsteps;
        //        MatrixDimD _plasticDistortionFromVelocities;
        MatrixDimD _plasticDistortionFromAreas;
        MatrixDimD _plasticDistortionRateFromVelocities;
        MatrixDimD _plasticDistortionRateFromAreas;
        int ddSolverType;
        bool computeDDinteractions;
        int crossSlipModel;
        //        bool use_boundary;
        //        unsigned int use_bvp;
        //        bool useVirtualExternalLoops;
        //        bool use_externalStress;
        //        bool use_extraStraightSegments;
        //        std::deque<StressStraight<dim>,Eigen::aligned_allocator<StressStraight<dim>>> ssdeq;
        //        std::deque<StressStraight<dim>,Eigen::aligned_allocator<StressStraight<dim>>> straightSegmentsDeq;
        int  outputFrequency;
        bool outputBinary;
        bool outputGlidePlanes;
        //        bool outputSpatialCells;
        //        bool outputPKforce;
        bool outputElasticEnergy;
        bool outputMeshDisplacement;
        bool outputFEMsolution;
        bool outputDislocationLength;
        bool outputPlasticDistortion;
        bool outputPlasticDistortionRate;
        bool outputQuadraturePoints;
        bool outputLinkingNumbers;
        bool outputLoopLength;
        bool outputSegmentPairDistances;
        //        unsigned int _userOutputColumn;
        bool use_stochasticForce;
        double surfaceAttractionDistance;
        //        bool computePlasticDistortionRateFromVelocities;
        std::string folderSuffix;
        
        /**********************************************************************/
        DislocationNetwork(int& argc, char* argv[],
                           const DefectiveCrystalParameters& _simulationParameters,
                           const SimplicialMesh<dim>& _mesh,
                           const Polycrystal<dim>& _poly,
                           const std::unique_ptr<BVPsolver<dim,2>>& _bvpSolver,
                           const std::unique_ptr<ExternalLoadControllerBase<dim>>& _externalLoadController,
                           const std::vector<VectorDim,Eigen::aligned_allocator<VectorDim>>& _periodicShifts,
                           long int& runID) :
        /* init */ simulationParameters(_simulationParameters)
        /* init */,mesh(_mesh)
        /* init */,poly(_poly)
        /* init */,bvpSolver(_bvpSolver)
        /* init */,externalLoadController(_externalLoadController)
        /* init */,periodicShifts(_periodicShifts)
        /* init */,networkRemesher(*this)
        /* init */,junctionsMaker(*this)
        /* init */,nodeContractor(*this)
        /* init */,gbTransmission(*this)
        //        /* init */,timeIntegrationMethod(TextFileParser("inputFiles/DD.txt").readScalar<int>("timeIntegrationMethod",true))
        ///* init */,maxJunctionIterations(TextFileParser("inputFiles/DD.txt").readScalar<int>("maxJunctionIterations",true))
        //        /* init */,runID(TextFileParser("inputFiles/DD.txt").readScalar<int>("startAtTimeStep",true)),
        //        /* init */,totalTime(0.0),
        //        /* init */ dt(0.0),
        //        /* init */ vMax(0.0),
        //        /* init */ Nsteps(TextFileParser("inputFiles/DD.txt").readScalar<size_t>("Nsteps",true)),
        //        /* init */,_plasticDistortionFromVelocities(MatrixDimD::Zero())
        /* init */,_plasticDistortionFromAreas(MatrixDimD::Zero())
        /* init */,_plasticDistortionRateFromVelocities(MatrixDimD::Zero())
        /* init */,_plasticDistortionRateFromAreas(MatrixDimD::Zero())
        /* init */,ddSolverType(TextFileParser("inputFiles/DD.txt").readScalar<int>("ddSolverType",true))
        /* init */,computeDDinteractions(TextFileParser("inputFiles/DD.txt").readScalar<int>("computeDDinteractions",true))
        /* init */,crossSlipModel(TextFileParser("inputFiles/DD.txt").readScalar<int>("crossSlipModel",true))
        //        /* init */ use_boundary(true),
        //        /* init */ use_bvp(TextFileParser("inputFiles/DD.txt").readScalar<int>("use_bvp",true)),
        //        /* init */,useVirtualExternalLoops(TextFileParser("inputFiles/DD.txt").readScalar<int>("useVirtualExternalLoops",true))
        //        /* init */,use_externalStress(false)
        //        /* init */,use_extraStraightSegments(false)
        /* init */,outputFrequency(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputFrequency",true))
        /* init */,outputBinary(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputBinary",true))
        /* init */,outputGlidePlanes(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputGlidePlanes",true))
        //        /* init */ outputSpatialCells(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputSpatialCells",true))
        //        /* init */ outputPKforce(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputPKforce",true))
        /* init */,outputElasticEnergy(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputElasticEnergy",true))
        /* init */,outputMeshDisplacement(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputMeshDisplacement",true))
        /* init */,outputFEMsolution(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputFEMsolution",true))
        /* init */,outputDislocationLength(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputDislocationLength",true))
        /* init */,outputPlasticDistortion(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputPlasticDistortion",true))
        /* init */,outputPlasticDistortionRate(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputPlasticDistortionRate",true))
        /* init */,outputQuadraturePoints(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputQuadraturePoints",true))
        /* init */,outputLinkingNumbers(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputLinkingNumbers",true))
        /* init */,outputLoopLength(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputLoopLength",true))
        /* init */,outputSegmentPairDistances(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputSegmentPairDistances",true))
        //        /* init */ _userOutputColumn(3)
        /* init */,use_stochasticForce(TextFileParser("inputFiles/DD.txt").readScalar<int>("use_stochasticForce",true))
        /* init */,surfaceAttractionDistance(TextFileParser("inputFiles/DD.txt").readScalar<double>("surfaceAttractionDistance",true))
        //        /* init */,computePlasticDistortionRateFromVelocities(TextFileParser("inputFiles/DD.txt").readScalar<int>("computePlasticDistortionRateFromVelocities",true))
        /* init */,folderSuffix("")
        {
            
            // Some sanity checks
            //            assert(Nsteps>=0 && "Nsteps MUST BE >= 0");
            
            // Initialize static variables
            LinkType::initFromFile("inputFiles/DD.txt");
            NodeType::initFromFile("inputFiles/DD.txt");
            LoopType::initFromFile("inputFiles/DD.txt");
            DislocationNetworkComponentType::initFromFile("inputFiles/DD.txt");
            DislocationStressBase<dim>::initFromFile("inputFiles/DD.txt");
            DDtimeIntegrator<0>::initFromFile("inputFiles/DD.txt");
            DislocationCrossSlip<DislocationNetworkType>::initFromFile("inputFiles/DD.txt");
            int stochasticForceSeed=TextFileParser("inputFiles/DD.txt").readScalar<int>("stochasticForceSeed",true);
            if(stochasticForceSeed<0)
            {
                StochasticForceGenerator::init(std::chrono::system_clock::now().time_since_epoch().count());
            }
            else
            {
                StochasticForceGenerator::init(stochasticForceSeed);
            }
            
            if(argc>1)
            {
                folderSuffix=argv[1];
                //                std::cout<<"folderSuffix="<<folderSuffix<<std::endl;
            }
            //            ParticleSystemType::initMPI(argc,argv);
            
            //
            // Read Vertex and Edge information
            EVLio<dim> evl;
            evl.read(runID,folderSuffix);
            createVertices(evl);
            createEdges(evl);
            updatePlasticDistortionFromAreas(simulationParameters.dt);
            createEshelbyInclusions();
            
            // IO
            //            io().read("./","DDinput.txt",runID);
#ifdef _MODEL_MPI_
            // Avoid that a processor starts writing before other are reading
            MPI_Barrier(MPI_COMM_WORLD);
#endif
            // Initializing configuration
            moveGlide(0.0);	// initial configuration
        }
        
        /* readVertices *******************************************************/
        void createVertices(const EVLio<dim>& evl)
        {/*!Creates DislocationNode(s) based on the data read by the EVLio<dim>
          * object.
          */
            size_t kk(1);
            for (const auto& node : evl.nodes())
            {
                const size_t nodeIDinFile(node.sID);
                NodeType::set_count(nodeIDinFile);
                if(node.sID==node.masterID)
                {// a regular node is created
                    model::cout<<"Creating DislocationNode "<<nodeIDinFile<<" ("<<kk<<" of "<<evl.nodes().size()<<")"<<std::endl;
                    const size_t nodeID=this->insertDanglingNode(node.P,node.V,node.velocityReduction).first->first;
                    assert(nodeID==nodeIDinFile);
                }
                else
                {
                    model::cout<<"Creating DislocationNode "<<nodeIDinFile<<" ("<<kk<<" of "<<evl.nodes().size()<<"), virtual of "<<node.masterID<<std::endl;
                    const auto isNode(this->node(node.masterID));
                    assert(isNode.first);
                    isNode.second->resetVirtualBoundaryNode();
                }
                kk++;
            }
        }
        
        /* readEdges **********************************************************/
        void createEdges(const EVLio<dim>& evl)
        {/*!
          */
            
            std::map<size_t,std::map<size_t,size_t>> loopMap;
            for(const auto& looplink : evl.links())
            {// Collect LoopLinks by loop IDs
                loopMap[looplink.loopID].emplace(looplink.sourceID,looplink.sinkID);
            }
            
            assert(loopMap.size()==evl.loops().size());
            
            size_t loopLumber=1;
            for(const auto& loop : evl.loops())
            {// for each loop in the EVLio<dim> object
                
                const auto loopFound=loopMap.find(loop.sID); // there must be an entry with key loopID in loopMap
                assert(loopFound!=loopMap.end());
                std::vector<size_t> nodeIDs;
                nodeIDs.push_back(loopFound->second.begin()->first);
                for(size_t k=0;k<loopFound->second.size();++k)
                {
                    const auto nodeFound=loopFound->second.find(*nodeIDs.rbegin());
                    if(k<loopFound->second.size()-1)
                    {
                        nodeIDs.push_back(nodeFound->second);
                    }
                    else
                    {
                        assert(nodeFound->second==nodeIDs[0]);
                    }
                }
                
                LoopType::set_count(loop.sID);
                
                model::cout<<"Creating Dislocation Loop "<<loop.sID<<" ("<<loopLumber<<" of "<<evl.loops().size()<<"), type="<<loop.loopType<<std::endl;
                switch (loop.loopType)
                {
                    case DislocationLoopIO<dim>::GLISSILELOOP:
                    {
                        assert(loop.N.squaredNorm()>FLT_EPSILON && "GLISSILELOOP must have non-zero normal");
                        const size_t newLoopID=this->insertLoop(nodeIDs,loop.B,loop.N,loop.P,loop.grainID)->sID;
                        assert(loop.sID==newLoopID);
                        break;
                    }
                        
                    case DislocationLoopIO<dim>::SESSILELOOP:
                    {
                        const size_t newLoopID=this->insertLoop(nodeIDs,loop.B,loop.grainID,loop.loopType)->sID;
                        assert(loop.sID==newLoopID);
                        break;
                    }
                        
                    case DislocationLoopIO<dim>::VIRTUALLOOP:
                    {
                        std::vector<std::shared_ptr<NodeType>> sharedNodes;
                        for(const size_t nodeID : nodeIDs)
                        {// collect shared_ptrs to nodes
                            const auto isNode(this->node(nodeID));
                            assert(isNode.first);
                            if(isNode.second->masterNode)
                            {// a virtual node
                                sharedNodes.push_back(isNode.second->masterNode->virtualBoundaryNode());
                            }
                            else
                            {
                                const auto isSharedNode(this->danglingNode(nodeID));
                                if(!isSharedNode.first)
                                {
                                    model::cout<<"node "<<nodeID<<" not found"<<std::endl;
                                    assert(false && "node shared pointer not found");
                                }
                                sharedNodes.push_back(isSharedNode.second);
                            }
                        }
                        const size_t newLoopID=this->insertLoop(sharedNodes,loop.B,loop.grainID,loop.loopType)->sID;
                        assert(loop.sID==newLoopID);
                        break;
                    }
                        
                    default:
                        assert(false && "Unknown DislocationLoop type");
                        break;
                }
                
                //                if(loop.N.squaredNorm()>FLT_EPSILON)
                //                {// Regular loop
                ////                    model::cout<<"Creating Dislocation Loop "<<loop.sID<<" ("<<loopLumber<<" of "<<evl.loops().size()<<")"<<std::endl;
                //                    const size_t newLoopID=this->insertLoop(nodeIDs,loop.B,loop.N,loop.P,loop.grainID)->sID;
                //                    assert(loop.sID==newLoopID);
                //                }
                //                else
                //                {// Virtual loop. Note that virtual nodes are not in the danglingNodes, but they are instead store in their master nodes
                ////                    model::cout<<"Creating Dislocation Loop "<<loop.sID<<" ("<<loopLumber<<" of "<<evl.loops().size()<<") virtual"<<std::endl;
                //                    std::vector<std::shared_ptr<NodeType>> sharedNodes;
                //                    for(const size_t nodeID : nodeIDs)
                //                    {// collect shared_ptrs to nodes
                //                        const auto isNode(this->node(nodeID));
                //                        assert(isNode.first);
                //                        if(isNode.second->masterNode)
                //                        {// a virtual node
                //                            sharedNodes.push_back(isNode.second->masterNode->virtualBoundaryNode());
                //                        }
                //                        else
                //                        {
                //                            const auto isSharedNode(this->danglingNode(nodeID));
                //                            if(!isSharedNode.first)
                //                            {
                //                                model::cout<<"node "<<nodeID<<" not found"<<std::endl;
                //                                assert(false && "node shared pointer not found");
                //                            }
                //                            sharedNodes.push_back(isSharedNode.second);
                //                        }
                //                    }
                //                    const size_t newLoopID=this->insertLoop(sharedNodes,loop.B,loop.grainID,loop.loopType)->sID;
                //                    assert(loop.sID==newLoopID);
                //                }
                loopLumber++;
            }
            
            this->clearDanglingNodes();
            
            
            
            model::cout<<std::endl;
        }
        
        /**********************************************************************/
        void createEshelbyInclusions()
        {
            IDreader<'E',1,14,double> inclusionsReader;
            inclusionsReader.read(0,true);
            
            const std::vector<double> inclusionsMobilityReduction(TextFileParser("./inputFiles/initialMicrostructure.txt").readArray<double>("inclusionsMobilityReduction",true));
            for(const auto& pair : inclusionsReader)
            {
                
                const size_t& inclusionID(pair.first);
                Eigen::Map<const Eigen::Matrix<double,1,14>> row(pair.second.data());
                
                const VectorDim C(row.template segment<dim>(0));
                const double a(row(dim+0));
                MatrixDimD eT(MatrixDimD::Zero());
                const int typeID(row(13));
                int k=dim+1;
                for(int i=0;i<dim;++i)
                {
                    for(int j=0;j<dim;++j)
                    {
                        eT(i,j)=row(k);
                        k++;
                    }
                }
                
                
                
                EshelbyInclusion<dim>::set_count(inclusionID);
                eshelbyInclusions().emplace(std::piecewise_construct,
                                            std::make_tuple(inclusionID),
                                            std::make_tuple(C,a,eT,poly.nu,poly.mu,inclusionsMobilityReduction[typeID],typeID) );
            }
        }
        
        /**********************************************************************/
        bool remove(const size_t& nodeID)
        {
            const auto isNode=this->node(nodeID);
            if(isNode.first)
            {
                switch (simulationParameters.simulationType)
                {
                    case DefectiveCrystalParameters::FINITE_FEM:
                    {
                        if(isNode.second->virtualBoundaryNode())
                        {
                            LoopNetworkType::remove(isNode.second->virtualBoundaryNode()->sID);
                        }
                        break;
                    }
                        
                    case DefectiveCrystalParameters::PERIODIC_IMAGES:
                    {
                        // FINISH HERE
                        break;
                    }
                }
            }
            return LoopNetworkType::remove(nodeID);
        }
        
        /**********************************************************************/
        void updateGeometry(const double& dt)
        {
            for(auto& loop : this->loops())
            {// copmute slipped areas and right-handed normal // TODO: PARALLELIZE THIS LOOP
                loop.second->updateGeometry();
            }
            updatePlasticDistortionFromAreas(dt);
        }
        
        //        /**********************************************************************/
        //        double get_dt() const
        //        {
        //            switch (timeIntegrationMethod)
        //            {
        //                case 0:
        //                    return DDtimeIntegrator<0>::get_dt(*this);
        //                    break;
        //
        //                default:
        //                    assert(0 && "time integration method not implemented");
        //                    return 0;
        //                    break;
        //            }
        //        }
        
        /**********************************************************************/
        void removeZeroAreaLoops()
        {
            const auto t0= std::chrono::system_clock::now();
            model::cout<<"        Removing zero-area loops "<<std::flush;
            std::deque<size_t> loopIDs;
            for(const auto& loop : this->loops())
            {
                if(loop.second->slippedArea()<FLT_EPSILON)
                {
                    loopIDs.push_back(loop.second->sID);
                }
            }
            
            for(const auto& loopID : loopIDs)
            {
                this->deleteLoop(loopID);
            }
            model::cout<<"("<<loopIDs.size()<<" removed)"<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<defaultColor<<std::endl;
        }
        
        /**********************************************************************/
        void singleGlideStepDiscreteEvents(const long int& runID)
        {
            //! A simulation step consists of the following:
            //            model::cout<<blueBoldColor<< "runID="<<runID<<" (of "<<Nsteps<<")"
            //            /*                    */<< ", time="<<totalTime
            //            /*                    */<< ": nodes="<<this->nodes().size()
            //            /*                    */<< ", segments="<<this->links().size()
            //            /*                    */<< ", loopSegments="<<this->loopLinks().size()
            //            /*                    */<< ", loops="<<this->loops().size()
            //            /*                    */<< ", components="<<this->components().size()
            //            /*                    */<< defaultColor<<std::endl;
            
            //! 1- Check that all nodes are balanced
            //            checkBalance();
            
            //! 2 - Update quadrature points
            //            updateQuadraturePoints();
            //            updateStressStraightSegments();
            
            //            for(auto& loop : this->loops()) // TODO: PARALLELIZE THIS LOOP
            //            {// copmute slipped areas and right-handed normal
            //                loop.second->updateGeometry();
            //            }
            //            updatePlasticDistortionFromAreas();
            
            //! 3- Calculate BVP correction
            //            updateLoadControllers(runID);
            
            //#ifdef DislocationNucleationFile
            //            if(use_bvp && !(runID%use_bvp))
            //            {
            //                nucleateDislocations(); // needs to be called before updateQuadraturePoints()
            //                updateQuadraturePoints();
            //            }
            //#endif
            //            assembleAndSolve(runID,straightSegmentsDeq);
            //            computeNodaVelocities(runID);
            
            
            //! 4- Solve the equation of motion
            
            
            //! 5- Compute time step dt (based on max nodal velocity) and increment totalTime
            // make_dt();
            
            
            
            //! 6- Output the current configuration before changing it
            //            output(runID);
            //            io().output(runID);
            
            
            //! 7- Moves DislocationNodes(s) to their new configuration using stored velocity and dt
            //            move(dt);
            
            //! 8- Update accumulated quantities (totalTime and plasticDistortion)
            //            totalTime+=dt;
            //            updatePlasticDistortionRateFromVelocities();
            
            
            //! 9- Contract segments of zero-length
            //            DislocationNetworkRemesh<DislocationNetworkType>(*this).contract0chordSegments();
            
            //            if(runID>0)
            //            {
            //                removeZeroAreaLoops();
            //            }
            
            //! 10- Cross Slip (needs upated PK force)
            DislocationCrossSlip<DislocationNetworkType>(*this);
            
            
            //                        gbTransmission.transmit();
            //gbTransmission.directTransmit();
            
            //
            //            GrainBoundaryDissociation<DislocationNetworkType>(*this).dissociate();
            
            //            poly.reConnectGrainBoundarySegments(*this); // this makes stressGauss invalid, so must follw other GB operations
            
            
            //! 11- detect loops that shrink to zero and expand as inverted loops
            //            DislocationNetworkRemesh<DislocationNetworkType>(*this).loopInversion(dt);
            
            //! 12- Form Junctions
            junctionsMaker.formJunctions(DDtimeIntegrator<0>::dxMax);
            
            //            // Remesh may contract juncitons to zero lenght. Remove those juncitons:
            //            DislocationJunctionFormation<DislocationNetworkType>(*this).breakZeroLengthJunctions();
            
            
            
            //! 13- Node redistribution
            networkRemesher.remesh(runID);
            
            createImageLoops();
//            updateImageLoops();
            updateVirtualBoundaryLoops();
            
            //            mergeLoopsAtNodes();
            
            //            DislocationInjector<DislocationNetworkType>(*this).insertRandomStraightDislocation();
            
            //! 9- Contract segments of zero-length
            //            DislocationNetworkRemesh<DislocationNetworkType>(*this).contract0chordSegments();
            
            //! 14- If BVP solver is not used, remove DislocationSegment(s) that exited the boundary
            //            removeBoundarySegments();
            
            //            removeSmallComponents(3.0*dx,4);
            
            //            make_bndNormals();
            
            //! 16 - Increment runID counter
            //            ++runID;     // increment the runID counter
        }
        
        //        const VectorDim periodicVector(const Eigen::Array<int,dim,1>& cellID) const
        //        {
        //            return (meshDimensions.array()*cellID.template cast<double>()).matrix();
        //        }
        
        /**********************************************************************/
        const EshelbyInclusionContainerType& eshelbyInclusions() const
        {
            return *this;
        }
        
        EshelbyInclusionContainerType& eshelbyInclusions()
        {
            return *this;
        }
        
        
        /**********************************************************************/
        bool contract(std::shared_ptr<NodeType> nA,
                      std::shared_ptr<NodeType> nB)
        {
            return nodeContractor.contract(nA,nB);
        }
        
        
        
        
        
        /**********************************************************************/
        void assembleAndSolveGlide(const long int& runID)
        {/*! Performs the following operatons:
          */
#ifdef _OPENMP
            const size_t nThreads = omp_get_max_threads();
#else
            const size_t nThreads = 1;
#endif
            
            //! -1 Compute the interaction StressField between dislocation particles
            
            
            if(corder==0)
            {// For straight segments use analytical expression of stress field
                //                const auto t0= std::chrono::system_clock::now();
                //                model::cout<<"		Collecting StressStraight objects: "<<std::flush;
                //
                ////                std::deque<StressStraight<dim>,Eigen::aligned_allocator<StressStraight<dim>>> straightSegmentsDeq;
                //                size_t currentSize=0;
                //                if(computeDDinteractions)
                //                {
                //                    for(const auto& link : this->networkLinks())
                //                    {
                //                        link.second->addToStressStraight(straightSegmentsDeq);
                //                    }
                //
                //                    currentSize=straightSegmentsDeq.size();
                //
                //                    const VectorDim meshSize(this->mesh.xMax()-this->mesh.xMin());
                //
                //                    for(int i=-dislocationImages_x;i<=dislocationImages_x;++i)
                //                    {
                //                        for(int j=-dislocationImages_y;j<=dislocationImages_y;++j)
                //                        {
                //                            for(int k=-dislocationImages_z;k<=dislocationImages_z;++k)
                //                            {
                //
                //                                const Eigen::Matrix<int,3,1> cellID((Eigen::Matrix<int,3,1>()<<i,j,k).finished());
                //
                //                                if( cellID.squaredNorm()!=0) //skip current cell
                //                                {
                //                                    for (size_t c=0;c<currentSize;++c)
                //                                    {
                //                                        const VectorDim P0=straightSegmentsDeq[c].P0+(meshSize.array()*cellID.cast<double>().array()).matrix();
                //                                        const VectorDim P1=straightSegmentsDeq[c].P1+(meshSize.array()*cellID.cast<double>().array()).matrix();
                //
                //                                        straightSegmentsDeq.emplace_back(P0,P1,straightSegmentsDeq[c].b);
                //                                    }
                //                                }
                //                            }
                //                        }
                //                    }
                //
                //
                //                }
                //                model::cout<< straightSegmentsDeq.size()<<" straight segments ("<<currentSize<<"+"<<straightSegmentsDeq.size()-currentSize<<" images)"<<std::flush;
                //                model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<defaultColor<<std::endl;
                //
                const auto t1= std::chrono::system_clock::now();
                model::cout<<"		Computing analytical stress field at quadrature points ("<<nThreads<<" threads) "<<std::flush;
#ifdef _OPENMP
                //                const size_t nThreads = omp_get_max_threads();
                EqualIteratorRange<typename NetworkLinkContainerType::iterator> eir(this->links().begin(),this->links().end(),nThreads);
                //                SOMETHING WRONG HERE? CPU USE SAYS 100% ONLY?
                
#pragma omp parallel for
                for(size_t thread=0;thread<eir.size();thread++)
                {
                    for (auto& linkIter=eir[thread].first;linkIter!=eir[thread].second;linkIter++)
                    {
                        linkIter->second->assembleGlide();
                    }
                }
#else
                for (auto& linkIter : this->links())
                {
                    linkIter.second->assembleGlide();
                }
#endif
                model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t1)).count()<<" sec]."<<defaultColor<<std::endl;
                
                
                
                
                
                
                //#ifdef _OPENMP
                //#pragma omp parallel for
                //#endif
                //                for (unsigned int k=0; k<this->particleSystem().size();++k)
                //                {
                //                    if(this->particleSystem()[k].template fieldPointBase<StressField>().enabled)
                //                    {
                //                        this->particleSystem()[k].template fieldPointBase<StressField>() += StressField::addSourceContribution(this->particleSystem()[k],straightSegmentsDeq);
                //                    }
                //                }
                //                model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<defaultColor<<std::endl;
            }
            else
            {// For curved segments use quandrature integration of stress field
                //                assert(0 && "RE-ENABLE THIS. New Material class is incompatible with old implmentation using static objects");
                
                assert(0 && "ALL THIS MUST BE RE-IMPLEMENTED FOR CURVED SEGMENTS");
                
                //
                //                const auto t0= std::chrono::system_clock::now();
                //
                //                if(computeDDinteractions)
                //                {
                //
                //                    if(dislocationImages_x!=0 || dislocationImages_y!=0 || dislocationImages_z!=0)
                //                    {
                //                        assert(0 && "FINISH HERE");
                //                    }
                //
                //                    model::cout<<"		Computing numerical stress field at quadrature points ("<<nThreads<<" threads)..."<<std::flush;
                //                    if (use_extraStraightSegments)
                //                    {
                //                        this->template computeNeighborField<StressField>(ssdeq);
                //                    }
                //                    else
                //                    {
                //                        this->template computeNeighborField<StressField>();
                //                    }
                //                    model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<defaultColor<<std::endl;
                //                }
                //
                //                //! -2 Loop over DislocationSegments and assemble stiffness matrix and force vector
                //                const auto t1= std::chrono::system_clock::now();
                //                model::cout<<"		Computing segment stiffness matrices and force vectors ("<<nThreads<<" threads)..."<<std::flush;
                //                typedef void (LinkType::*LinkMemberFunctionPointerType)(void); // define type of Link member function
                //                assert(0 && "THIS CASE MUST BE REWORKED, since LinkType::assemble DOES NOT EXIST ANYMORE");
                //                //LinkMemberFunctionPointerType Lmfp(&LinkType::assemble); // Lmfp is a member function pointer to Link::assemble
                //                //this->parallelExecute(Lmfp);
                //                model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t1)).count()<<" sec]."<<defaultColor<<std::endl;
                //
                
            }
            
            //! -3 Loop over DislocationSubNetworks, assemble subnetwork stiffness matrix and force vector, and solve
            const auto t2= std::chrono::system_clock::now();
            model::cout<<"		Assembling NetworkComponents and solving "<<std::flush;
            
            
            
            switch (ddSolverType)
            {
                case 1: // iterative solver
                {
                    model::cout<<"(MINRES "<<nThreads<<" threads)..."<<std::flush;
#ifdef _OPENMP // MINRES is not multi-threaded. So parallelize loop over NetworkComponents.
#pragma omp parallel for
                    for (unsigned int k=0;k<this->components().size();++k)
                    {
                        auto snIter(this->components().begin());
                        std::advance(snIter,k);
                        DislocationNetworkComponentType(*snIter->second).iterativeSolve();
                    }
#else
                    for (const auto& networkComponent : this->components())
                    {
                        DislocationNetworkComponentType(*networkComponent.second).iterativeSolve();
                    }
#endif
                    break;
                }
                    
                case 2: // direct solver
                {
#ifdef _MODEL_PARDISO_SOLVER_
                    model::cout<<"(PardisoLDLT "<<nThreads<<" threads)..."<<std::flush;
                    for (const auto& networkComponent : this->components())
                    {
                        DislocationNetworkComponentType(*networkComponent.second).directSolve();
                    }
#else
                    model::cout<<"(SimplicialLDLT "<<nThreads<<" threads)..."<<std::flush;
#ifdef _OPENMP // SimplicialLDLT is not multi-threaded. So parallelize loop over NetworkComponents.
#pragma omp parallel for
                    for (unsigned int k=0;k<this->components().size();++k)
                    {
                        auto snIter(this->components().begin());
                        std::advance(snIter,k);
                        DislocationNetworkComponentType(*snIter->second).directSolve();
                    }
#else
                    for (const auto& networkComponent : this->components())
                    {
                        DislocationNetworkComponentType(*networkComponent.second).directSolve();
                    }
#endif
#endif
                    break;
                }
                    
                default: // lumped solver
                {
                    model::cout<<"(lumpedSolver "<<nThreads<<" threads)..."<<std::flush;
#ifdef _OPENMP
#pragma omp parallel for
                    for (unsigned int k=0;k<this->components().size();++k)
                    {
                        auto snIter(this->components().begin());
                        std::advance(snIter,k);
                        DislocationNetworkComponentType(*snIter->second).lumpedSolve(runID);
                    }
#else
                    for (const auto& networkComponent : this->components())
                    {
                        DislocationNetworkComponentType(*networkComponent.second).lumpedSolve(runID);
                    }
#endif
                    break;
                }
            }
            
            //            if(DislocationNetworkComponentType::use_directSolver)
            //            {
            //
            //
            //            }
            //            else // iterative solver
            //            {
            //
            //            }
            
            model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t2)).count()<<" sec]."<<defaultColor<<std::endl;
        }
        
        
        
        
        /**********************************************************************/
        void moveGlide(const double & dt_in)
        {/*! Moves all nodes in the DislocationNetwork using the stored velocity and current dt
          */
            model::cout<<"		Moving DislocationNodes (dt="<<dt_in<< ")... "<<std::flush;
            const auto t0= std::chrono::system_clock::now();
            for (auto& nodeIter : this->nodes())
            {
                nodeIter.second->moveGlide(dt_in);
            }
            model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<defaultColor<<std::endl;
        }
        
        
        
        /**********************************************************************/
        void updatePlasticDistortionFromAreas(const double& dt)
        {
            //            if(!computePlasticDistortionRateFromVelocities)
            //            {
            const MatrixDimD old(_plasticDistortionFromAreas);
            _plasticDistortionFromAreas.setZero();
            for(const auto& loop : this->loops())
            {
                _plasticDistortionFromAreas+= loop.second->plasticDistortion();
            }
            _plasticDistortionRateFromAreas=(_plasticDistortionFromAreas-old)/dt;
            //            }
        }
        
        //        /**********************************************************************/
        //        const MatrixDimD& plasticDistortionRate() const
        //        {
        //            return computePlasticDistortionRateFromVelocities? _plasticDistortionRateFromVelocities : _plasticDistortionRateFromAreas;
        //        }
        
        /**********************************************************************/
        const MatrixDimD& plasticDistortionRate() const
        {
            return  _plasticDistortionRateFromAreas;
        }
        
        /**********************************************************************/
        const MatrixDimD& plasticDistortion() const
        {
            return _plasticDistortionFromAreas;
        }
        
        //        /**********************************************************************/
        //        const MatrixDimD& plasticDistortion() const
        //        {
        //            return computePlasticDistortionRateFromVelocities? _plasticDistortionFromVelocities : _plasticDistortionFromAreas;
        //        }
        
        /**********************************************************************/
        MatrixDimD plasticStrainRate() const
        {/*!\returns the plastic strain rate tensor generated during the last time step.
          */
            //const MatrixDimD temp(plasticDistortionRate());
            return (plasticDistortionRate()+plasticDistortionRate().transpose())*0.5;
        }
        
        /**********************************************************************/
        std::tuple<double,double,double,double> networkLength() const
        {/*!\returns the total line length of the DislocationNetwork. The return
          * value is a tuple, where the first value is the length of bulk glissile
          * dislocations, the second value is the length of bulk sessile
          * dislocations, and the third value is the length accumulated on
          * the mesh boundary.
          */
            double bulkGlissileLength(0.0);
            double bulkSessileLength(0.0);
            double boundaryLength(0.0);
            double grainBoundaryLength(0.0);
            
            for(auto& loop : this->loops())
            {
                for(const auto& loopLink : loop.second->links())
                {
                    if(!loopLink.second->pLink->hasZeroBurgers())
                    {
                        if(loopLink.second->pLink->isBoundarySegment())
                        {
                            boundaryLength+=loopLink.second->pLink->chord().norm();
                        }
                        else if(loopLink.second->pLink->isGrainBoundarySegment())
                        {
                            grainBoundaryLength+=loopLink.second->pLink->chord().norm();
                        }
                        else
                        {
                            if(loopLink.second->pLink->isSessile())
                            {
                                bulkSessileLength+=loopLink.second->pLink->chord().norm()/loopLink.second->pLink->loopLinks().size();
                            }
                            else
                            {
                                bulkGlissileLength+=loopLink.second->pLink->chord().norm()/loopLink.second->pLink->loopLinks().size();
                            }
                        }
                    }
                }
            }
            return std::make_tuple(bulkGlissileLength,bulkSessileLength,boundaryLength,grainBoundaryLength);
        }
        
        /**********************************************************************/
        MatrixDimD stress(const VectorDim& x) const
        {/*!\param[in] P position vector
          * \returns The stress field generated by the DislocationNetwork at P
          *
          * Note:
          */
            MatrixDimD temp(MatrixDimD::Zero());
            for(const auto& link : this->links())
            {// sum stress field per segment
                if(   !link.second->hasZeroBurgers()
                   && !link.second->isBoundarySegment()
                   && !(link.second->isVirtualBoundarySegment() && simulationParameters.simulationType==DefectiveCrystalParameters::PERIODIC_IMAGES)
                   )
                {
                    for(const auto& shift : periodicShifts)
                    {
                        temp+=link.second->straight.stress(x+shift);
                    }
                }
            }
            return temp;
        }
        
        /**********************************************************************/
        void stress(std::deque<FEMfaceEvaluation<ElementType,dim,dim>,Eigen::aligned_allocator<FEMfaceEvaluation<ElementType,dim,dim>>>& fieldPoints) const
        {
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for(size_t k=0;k<fieldPoints.size();++k)
            {
                fieldPoints[k]=stress(fieldPoints[k].P);
            }
        }
        
        /**********************************************************************/
        VectorDim displacement(const VectorDim& x) const
        {/*!\param[in] P position vector
          * \returns The stress field generated by the DislocationNetwork at P
          *
          * Note:
          */
            
            VectorDim temp(VectorDim::Zero());
            
            for(const auto& loop : this->loops())
            {// sum solid angle of each loop
                if(!(loop.second->isVirtualBoundaryLoop() && simulationParameters.simulationType==DefectiveCrystalParameters::PERIODIC_IMAGES))
                {
                    for(const auto& shift : periodicShifts)
                    {
                        temp-=loop.second->solidAngle(x+shift)/4.0/M_PI*loop.second->burgers();
                    }
                }
            }
            
            for(const auto& link : this->links())
            {// sum line-integral part of displacement field per segment
                if(   !link.second->hasZeroBurgers()
                   && (!link.second->isBoundarySegment() || simulationParameters.simulationType==DefectiveCrystalParameters::FINITE_NO_FEM)
                   && !(link.second->isVirtualBoundarySegment() && simulationParameters.simulationType==DefectiveCrystalParameters::PERIODIC_IMAGES)
                   )
                {
                    for(const auto& shift : periodicShifts)
                    {
                        temp+=link.second->straight.displacement(x+shift);
                    }
                }
            }
            
            return temp;
        }
        
        /**********************************************************************/
        void displacement(std::vector<FEMnodeEvaluation<ElementType,dim,1>,Eigen::aligned_allocator<FEMnodeEvaluation<ElementType,dim,1>>>& fieldPoints) const
        {
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for(size_t k=0;k<fieldPoints.size();++k)
            {
                fieldPoints[k]=displacement(fieldPoints[k].P);
            }
        }
        
        
        /**********************************************************************/
        DislocationNetworkIOType io()
        {
            return DislocationNetworkIOType(*this,folderSuffix);
        }
        
        /**********************************************************************/
        DislocationNetworkIOType io() const
        {
            return DislocationNetworkIOType(*this,folderSuffix);
        }
        
    };
    
}
#endif

//        /**********************************************************************/
//        double vacancyConcentration(const VectorDim& x) const
//        {
//            double temp(0.0);
//            for(const auto& link : this->links())
//            {// sum line-integral part of displacement field per segment
//                if(   !link.second->hasZeroBurgers()
//                   && !link.second->isBoundarySegment()
//                   && !link.second->isGrainBoundarySegment()
//                   )
//                {
//                    for(const auto& shift : periodicShifts)
//                    {
//                        temp+=link.second->vacancyConcentration(x+shift);
//                    }
//                }
//            }
//            return temp;
//        }
//
//        /**********************************************************************/
//        void vacancyConcentrationAssembly() const
//        {
//            std::vector<Eigen::Triplet<double>> lhsT;
//            std::vector<Eigen::Triplet<double>> rhsT;
//            size_t globalRow(0);
//            for(const auto& fieldLink : this->links())
//            {// sum line-integral part of displacement field per segment
//                if(   !fieldLink.second->hasZeroBurgers()
//                   && !fieldLink.second->isBoundarySegment()
//                   && !fieldLink.second->isGrainBoundarySegment())
//                {
//                    for(const auto& qPoint : fieldLink.second->quadraturePoints())
//                    {
//                        for(const auto& sourceLink : this->links())
//                        {// sum line-integral part of displacement field per segment
//                            if(   !sourceLink.second->hasZeroBurgers()
//                               && !sourceLink.second->isBoundarySegment()
//                               && !sourceLink.second->isGrainBoundarySegment())
//                            {
//                                for(const auto& shift : periodicShifts)
//                                {
//                                    sourceLink.second->addToVacancyConcentrationAssembly(qPoint.r+shift,
//                                                                                   lhsT,
//                                                                                   globalRow);
//                                }
//                            }
//                        }
////                        rhsT.emplace_back(globalRow,1,qPoint.equilibriumVacancyConcentration(*fieldLink.second)-bvpSolver->vacancyConcentration(qPoint.r,fieldLink.second->source->includingSimplex()));
//                        rhsT.emplace_back(globalRow,0,qPoint.equilibriumVacancyConcentration(*fieldLink.second)-1.0);
//                        globalRow++;
//                    }
//                }
//            }
//
////            for(const auto& trip : lhsT)
////            {
////                std::cout<<trip.row()<<","<<trip.col()<<" "<<trip.value()<<std::endl;
////            }
//            Eigen::SparseMatrix<double> K(globalRow,this->nodes().size()*dim);
//            K.setFromTriplets(lhsT.begin(),lhsT.end());
//            std::ofstream kFile("K.txt");
//            kFile<<K.toDense()<<std::endl;
//
////            for(const auto& trip : rhsT)
////            {
////                std::cout<<trip.row()<<","<<trip.col()<<" "<<trip.value()<<std::endl;
////            }
//            Eigen::SparseMatrix<double> F(globalRow,1);
//            F.setFromTriplets(rhsT.begin(),rhsT.end());
//            std::ofstream fFile("F.txt");
//            fFile<<F.toDense()<<std::endl;
//
//
//
//        }
//
//        /**********************************************************************/
//        void vacancyConcentration(std::vector<FEMnodeEvaluation<ElementType,1,1>,Eigen::aligned_allocator<FEMnodeEvaluation<ElementType,1,1>>>& fieldPoints) const
//        {
//#ifdef _OPENMP
//#pragma omp parallel for
//#endif
//            for(size_t k=0;k<fieldPoints.size();++k)
//            {
//                fieldPoints[k](0,0)=vacancyConcentration(fieldPoints[k].P);
//            }
//        }



//        /**********************************************************************/
//        void runSteps() __attribute__ ((deprecated))
//        {/*! Runs Nsteps simulation steps
//          */
//            const auto t0= std::chrono::system_clock::now();
//            while (runID<Nsteps)
//            {
//                model::cout<<std::endl; // leave a blank line
//                singleStep();
//            }
//            updateQuadraturePoints(); // necessary if quadrature data are necessary in main
//            model::cout<<greenBoldColor<<std::setprecision(3)<<std::scientific<<Nsteps<< " simulation steps completed in "<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" [sec]"<<defaultColor<<std::endl;
//        }

//        /**********************************************************************/
//        void updatePlasticDistortionRateFromVelocities()
//        {
//            if(computePlasticDistortionRateFromVelocities)
//            {
//                _plasticDistortionRateFromVelocities.setZero();
//                for (const auto& linkIter : this->links())
//                {
//                    _plasticDistortionRateFromVelocities+= linkIter.second->plasticDistortionRate();
//                }
//                _plasticDistortionFromVelocities += _plasticDistortionRateFromVelocities*dt;
//            }
//        }

//        /**********************************************************************/
//        void updateLoadControllers(const long int& runID)
//        {/*! Updates bvpSolver using the stress and displacement fields of the
//          *  current DD configuration.
//          */
//            const int quadraturePerTriangle=37;
//            if(use_bvp)
//            {
//                if (!(runID%use_bvp))
//                {// enter the if statement if use_bvp!=0 and runID is a multiple of use_bvp
//                    model::cout<<"		Updating elastic bvp... "<<std::endl;
//                    this->bvpSolver.template assembleAndSolve<DislocationNetworkType,quadraturePerTriangle>(*this);
//                }
//            }
//            if (use_externalStress)
//            {
//                externalLoadController.update(*this,runID);
//            }
//        }

//        /**********************************************************************/
//        const unsigned int& userOutputColumn()  const
//        {
//            return _userOutputColumn;
//        }

//        /**********************************************************************/
//        std::pair<bool,const Simplex<dim,dim>*> pointIsInsideMesh(const VectorDim& P0, const Simplex<dim,dim>* const guess) const
//        {/*!\param[in] P0 position vector
//          * \param[in] guess pointer of the Simplex where the search starts
//          * \returns true if P0 is inside the mesh
//          */
////            std::pair<bool,const Simplex<dim,dim>*> temp(true,NULL);
////            if (use_boundary)
////            {
////                temp=this->mesh.searchWithGuess(P0,guess);
////            }
//            return this->mesh.searchWithGuess(P0,guess);
//        }

//        /**********************************************************************/
//        const ParticleSystemType& particleSystem() const
//        {
//            return *this;
//        }
//
//        /**********************************************************************/
//        ParticleSystemType& particleSystem()
//        {
//            return *this;
//        }

//        /**********************************************************************/
//        void updateStressStraightSegments()
//        {
//
//            const auto t0= std::chrono::system_clock::now();
//            model::cout<<"		Collecting StressStraight objects: "<<std::flush;
//
//            straightSegmentsDeq.clear();
//            //                std::deque<StressStraight<dim>,Eigen::aligned_allocator<StressStraight<dim>>> straightSegmentsDeq;
//            size_t currentSize=0;
////            if(computeDDinteractions)
////            {
//                for(const auto& link : this->networkLinks())
//                {
//                    link.second->addToStressStraight(straightSegmentsDeq);
//                }
//
//                currentSize=straightSegmentsDeq.size();
//
//                const VectorDim meshSize(this->mesh.xMax()-this->mesh.xMin());
//
//                for(int i=-dislocationImages_x;i<=dislocationImages_x;++i)
//                {
//                    for(int j=-dislocationImages_y;j<=dislocationImages_y;++j)
//                    {
//                        for(int k=-dislocationImages_z;k<=dislocationImages_z;++k)
//                        {
//
//                            const Eigen::Matrix<int,3,1> cellID((Eigen::Matrix<int,3,1>()<<i,j,k).finished());
//
//                            if( cellID.squaredNorm()!=0) //skip current cell
//                            {
//                                for (size_t c=0;c<currentSize;++c)
//                                {
//                                    const VectorDim P0=straightSegmentsDeq[c].P0+(meshSize.array()*cellID.cast<double>().array()).matrix();
//                                    const VectorDim P1=straightSegmentsDeq[c].P1+(meshSize.array()*cellID.cast<double>().array()).matrix();
//
//                                    straightSegmentsDeq.emplace_back(P0,P1,straightSegmentsDeq[c].b);
//                                }
//                            }
//                        }
//                    }
//                }
//
//
////            }
//            model::cout<< straightSegmentsDeq.size()<<" straight segments ("<<currentSize<<"+"<<straightSegmentsDeq.size()-currentSize<<" images)"<<std::flush;
//            model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<defaultColor<<std::endl;
//
//        }

//        /**********************************************************************/
//        const double& get_dt() const
//        {/*!\returns the current time step increment dt
//          */
//            return dt;
//        }
//
//        /**********************************************************************/
//        void set_dt(const double& dt_in,const double& vMax_in)
//        {/*!\param[in] dt_in
//          * Sets the time increment dt to dt_in
//          */
//            dt=dt_in;
//            vMax=vMax_in;
//        }

//        /**********************************************************************/
//        const double& get_totalTime() const
//        {/*! The elapsed simulation time step in dimensionless units
//          */
//            return totalTime;
//        }

//        /**********************************************************************/
//        void updateQuadraturePoints()
//        {
//            if(corder>0)
//            {// quadrature points only used for curved segments
//                model::cout<<"        Updating quadrature points... "<<std::flush;
//                const auto t0=std::chrono::system_clock::now();
//
//                // Clear DislocationParticles
//                this->clearParticles(); // this also destroys all cells
//
//                // Populate DislocationParticles
//                for(auto& linkIter : this->links())
//                {
//                    linkIter.second->updateQuadraturePoints(*this);
//                }
//                model::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
//            }
//        }
