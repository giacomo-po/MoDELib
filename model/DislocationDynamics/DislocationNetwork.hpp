/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationNetwork_cpp_
#define model_DislocationNetwork_cpp_


#include <DislocationNetwork.h>


namespace model
{
    
    /**********************************************************************/
    template <int dim, short unsigned int corder, typename InterpolationType>
    DislocationNetwork<dim,corder,InterpolationType>::DislocationNetwork(int& argc, char* argv[],
                                                                         const DefectiveCrystalParameters& _simulationParameters,
                                                                         const SimplicialMesh<dim>& _mesh,
                                                                         const Polycrystal<dim>& _poly,
                                                                         const std::unique_ptr<BVPsolver<dim,2>>& _bvpSolver,
                                                                         const std::unique_ptr<ExternalLoadControllerBase<dim>>& _externalLoadController,
                                                                         const std::vector<VectorDim>& _periodicShifts,
                                                                         long int& runID) :
    //    /* init */ LoopNetworkType(_simulationParameters.isPeriodicSimulation()? std::shared_ptr<NetworkComponentType>(new NetworkComponentType()) : nullptr)
    /* init */ simulationParameters(_simulationParameters)
    /* init */,mesh(_mesh)
    /* init */,poly(_poly)
    /* init */,glidePlaneFactory(poly)
    /* init */,periodicGlidePlaneFactory(simulationParameters.isPeriodicSimulation()? new PeriodicGlidePlaneFactory<dim>(poly, glidePlaneFactory) : nullptr)
    //    /* init */,periodicDislocationLoopFactory(simulationParameters.isPeriodicSimulation()? new PeriodicDislocationLoopFactory<DislocationNetworkType>(poly,glidePlaneFactory) : nullptr)
    /* init */,bvpSolver(_bvpSolver)
    /* init */,externalLoadController(_externalLoadController)
    /* init */,periodicShifts(_periodicShifts)
        /* init */,networkRemesher(*this)
        /* init */,junctionsMaker(*this)
        /* init */,nodeContractor(*this)
    //    /* init */,gbTransmission(*this)
    //        /* init */,timeIntegrationMethod(TextFileParser("inputFiles/DD.txt").readScalar<int>("timeIntegrationMethod",true))
    ///* init */,maxJunctionIterations(TextFileParser("inputFiles/DD.txt").readScalar<int>("maxJunctionIterations",true))
    //        /* init */,runID(TextFileParser("inputFiles/DD.txt").readScalar<int>("startAtTimeStep",true)),
    //        /* init */,totalTime(0.0),
    //        /* init */ dt(0.0),
    //        /* init */ vMax(0.0),
    //        /* init */ Nsteps(TextFileParser("inputFiles/DD.txt").readScalar<size_t>("Nsteps",true)),
    //        /* init */,_plasticDistortionFromVelocities(MatrixDim::Zero())
    /* init */,oldPlasticDistortionFromAreas(std::make_pair(0.0,MatrixDim::Zero()))
    /* init */,_plasticDistortionRateFromVelocities(MatrixDim::Zero())
    /* init */,_plasticDistortionRateFromAreas(MatrixDim::Zero())
    /* init */,ddSolverType(TextFileParser("inputFiles/DD.txt").readScalar<int>("ddSolverType",true))
    /* init */,computeDDinteractions(TextFileParser("inputFiles/DD.txt").readScalar<int>("computeDDinteractions",true))
    /* init */,crossSlipModel(TextFileParser("inputFiles/DD.txt").readScalar<int>("crossSlipModel",true))
    /* init */,outputFrequency(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputFrequency",true))
    /* init */,outputBinary(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputBinary",true))
    /* init */,outputGlidePlanes(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputGlidePlanes",true))
    /* init */,outputElasticEnergy(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputElasticEnergy",true))
    /* init */,outputMeshDisplacement(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputMeshDisplacement",true))
    /* init */,outputFEMsolution(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputFEMsolution",true))
    /* init */,outputDislocationLength(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputDislocationLength",true))
    //        /* init */,outputPlasticDistortion(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputPlasticDistortion",true))
    /* init */,outputPlasticDistortionRate(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputPlasticDistortionRate",true))
    /* init */,outputQuadraturePoints(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputQuadraturePoints",true))
    /* init */,outputLinkingNumbers(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputLinkingNumbers",true))
    /* init */,outputLoopLength(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputLoopLength",true))
    /* init */,outputSegmentPairDistances(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputSegmentPairDistances",true))
    /* init */,computeElasticEnergyPerLength(TextFileParser("inputFiles/DD.txt").readScalar<int>("computeElasticEnergyPerLength",true))
    //    /* init */,outputPeriodicConfiguration(simulationParameters.isPeriodicSimulation()? TextFileParser("inputFiles/DD.txt").readScalar<int>("outputPeriodicConfiguration",true) : false)
    //        /* init */ _userOutputColumn(3)
    /* init */,use_stochasticForce(TextFileParser("inputFiles/DD.txt").readScalar<int>("use_stochasticForce",true))
    /* init */,surfaceAttractionDistance(TextFileParser("inputFiles/DD.txt").readScalar<double>("surfaceAttractionDistance",true))
    //        /* init */,computePlasticDistortionRateFromVelocities(TextFileParser("inputFiles/DD.txt").readScalar<int>("computePlasticDistortionRateFromVelocities",true))
    /* init */,folderSuffix("")
    {
        
        // Some sanity checks
        //            assert(Nsteps>=0 && "Nsteps MUST BE >= 0");
        
        // Initialize static variables
        LoopNetworkType::verboseLevel=TextFileParser("inputFiles/DD.txt").readScalar<int>("verboseLoopNetwork",true);
        verboseDislocationNetwork=TextFileParser("inputFiles/DD.txt").readScalar<int>("verboseDislocationNetwork",true);
        LoopType::initFromFile("inputFiles/DD.txt");
        LoopNodeType::initFromFile("inputFiles/DD.txt");
        LoopLinkType::initFromFile("inputFiles/DD.txt");
        NetworkLinkType::initFromFile("inputFiles/DD.txt");
        NetworkNodeType::initFromFile("inputFiles/DD.txt");
        //        PeriodicDislocationBase::initFromFile("inputFiles/DD.txt");
        //        DislocationNetworkComponentType::initFromFile("inputFiles/DD.txt");
        DislocationStressBase<dim>::initFromFile("inputFiles/DD.txt");
        //        DDtimeIntegrator<0>::initFromFile("inputFiles/DD.txt");
        //        DislocationCrossSlip<DislocationNetworkType>::initFromFile("inputFiles/DD.txt");
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
        
        DDconfigIO<dim> evl(folderSuffix);
        evl.read(runID);
        setConfiguration(evl);
        createEshelbyInclusions();
    }
    
    /**********************************************************************/
    template <int dim, short unsigned int corder, typename InterpolationType>
    void DislocationNetwork<dim,corder,InterpolationType>::setConfiguration(const DDconfigIO<dim>& evl)
    {
        this->loopLinks().clear(); // erase base network to clear current config
        
        // Create Loops
        std::deque<std::shared_ptr<LoopType>> tempLoops; // keep loops alive during setConfiguration
        size_t loopNumber=1;
        for(const auto& loop : evl.loops())
        {
            const bool faulted= poly.grain(loop.grainID).rationalLatticeDirection(loop.B).rat.asDouble()!=1.0? true : false;
            model::cout<<"Creating DislocationLoop "<<loop.sID<<" ("<<loopNumber<<" of "<<evl.loops().size()<<"), type="<<loop.loopType<<", faulted="<<faulted<<", |b|="<<loop.B.norm()<<std::endl;
            const size_t loopIDinFile(loop.sID);
            LoopType::set_count(loopIDinFile);
            GlidePlaneKey<dim> loopPlaneKey(loop.P,poly.grain(loop.grainID).reciprocalLatticeDirection(loop.N));
            tempLoops.push_back(this->loops().create(loop.B,glidePlaneFactory.getFromKey(loopPlaneKey)));
            assert(this->loops().get(loopIDinFile)->sID==loopIDinFile);
            loopNumber++;
        }
        
        // Create NetworkNodes
        std::deque<std::shared_ptr<NetworkNodeType>> tempNetNodes; // keep loops alive during setConfiguration
        size_t netNodeNumber=1;
        for(const auto& node : evl.nodes())
        {
            model::cout<<"Creating DislocationNode "<<node.sID<<" ("<<netNodeNumber<<" of "<<evl.nodes().size()<<")"<<std::endl;
            const size_t nodeIDinFile(node.sID);
            NetworkNodeType::set_count(nodeIDinFile);
            tempNetNodes.push_back(this->networkNodes().create(node.P,node.V,node.velocityReduction));
            assert(this->networkNodes().get(nodeIDinFile)->sID==nodeIDinFile);
            netNodeNumber++;
        }
        
        // Create LoopNodes
        std::deque<std::shared_ptr<LoopNodeType>> tempLoopNodes; // keep loops alive during setConfiguration
        size_t loopNodeNumber=1;
        for(const auto& node : evl.loopNodes())
        {
            model::cout<<"Creating DislocationLoopNode "<<node.sID<<" ("<<loopNodeNumber<<" of "<<evl.loopNodes().size()<<")"<<std::flush;
            const size_t nodeIDinFile(node.sID);
            LoopNodeType::set_count(nodeIDinFile);
            const auto loop(this->loops().get(node.loopID));
            const auto netNode(this->networkNodes().get(node.networkNodeID));
            assert(loop && "Loop does not exist");
            assert(netNode && "NetworkNode does not exist");
            const auto periodicPatch(loop->periodicGlidePlane? loop->periodicGlidePlane->patches().getFromKey(node.periodicShift) : nullptr);
            const auto periodicPatchEdge((periodicPatch && node.edgeID>=0)? periodicPatch->edges()[node.edgeID] : nullptr);
            if(periodicPatch)
            {
                model::cout<<", on patch "<<periodicPatch->shift.transpose()<<std::flush;
            }
            if(periodicPatchEdge)
            {
                model::cout<<" on edge "<<periodicPatchEdge->edgeID<<std::flush;
            }
            model::cout<<std::endl;
            tempLoopNodes.push_back(this->loopNodes().create(loop,netNode,node.P,periodicPatch,periodicPatchEdge));
            assert(this->loopNodes().get(nodeIDinFile)->sID==nodeIDinFile);
            loopNodeNumber++;
        }
        
        // Insert Loops
        std::map<size_t,std::map<size_t,size_t>> loopMap;
        for(const auto& looplink : evl.loopLinks())
        {// Collect LoopLinks by loop IDs
            loopMap[looplink.loopID].emplace(looplink.sourceID,looplink.sinkID);
        }
        assert(loopMap.size()==evl.loops().size());
        
        for(const auto& loop : evl.loops())
        {// for each loop in the DDconfigIO<dim> object
            
            const auto loopFound=loopMap.find(loop.sID); // there must be an entry with key loopID in loopMap
            assert(loopFound!=loopMap.end());
            std::vector<std::shared_ptr<LoopNodeType>> loopNodes;
            loopNodes.push_back(this->loopNodes().get(loopFound->second.begin()->first));
            for(size_t k=0;k<loopFound->second.size();++k)
            {
                const auto nodeFound=loopFound->second.find(loopNodes.back()->sID);
                if(k<loopFound->second.size()-1)
                {
                    loopNodes.push_back(this->loopNodes().get(nodeFound->second));
                }
                else
                {
                    assert(nodeFound->second==loopNodes[0]->sID);
                }
            }
            this->insertLoop(this->loops().get(loop.sID),loopNodes);
        }
        updateGeometry();
    }
    
    /**********************************************************************/
    template <int dim, short unsigned int corder, typename InterpolationType>
    const typename DislocationNetwork<dim,corder,InterpolationType>::EshelbyInclusionContainerType& DislocationNetwork<dim,corder,InterpolationType>::eshelbyInclusions() const
    {
        return *this;
    }
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    typename DislocationNetwork<dim,corder,InterpolationType>::EshelbyInclusionContainerType& DislocationNetwork<dim,corder,InterpolationType>::eshelbyInclusions()
    {
        return *this;
    }
    
    
    /**********************************************************************/
    template <int dim, short unsigned int corder, typename InterpolationType>
    void DislocationNetwork<dim,corder,InterpolationType>::createEshelbyInclusions()
    {
        for(const auto& grain : poly.grains())
        {
            EshelbyInclusion<dim>::addSlipSystems(grain.second.slipSystems());
        }
        
        
        IDreader<'E',1,14,double> inclusionsReader;
        inclusionsReader.read(0,true);
        
        const std::vector<double> inclusionsMobilityReduction(TextFileParser("./inputFiles/initialMicrostructure.txt").readArray<double>("inclusionsMobilityReduction",true));
        for(const auto& pair : inclusionsReader)
        {
            
            const size_t& inclusionID(pair.first);
            Eigen::Map<const Eigen::Matrix<double,1,14>> row(pair.second.data());
            
            const VectorDim C(row.template segment<dim>(0));
            const double a(row(dim+0));
            MatrixDim eT(MatrixDim::Zero());
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
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    void DislocationNetwork<dim,corder,InterpolationType>::updateGeometry()
    {
        VerboseDislocationNetwork(2,"DislocationNetwork::updateGeometry"<<std::endl;);
        for(auto& loop : this->loops())
        {// copmute slipped areas and right-handed normal // TODO: PARALLELIZE THIS LOOP
            loop.second.lock()->updateGeometry();
        }
        updatePlasticDistortionRateFromAreas();
        VerboseDislocationNetwork(3,"DislocationNetwork::updateGeometry DONE"<<std::endl;);
    }
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    void DislocationNetwork<dim,corder,InterpolationType>::updatePlasticDistortionRateFromAreas()
    {
        VerboseDislocationNetwork(2,"DislocationNetwork::updatePlasticDistortionRateFromAreas"<<std::endl;);
        const double dt(simulationParameters.totalTime-oldPlasticDistortionFromAreas.first);
        if(dt>=0.0)
        {
            const MatrixDim pd(plasticDistortion());
            if(dt>0.0)
            {
                _plasticDistortionRateFromAreas=(pd-oldPlasticDistortionFromAreas.second)/dt;
            }
            oldPlasticDistortionFromAreas=std::make_pair(simulationParameters.totalTime,pd); // update old values
        }
        VerboseDislocationNetwork(3,"DislocationNetwork::updatePlasticDistortionRateFromAreas DONE"<<std::endl;);
    }
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    typename DislocationNetwork<dim,corder,InterpolationType>::DislocationNetworkIOType DislocationNetwork<dim,corder,InterpolationType>::io()
    {
        return DislocationNetworkIOType(*this,folderSuffix);
    }
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    typename DislocationNetwork<dim,corder,InterpolationType>::DislocationNetworkIOType DislocationNetwork<dim,corder,InterpolationType>::io() const
    {
        return DislocationNetworkIOType(*this,folderSuffix);
    }
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    typename DislocationNetwork<dim,corder,InterpolationType>::MatrixDim DislocationNetwork<dim,corder,InterpolationType>::plasticDistortion() const
    {
        MatrixDim temp(MatrixDim::Zero());
        for(const auto& loop : this->loops())
        {
            temp+= loop.second.lock()->plasticDistortion();
        }
        return temp;
    }
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    const typename DislocationNetwork<dim,corder,InterpolationType>::MatrixDim& DislocationNetwork<dim,corder,InterpolationType>::plasticDistortionRate() const
    {
        return  _plasticDistortionRateFromAreas;
    }
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    typename DislocationNetwork<dim,corder,InterpolationType>::MatrixDim DislocationNetwork<dim,corder,InterpolationType>::plasticStrainRate() const
    {/*!\returns the plastic strain rate tensor generated during the last time step.
      */
        return (plasticDistortionRate()+plasticDistortionRate().transpose())*0.5;
    }
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    std::tuple<double,double,double,double> DislocationNetwork<dim,corder,InterpolationType>::networkLength() const
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
            for(const auto& loopLink : loop.second.lock()->loopLinks())
            {
                if(loopLink->networkLink())
                {
                    if(!loopLink->networkLink()->hasZeroBurgers())
                    {
                        if(loopLink->networkLink()->isBoundarySegment())
                        {
                            boundaryLength+=loopLink->networkLink()->chord().norm();
                        }
                        else if(loopLink->networkLink()->isGrainBoundarySegment())
                        {
                            grainBoundaryLength+=loopLink->networkLink()->chord().norm();
                        }
                        else
                        {
                            if(loopLink->networkLink()->isSessile())
                            {
                                bulkSessileLength+=loopLink->networkLink()->chord().norm()/loopLink->networkLink()->loopLinks().size();
                            }
                            else
                            {
                                bulkGlissileLength+=loopLink->networkLink()->chord().norm()/loopLink->networkLink()->loopLinks().size();
                            }
                        }
                    }
                }
            }
        }
        return std::make_tuple(bulkGlissileLength,bulkSessileLength,boundaryLength,grainBoundaryLength);
    }

    
    template <int dim, short unsigned int corder, typename InterpolationType>
    bool DislocationNetwork<dim,corder,InterpolationType>::contract(std::shared_ptr<NetworkNodeType> nA,
                  std::shared_ptr<NetworkNodeType> nB)
    {
        return nodeContractor.contract(nA,nB);
    }

    
    template <int dim, short unsigned int corder, typename InterpolationType>
    void DislocationNetwork<dim,corder,InterpolationType>::dummyMove(const int& runID)
    {/*!\returns the plastic strain rate tensor generated during the last time step.
      */
        std::cout<<"Dummy move: runID="<<runID<<std::endl;
        
        danglingBoundaryLoopNodes.clear();
        
        //        std::cout<<"Computing loop centers"<<std::endl;
        std::map<size_t,std::pair<VectorDim,size_t>> loopCenters;
        for(const auto& node : this->loopNodes())
        {
            const auto sharedNode(node.second.lock());
            const auto sharedLoop(sharedNode->loop());
            auto loopIter(loopCenters.find(sharedLoop->sID));
            if(loopIter==loopCenters.end())
            {
                loopCenters.emplace(sharedLoop->sID,std::make_pair(sharedNode->get_P(),1));
            }
            else
            {
                loopIter->second.first+=sharedNode->get_P();
                loopIter->second.second++;
            }
        }
        
        for(auto& pair : loopCenters)
        {
            pair.second.first/=pair.second.second;
        }
        
        //        std::cout<<"Set_P loop"<<std::endl;
        
        
        
        const double dt=10;
        for(auto& node : this->networkNodes())
        {// Expansion

            VectorDim v(VectorDim::Zero());
            for(const auto& loopNode : node.second.lock()->loopNodes())
            {
                auto loopIter(loopCenters.find(loopNode->loop()->sID));
                v+=(loopNode->get_P()-loopIter->second.first).normalized()/node.second.lock()->loopNodes().size();
            }

            if(node.second.lock()->loopNodes().size()>1)
            {
                const auto dir((node.second.lock()->glidePlaneIntersections()->P1-node.second.lock()->glidePlaneIntersections()->P0).normalized());
                v=v.dot(dir)*dir;
            }

            node.second.lock()->trySet_P(node.second.lock()->get_P()+v*dt);

        }
        
//        const auto loop0(this->loops().begin());
//        const auto loop1(this->loops().rbegin());
//
//                const auto v=loop0->second.lock()->glidePlane->unitNormal.cross(loop1->second.lock()->glidePlane->unitNormal).normalized();
//                for(auto& node : this->networkNodes())
//                {// translation
//                    node.second.lock()->trySet_P(node.second.lock()->get_P()+v*dt);
//                }
        
        
        
        std::cout<<"Removing Nodes"<<std::endl;
        for(const auto& node : danglingBoundaryLoopNodes)
        {
            //            std::cout<<"removing node="<<node->sID<<std::endl;
            this->removeLoopNode(node->sID);
        }
        
        // std::vector<std::tuple<typename PeriodicGlidePlane<dim>::VectorLowerDim,typename PeriodicGlidePlane<dim>::VectorDim,short int,const T* const>>
        
        std::cout<<"Inserting new boundary nodes"<<std::endl;
        //        std::map<std::pair<std::shared_ptr<NetworkLinkType>,float>,std::shared_ptr<NetworkNodeType>> newNetworkNodesMap;
        std::map<std::tuple<const NetworkNodeType*,const NetworkNodeType*,std::set<const PlanarMeshFace<dim>*>,float>,std::shared_ptr<NetworkNodeType>> newNetworkNodesMap;
        
        for(const auto& weakLoop : this->loops())
        {
            const auto loop(weakLoop.second.lock());
            
            //            std::cout<<"internal loop nodes "<<std::endl;
            
            std::vector<std::pair<VectorLowerDim,const LoopNodeType* const>> loopNodesPos;
            std::map<std::tuple<const LoopNodeType*,const LoopNodeType*,const PeriodicPlaneEdge<dim>*>,const LoopNodeType*> bndNodesMap;
            for(const auto& loopLink : loop->linkSequence())
            {
                if(!loopLink->source->periodicPlaneEdge)
                {
                    loopNodesPos.emplace_back(loop->periodicGlidePlane->referencePlane->localPosition(loopLink->source->get_P()),loopLink->source.get());
                    //                    std::cout<<loopNodesPos.back().first.transpose()<<std::endl;
                    
                }
                else
                {
                    bndNodesMap.emplace(std::make_tuple(loopLink->source->periodicPrev(),loopLink->source->periodicNext(),loopLink->source->periodicPlaneEdge.get()),loopLink->source.get());
                }
            }
            
            
            const auto polyInt(loop->periodicGlidePlane->polygonPatchIntersection(loopNodesPos)); // 2dPos,shift,edgeID,LoopNodeType*
            std::map<const LoopNodeType*,size_t> polyIntMap;
            for(size_t p=0;p<polyInt.size();++p)
            {
                if(std::get<3>(polyInt[p]))
                {
                    polyIntMap.emplace(std::get<3>(polyInt[p]),p);
                }
            }
            
            
            for(size_t k=0;k<loopNodesPos.size();++k)
            {
                const auto periodicPrev(loopNodesPos[k].second);
                const auto periodicPrevNetwork(periodicPrev->networkNode.get());
                const auto polyIter(polyIntMap.find(periodicPrev));
                assert(polyIter!=polyIntMap.end());
                const size_t p(polyIter->second);
                
                const size_t k1(k<loopNodesPos.size()-1? k+1 : 0);
                const auto periodicNext(loopNodesPos[k1].second);
                const auto periodicNextNetwork(periodicNext->networkNode.get());
                const auto polyIter1(polyIntMap.find(periodicNext));
                assert(polyIter1!=polyIntMap.end());
                const size_t p1(polyIter1->second);
                
                const auto periodicNetworkSource(periodicPrevNetwork->sID<periodicNextNetwork->sID? periodicPrevNetwork : periodicNextNetwork);
                const auto periodicNetworkSink  (periodicPrevNetwork->sID<periodicNextNetwork->sID? periodicNextNetwork : periodicPrevNetwork);
                
                const LoopNodeType* currentSource(periodicPrev);
                //                std::cout<<"periodicPrev="<<periodicPrev->sID<<std::endl;
                //                std::cout<<"periodicNext="<<periodicNext->sID<<std::endl;
                size_t p2 = (p+1)%polyInt.size();
                while(p2<p1)
                {
                    const auto periodicPatch(loop->periodicGlidePlane->getPatch(std::get<1>(polyInt[p2])));
                    const auto periodicPatchEdge(periodicPatch->edges()[std::get<2>(polyInt[p2])]);
                    const auto bndIter(bndNodesMap.find(std::make_tuple(periodicPrev,periodicNext,periodicPatchEdge.get())));
                    
                    if(bndIter!=bndNodesMap.end())
                    {// exising bnd node found
                        currentSource=bndIter->second;
                        //                        std::cout<<"existing point="<<std::get<0>(polyInt[p2]).transpose()<<std::endl;
                        
                    }
                    else
                    {
                        const VectorDim    loopNodePos(loop->periodicGlidePlane->referencePlane->globalPosition(std::get<0>(polyInt[p2])));
                        const VectorDim networkNodePos(loopNodePos+std::get<1>(polyInt[p2]));
                        const auto currentLoopLink(currentSource->next.second);
                        const auto currentNetworkLink(currentLoopLink->networkLink());
                        
                        
                        //                        const float u((networkNodePos-currentNetworkLink->source->get_P()).dot(currentNetworkLink->chord())/currentNetworkLink->chordLengthSquared());
                        
                        std::tuple<const NetworkNodeType*,const NetworkNodeType*,std::set<const PlanarMeshFace<dim>*>,float> key(std::tuple<const NetworkNodeType* const,const NetworkNodeType* const,std::set<const PlanarMeshFace<dim>*>,float>(nullptr,nullptr,std::set<const PlanarMeshFace<dim>*>(),-1.0));
                        if(periodicPrevNetwork==periodicNetworkSource)
                        {
                            const auto chord(periodicNext->get_P()-periodicPrev->get_P());
                            const auto chordNorm2(chord.squaredNorm());
                            const float u((loopNodePos-periodicPrev->get_P()).dot(chord)/chordNorm2);
                            key=std::make_tuple(periodicNetworkSource,periodicNetworkSink,periodicPatchEdge->meshIntersection->faces,u);
                        }
                        else if(periodicNextNetwork==periodicNetworkSource)
                        {
                            const auto chord(periodicPrev->get_P()-periodicNext->get_P());
                            const auto chordNorm2(chord.squaredNorm());
                            const float u((loopNodePos-periodicNext->get_P()).dot(chord)/chordNorm2);
                            key=std::make_tuple(periodicNetworkSource,periodicNetworkSink,periodicPatchEdge->meshIntersection->faces,u);
                        }
                        else
                        {
                            assert(false && "periodicNetworkSource must be periodicPrevNetwork or periodicNextNetwork");
                        }
                        //
                        //
                        //
                        //                        const auto key(std::make_tuple(periodicNetworkSource,periodicNetworkSink,u));
                        
                        
                        //                        std::cout<<"currentLoopLink="<<currentLoopLink->tag()<<std::endl;
                        //                        std::cout<<"u="<<std::setprecision(15)<<std::scientific<<std::get<3>(key)<<std::endl;
                        //                        std::cout<<"faces:";
                        //                        for(const auto& face : std::get<2>(key))
                        //                        {
                        //                            std::cout<<" "<<face->sID;
                        //                        }
                        //                        std::cout<<std::endl;
                        //                        std::cout<<"currentNetworkLink="<<currentNetworkLink->tag()<<std::endl;
                        //                        std::cout<<"networkNodePos="<<networkNodePos.transpose()<<std::endl;
                        //                        std::cout<<"currentNetworkLink->source->get_P()="<<currentNetworkLink->source->get_P().transpose()<<std::endl;
                        //                        std::cout<<"currentNetworkLink->chordLengthSquared()="<<currentNetworkLink->chordLengthSquared()<<std::endl;
                        
                        
                        const auto networkNodeIter(newNetworkNodesMap.find(key));
                        if(networkNodeIter!=newNetworkNodesMap.end())
                        {
                            //                            std::cout<<"using networkNode "<<networkNodeIter->second->tag()<<std::endl;
                            const auto newLoopNode(this->loopNodes().create(loop,networkNodeIter->second,loopNodePos,periodicPatch,periodicPatchEdge));
                            currentSource=this->expandLoopLink(*currentLoopLink,newLoopNode).get();
                        }
                        else
                        {
                            const auto newNetNode(this->networkNodes().create(networkNodePos,VectorDim::Zero(),0.0)); // TODO compute velocity and velocityReduction by interpolation
                            //                            std::cout<<"emplacing "<<currentNetworkLink->tag()<<"@"<<std::setprecision(15)<<std::scientific<<std::get<3>(key)<<", newNetNode="<<newNetNode->tag()<<std::endl;
                            newNetworkNodesMap.emplace(key,newNetNode);
                            const auto newLoopNode(this->loopNodes().create(loop,newNetNode,loopNodePos,periodicPatch,periodicPatchEdge));
                            currentSource=this->expandLoopLink(*currentLoopLink,newLoopNode).get();
                        }
                    }
                    p2= (p2+1)%polyInt.size();
                }
            }
        }
        
        danglingBoundaryLoopNodes.clear();
        
//        std::cout<<"Ouputing"<<std::endl;
//        DDconfigIO<dim> configIO(*this);
//        configIO.writeTxt(runID);
        
    }
    

    template <int dim, short unsigned int corder, typename InterpolationType>
    typename DislocationNetwork<dim,corder,InterpolationType>::VectorDim DislocationNetwork<dim,corder,InterpolationType>::displacement(const VectorDim& x) const
    {/*!\param[in] P position vector
      * \returns The stress field generated by the DislocationNetwork at P
      *
      * Note:
      */
        
        assert(!simulationParameters.isPeriodicSimulation());
        
        VectorDim temp(VectorDim::Zero());
        
        for(const auto& loop : this->loops())
        {// sum solid angle of each loop
            for(const auto& shift : periodicShifts)
            {
                temp-=loop.second.lock()->solidAngle(x+shift)/4.0/M_PI*loop.second.lock()->burgers();
            }
            //                if(!(loop.second->isVirtualBoundaryLoop() && simulationParameters.simulationType==DefectiveCrystalParameters::PERIODIC_IMAGES))
            //                {
            //                    for(const auto& shift : periodicShifts)
            //                    {
            //                        temp-=loop.second->solidAngle(x+shift)/4.0/M_PI*loop.second->burgers();
            //                    }
            //                }
        }
        
        for(const auto& link : this->networkLinks())
        {// sum line-integral part of displacement field per segment
            if(   !link.second.lock()->hasZeroBurgers()
               //                   && (!link.second->isBoundarySegment() || simulationParameters.simulationType==DefectiveCrystalParameters::FINITE_NO_FEM)
               //                   && !(link.second->isVirtualBoundarySegment() && simulationParameters.simulationType==DefectiveCrystalParameters::PERIODIC_IMAGES)
               )
            {
                for(const auto& shift : periodicShifts)
                {
                    temp+=link.second.lock()->straight.displacement(x+shift);
                }
            }
        }
        
        return temp;
    }
    
    /**********************************************************************/
    template <int dim, short unsigned int corder, typename InterpolationType>
    void DislocationNetwork<dim,corder,InterpolationType>::displacement(std::vector<FEMnodeEvaluation<ElementType,dim,1>>& fieldPoints) const
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
    template <int dim, short unsigned int corder, typename InterpolationType>
    typename DislocationNetwork<dim,corder,InterpolationType>::MatrixDim DislocationNetwork<dim,corder,InterpolationType>::stress(const VectorDim& x) const
    {/*!\param[in] P position vector
      * \returns The stress field generated by the DislocationNetwork at P
      *
      * Note:
      */
        MatrixDim temp(MatrixDim::Zero());
        for(const auto& link : this->networkLinks())
        {// sum stress field per segment
            if(   !link.second.lock()->hasZeroBurgers()
               && !(link.second.lock()->isBoundarySegment() && simulationParameters.simulationType==DefectiveCrystalParameters::FINITE_NO_FEM) // exclude boundary segments even if they are non-zero Burgers
               //                   && !(link.second->isVirtualBoundarySegment() && simulationParameters.simulationType==DefectiveCrystalParameters::PERIODIC_IMAGES)
               )
            {
                for(const auto& shift : periodicShifts)
                {
                    temp+=link.second.lock()->straight.stress(x+shift);
                }
            }
        }
        return temp;
    }
    
    /**********************************************************************/
    template <int dim, short unsigned int corder, typename InterpolationType>
    void DislocationNetwork<dim,corder,InterpolationType>::stress(std::deque<FEMfaceEvaluation<ElementType,dim,dim>>& fieldPoints) const
    {
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for(size_t k=0;k<fieldPoints.size();++k)
        {
            fieldPoints[k]=stress(fieldPoints[k].P);
        }
    }

    template <int dim, short unsigned int corder, typename InterpolationType>
    void DislocationNetwork<dim,corder,InterpolationType>::assembleAndSolveGlide(const long int& runID)
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
            const auto t1= std::chrono::system_clock::now();
            model::cout<<"        Computing analytical stress field at quadrature points ("<<nThreads<<" threads) "<<std::flush;
#ifdef _OPENMP
            EqualIteratorRange<typename LoopNetworkType::NetworkLinkContainerType::iterator> eir(this->networkLinks().begin(),this->networkLinks().end(),nThreads);
#pragma omp parallel for
//            for(size_t thread=0;thread<eir.size();thread++)
//            {
//                for (auto& linkIter=eir[thread].first;linkIter!=eir[thread].second;linkIter++)
//                {
//                    linkIter->second.lock()->assembleGlide();
//                }
//            }
            for(size_t k=0;k<this->networkLinks().size();++k)
            {
                auto linkIter(this->networkLinks().begin());
                std::advance(linkIter,k);
                linkIter->second.lock()->assembleGlide();
            }
#else
            for (auto& linkIter : this->networkLinks())
            {
                linkIter.second.lock()->assembleGlide();
            }
#endif
            model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t1)).count()<<" sec]."<<defaultColor<<std::endl;
        }
        else
        {// For curved segments use quandrature integration of stress field
            assert(0 && "ALL THIS MUST BE RE-IMPLEMENTED FOR CURVED SEGMENTS");
        }
        
        //! -3 Loop over DislocationSubNetworks, assemble subnetwork stiffness matrix and force vector, and solve
        model::cout<<"        Assembling and solving "<<std::flush;
        DislocationGlideSolver<LoopNetworkType>(*this).solve(runID);
    }

    template <int dim, short unsigned int corder, typename InterpolationType>
    void DislocationNetwork<dim,corder,InterpolationType>::updateBoundaryNodes()
    {
        
        if(danglingBoundaryLoopNodes.size())
        {
            std::cout<<"Removing bnd Nodes"<<std::endl;
            for(const auto& node : danglingBoundaryLoopNodes)
            {
                this->removeLoopNode(node->sID);
            }
            
            std::cout<<"Inserting new boundary nodes"<<std::endl;
            std::map<std::tuple<const NetworkNodeType*,const NetworkNodeType*,std::set<const PlanarMeshFace<dim>*>,float>,std::shared_ptr<NetworkNodeType>> newNetworkNodesMap;
            for(const auto& weakLoop : this->loops())
            {
                const auto loop(weakLoop.second.lock());
                std::vector<std::pair<VectorLowerDim,const LoopNodeType* const>> loopNodesPos;
                std::map<std::tuple<const LoopNodeType*,const LoopNodeType*,const PeriodicPlaneEdge<dim>*>,const LoopNodeType*> bndNodesMap;
                for(const auto& loopLink : loop->linkSequence())
                {
                    if(!loopLink->source->periodicPlaneEdge)
                    {
                        loopNodesPos.emplace_back(loop->periodicGlidePlane->referencePlane->localPosition(loopLink->source->get_P()),loopLink->source.get());
                    }
                    else
                    {
                        bndNodesMap.emplace(std::make_tuple(loopLink->source->periodicPrev(),loopLink->source->periodicNext(),loopLink->source->periodicPlaneEdge.get()),loopLink->source.get());
                    }
                }
                
                const auto polyInt(loop->periodicGlidePlane->polygonPatchIntersection(loopNodesPos)); // 2dPos,shift,edgeID,LoopNodeType*
                std::map<const LoopNodeType*,size_t> polyIntMap;
                for(size_t p=0;p<polyInt.size();++p)
                {
                    if(std::get<3>(polyInt[p]))
                    {
                        polyIntMap.emplace(std::get<3>(polyInt[p]),p);
                    }
                }
                
                for(size_t k=0;k<loopNodesPos.size();++k)
                {
                    const auto periodicPrev(loopNodesPos[k].second);
                    const auto periodicPrevNetwork(periodicPrev->networkNode.get());
                    const auto polyIter(polyIntMap.find(periodicPrev));
                    assert(polyIter!=polyIntMap.end());
                    const size_t p(polyIter->second);
                    
                    const size_t k1(k<loopNodesPos.size()-1? k+1 : 0);
                    const auto periodicNext(loopNodesPos[k1].second);
                    const auto periodicNextNetwork(periodicNext->networkNode.get());
                    const auto polyIter1(polyIntMap.find(periodicNext));
                    assert(polyIter1!=polyIntMap.end());
                    const size_t p1(polyIter1->second);
                    
                    const auto periodicNetworkSource(periodicPrevNetwork->sID<periodicNextNetwork->sID? periodicPrevNetwork : periodicNextNetwork);
                    const auto periodicNetworkSink  (periodicPrevNetwork->sID<periodicNextNetwork->sID? periodicNextNetwork : periodicPrevNetwork);
                    
                    const LoopNodeType* currentSource(periodicPrev);
                    size_t p2 = (p+1)%polyInt.size();
                    while(p2<p1)
                    {
                        const auto periodicPatch(loop->periodicGlidePlane->getPatch(std::get<1>(polyInt[p2])));
                        const auto periodicPatchEdge(periodicPatch->edges()[std::get<2>(polyInt[p2])]);
                        const auto bndIter(bndNodesMap.find(std::make_tuple(periodicPrev,periodicNext,periodicPatchEdge.get())));
                        
                        if(bndIter!=bndNodesMap.end())
                        {// exising bnd node found
                            currentSource=bndIter->second;
                        }
                        else
                        {
                            const VectorDim    loopNodePos(loop->periodicGlidePlane->referencePlane->globalPosition(std::get<0>(polyInt[p2])));
                            const VectorDim networkNodePos(loopNodePos+std::get<1>(polyInt[p2]));
                            const auto currentLoopLink(currentSource->next.second);
                            const auto currentNetworkLink(currentLoopLink->networkLink());
                            
                            std::tuple<const NetworkNodeType*,const NetworkNodeType*,std::set<const PlanarMeshFace<dim>*>,float> key(std::tuple<const NetworkNodeType* const,const NetworkNodeType* const,std::set<const PlanarMeshFace<dim>*>,float>(nullptr,nullptr,std::set<const PlanarMeshFace<dim>*>(),-1.0));
                            if(periodicPrevNetwork==periodicNetworkSource)
                            {
                                const auto chord(periodicNext->get_P()-periodicPrev->get_P());
                                const auto chordNorm2(chord.squaredNorm());
                                const float u((loopNodePos-periodicPrev->get_P()).dot(chord)/chordNorm2);
                                key=std::make_tuple(periodicNetworkSource,periodicNetworkSink,periodicPatchEdge->meshIntersection->faces,u);
                            }
                            else if(periodicNextNetwork==periodicNetworkSource)
                            {
                                const auto chord(periodicPrev->get_P()-periodicNext->get_P());
                                const auto chordNorm2(chord.squaredNorm());
                                const float u((loopNodePos-periodicNext->get_P()).dot(chord)/chordNorm2);
                                key=std::make_tuple(periodicNetworkSource,periodicNetworkSink,periodicPatchEdge->meshIntersection->faces,u);
                            }
                            else
                            {
                                assert(false && "periodicNetworkSource must be periodicPrevNetwork or periodicNextNetwork");
                            }
                            
                            const auto networkNodeIter(newNetworkNodesMap.find(key));
                            if(networkNodeIter!=newNetworkNodesMap.end())
                            {
                                //                            std::cout<<"using networkNode "<<networkNodeIter->second->tag()<<std::endl;
                                const auto newLoopNode(this->loopNodes().create(loop,networkNodeIter->second,loopNodePos,periodicPatch,periodicPatchEdge));
                                currentSource=this->expandLoopLink(*currentLoopLink,newLoopNode).get();
                            }
                            else
                            {
                                const auto newNetNode(this->networkNodes().create(networkNodePos,VectorDim::Zero(),0.0)); // TODO compute velocity and velocityReduction by interpolation
                                //                            std::cout<<"emplacing "<<currentNetworkLink->tag()<<"@"<<std::setprecision(15)<<std::scientific<<std::get<3>(key)<<", newNetNode="<<newNetNode->tag()<<std::endl;
                                newNetworkNodesMap.emplace(key,newNetNode);
                                const auto newLoopNode(this->loopNodes().create(loop,newNetNode,loopNodePos,periodicPatch,periodicPatchEdge));
                                currentSource=this->expandLoopLink(*currentLoopLink,newLoopNode).get();
                            }
                        }
                        p2= (p2+1)%polyInt.size();
                    }
                }
            }
            danglingBoundaryLoopNodes.clear();
        }
    }
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    void DislocationNetwork<dim,corder,InterpolationType>::moveGlide(const double & dt_in)
    {/*! Moves all nodes in the DislocationNetwork using the stored velocity and current dt
      */
        
        
        if(simulationParameters.isPeriodicSimulation())
        {
            

            danglingBoundaryLoopNodes.clear();
            
            for(auto& node : this->networkNodes())
            {// Expansion
                node.second.lock()->trySet_P(node.second.lock()->get_P()+node.second.lock()->get_V()*dt_in);
            }
            updateBoundaryNodes();
            
        }
        else
        {
            model::cout<<"        Moving DislocationNodes by glide (dt="<<dt_in<< ")... "<<std::flush;
            const auto t0= std::chrono::system_clock::now();
            for (auto& nodeIter : this->networkNodes())
            {
//                nodeIter.second.lock()->moveGlide(dt_in);
                nodeIter.second.lock()->set_P(nodeIter.second.lock()->get_P()+nodeIter.second.lock()->get_V()*dt_in);
            }
            model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<defaultColor<<std::endl;
        }
        
    }

    /**********************************************************************/
    template <int dim, short unsigned int corder, typename InterpolationType>
    void DislocationNetwork<dim,corder,InterpolationType>::singleGlideStepDiscreteEvents(const long int& runID)
    {
        
        //#ifdef DislocationNucleationFile
        //            if(use_bvp && !(runID%use_bvp))
        //            {
        //                nucleateDislocations(); // needs to be called before updateQuadraturePoints()
        //                updateQuadraturePoints();
        //            }
        //#endif
        //            assembleAndSolve(runID,straightSegmentsDeq);
        //            computeNodaVelocities(runID);
        
        //! 10- Cross Slip (needs upated PK force)
//        DislocationCrossSlip<DislocationNetworkType> crossSlip(*this);
//        crossSlip.execute();
//
//
//        //                        gbTransmission.transmit();
//        //gbTransmission.directTransmit();
//        //            GrainBoundaryDissociation<DislocationNetworkType>(*this).dissociate();
//        //            poly.reConnectGrainBoundarySegments(*this); // this makes stressGauss invalid, so must follw other GB operations
//
//
//        //! 12- Form Junctions
        junctionsMaker.formJunctions(DDtimeIntegrator<0>::dxMax);
//
//        //            // Remesh may contract juncitons to zero lenght. Remove those juncitons:
//        //            DislocationJunctionFormation<DislocationNetworkType>(*this).breakZeroLengthJunctions();
//
//
//
//        //! 13- Node redistribution
        networkRemesher.remesh(runID);
//
//        updateVirtualBoundaryLoops();
        
    }


    template <int dim, short unsigned int corder, typename InterpolationType>
    int DislocationNetwork<dim,corder,InterpolationType>::verboseDislocationNetwork=0;

    
}
#endif
