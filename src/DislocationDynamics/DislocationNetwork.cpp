/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationNetwork_CPP_
#define model_DislocationNetwork_CPP_

#include <numbers>

#include <DislocationNetwork.h>
//#include <PolyhedronInclusionNodes.h>

namespace model
{

/**********************************************************************/
template <int dim, short unsigned int corder>
DislocationNetwork<dim,corder>::DislocationNetwork(const DefectiveCrystalParameters& _simulationParameters,
                                                   const SimplicialMesh<dim>& _mesh,
                                                   const Polycrystal<dim>& _poly,
                                                   const std::unique_ptr<BVPsolver<dim,2>>& _bvpSolver,
                                                   const std::unique_ptr<ExternalLoadControllerBase<dim>>& _externalLoadController,
                                                   const std::vector<VectorDim>& _periodicShifts,
                                                   long int& runID) :
/* init */ simulationParameters(_simulationParameters)
/* init */,mesh(_mesh)
/* init */,poly(_poly)
/* init */,glidePlaneFactory(poly)
/* init */,periodicGlidePlaneFactory(simulationParameters.isPeriodicSimulation()? new PeriodicGlidePlaneFactory<dim>(poly, glidePlaneFactory) : nullptr)
/* init */,bvpSolver(_bvpSolver)
/* init */,externalLoadController(_externalLoadController)
/* init */,periodicShifts(_periodicShifts)
/* init */,networkRemesher(*this)
/* init */,junctionsMaker(*this)
/* init */,crossSlipModel(DislocationCrossSlip<DislocationNetwork<dim,corder>>::getModel(poly,simulationParameters.traitsIO))
/* init */,crossSlipMaker(*this)
/* init */,nodeContractor(*this)
/* init */,timeIntegrator(simulationParameters.traitsIO.ddFile)
/* init */,stochasticForceGenerator(simulationParameters.use_stochasticForce? new StochasticForceGenerator(simulationParameters.traitsIO) : nullptr)
/* init */,networkIO(*this)
/* init */,ddSolverType(TextFileParser(simulationParameters.traitsIO.ddFile).readScalar<int>("ddSolverType",true))
/* init */,computeDDinteractions(TextFileParser(simulationParameters.traitsIO.ddFile).readScalar<int>("computeDDinteractions",true))
/* init */,outputFrequency(TextFileParser(simulationParameters.traitsIO.ddFile).readScalar<int>("outputFrequency",true))
/* init */,outputBinary(TextFileParser(simulationParameters.traitsIO.ddFile).readScalar<int>("outputBinary",true))
/* init */,outputMeshDisplacement(TextFileParser(simulationParameters.traitsIO.ddFile).readScalar<int>("outputMeshDisplacement",true))
/* init */,outputFEMsolution(TextFileParser(simulationParameters.traitsIO.ddFile).readScalar<int>("outputFEMsolution",true))
/* init */,outputQuadraturePoints(TextFileParser(simulationParameters.traitsIO.ddFile).readScalar<int>("outputQuadraturePoints",true))
/* init */,outputLinkingNumbers(TextFileParser(simulationParameters.traitsIO.ddFile).readScalar<int>("outputLinkingNumbers",true))
/* init */,outputLoopLength(TextFileParser(simulationParameters.traitsIO.ddFile).readScalar<int>("outputLoopLength",true))
/* init */,outputSegmentPairDistances(TextFileParser(simulationParameters.traitsIO.ddFile).readScalar<int>("outputSegmentPairDistances",true))
/* init */,computeElasticEnergyPerLength(TextFileParser(simulationParameters.traitsIO.ddFile).readScalar<int>("computeElasticEnergyPerLength",true))
/* init */,alphaLineTension(TextFileParser(simulationParameters.traitsIO.ddFile).readScalar<double>("alphaLineTension",true))
/* init */,use_velocityFilter(TextFileParser(simulationParameters.traitsIO.ddFile).readScalar<double>("use_velocityFilter",true))
/* init */,velocityReductionFactor(TextFileParser(simulationParameters.traitsIO.ddFile).readScalar<double>("velocityReductionFactor",true))
/* init */,verboseDislocationNode(TextFileParser(simulationParameters.traitsIO.ddFile).readScalar<int>("verboseDislocationNode",true))
/* init */,capMaxVelocity(TextFileParser(simulationParameters.traitsIO.ddFile).readScalar<int>("capMaxVelocity",true))
{
    
    assert(velocityReductionFactor>0.0 && velocityReductionFactor<=1.0);
        
    // Initialize static variables
    LoopNetworkType::verboseLevel=TextFileParser(simulationParameters.traitsIO.ddFile).readScalar<int>("verboseLoopNetwork",true);
    verboseDislocationNetwork=TextFileParser(simulationParameters.traitsIO.ddFile).readScalar<int>("verboseDislocationNetwork",true);
    LoopType::initFromFile(simulationParameters.traitsIO.ddFile);
    LoopNodeType::initFromFile(simulationParameters.traitsIO.ddFile);
    LoopLinkType::initFromFile(simulationParameters.traitsIO.ddFile);
    NetworkLinkType::initFromFile(simulationParameters.traitsIO.ddFile);
    DislocationFieldBase<dim>::initFromFile(simulationParameters.traitsIO.ddFile);
    
    DDconfigIO<dim> evl(simulationParameters.traitsIO.evlFolder);
    evl.read(runID);
    setConfiguration(evl);
}


//New Version
template <int dim, short unsigned int corder>
void DislocationNetwork<dim,corder>::setConfiguration(const DDconfigIO<dim>& evl)
{
    this->loopLinks().clear(); // erase base network to clear current config
    
    // Create Loops
    std::deque<std::shared_ptr<LoopType>> tempLoops; // keep loops alive during setConfiguration
    size_t loopNumber=1;
    for(const auto& loop : evl.loops())
    {
        const bool faulted(poly.grain(loop.grainID).singleCrystal->rationalLatticeDirection(loop.B).rat.asDouble()!=1.0? true : false);
        std::cout<<"Creating DislocationLoop "<<loop.sID<<" ("<<loopNumber<<" of "<<evl.loops().size()<<"), type="<<loop.loopType<<", faulted="<<faulted<<", |b|="<<loop.B.norm()<<std::endl;
        const size_t loopIDinFile(loop.sID);
        LoopType::set_count(loopIDinFile);
        
        GlidePlaneKey<dim> loopPlaneKey(loop.P, poly.grain(loop.grainID).singleCrystal->reciprocalLatticeDirection(loop.N));
        tempLoops.push_back(this->loops().create(loop.B, glidePlaneFactory.getFromKey(loopPlaneKey)));
        assert(this->loops().get(loopIDinFile)->sID == loopIDinFile);
        loopNumber++;
        
        
        //            switch (loop.loopType)
        //            {
        //                case DislocationLoopIO<dim>::GLISSILELOOP:
        //                {
        //                    GlidePlaneKey<dim> loopPlaneKey(loop.P, poly.grain(loop.grainID).reciprocalLatticeDirection(loop.N));
        //                    tempLoops.push_back(this->loops().create(loop.B, glidePlaneFactory.getFromKey(loopPlaneKey)));
        //                    assert(this->loops().get(loopIDinFile)->sID == loopIDinFile);
        //                    loopNumber++;
        //                    break;
        //                }
        //                case DislocationLoopIO<dim>::SESSILELOOP:
        //                {
        //                    tempLoops.push_back(this->loops().create(loop.B,loop.grainID,loop.loopType ));
        //                    assert(this->loops().get(loopIDinFile)->sID == loopIDinFile);
        //                    loopNumber++;
        //                    break;
        //                }
        //                default:
        //                    assert(false && "Unknown DislocationLoop type");
        //                    break;
        //            }
    }
    
    // Create NetworkNodes
    std::deque<std::shared_ptr<NetworkNodeType>> tempNetNodes; // keep loops alive during setConfiguration
    size_t netNodeNumber=1;
    for(const auto& node : evl.nodes())
    {
        std::cout<<"Creating DislocationNode "<<node.sID<<" ("<<netNodeNumber<<" of "<<evl.nodes().size()<<")"<<std::endl;
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
        std::cout<<"Creating DislocationLoopNode "<<node.sID<<" ("<<loopNodeNumber<<" of "<<evl.loopNodes().size()<<")"<<std::flush;
        // std::cout<<"Printing stuff for loopNodes "<<std::endl;
        
        // std::cout << node.loopID << std::endl;
        // std::cout << node.sID << std::endl;
        // std::cout << std::setprecision(15) << std::scientific << node.P.transpose() << std::endl;
        // std::cout << node.networkNodeID << std::endl;
        // std::cout << std::setprecision(15) << std::scientific << node.periodicShift.transpose() << std::endl;
        // std::cout << node.edgeID << std::endl;
        
        const size_t nodeIDinFile(node.sID);
        LoopNodeType::set_count(nodeIDinFile);
        const auto loop(this->loops().get(node.loopID));
        const auto netNode(this->networkNodes().get(node.networkNodeID));
        assert(loop && "Loop does not exist");
        assert(netNode && "NetworkNode does not exist");
        //            std::set<std::shared_ptr<PeriodicPlanePatch<dim>>> auxPatchSets; //For inserting patches corresponding to diagonally opposite itnersections
        const auto periodicPatch(loop->periodicGlidePlane? loop->periodicGlidePlane->patches().getFromKey(node.periodicShift) : nullptr);
        const auto periodicPatchEdge((periodicPatch && node.edgeIDs.first>=0)? (node.edgeIDs.second>=0 ? std::make_pair(periodicPatch->edges()[node.edgeIDs.first],
                                                                                                                        periodicPatch->edges()[node.edgeIDs.second]): std::make_pair(periodicPatch->edges()[node.edgeIDs.first],nullptr)):std::make_pair(nullptr,nullptr));
        // std::cout<<"PeriodicPlane edge created "<<std::endl;
        if(periodicPatch)
        {
            std::cout<<", on patch "<<periodicPatch->shift.transpose()<<std::flush;
        }
        if(periodicPatchEdge.first)
        {
            std::cout<<" on edge "<<periodicPatchEdge.first->edgeID<<std::flush;
        }
        if(periodicPatchEdge.second)
        {
            std::cout<<" and "<<periodicPatchEdge.second->edgeID<<std::flush;
            //                auxPatchSets.emplace(loop->periodicGlidePlane->patches().getFromKey(periodicPatch->shift + periodicPatchEdge.first->deltaShift));
            //                auxPatchSets.emplace(loop->periodicGlidePlane->patches().getFromKey(periodicPatch->shift + periodicPatchEdge.second->deltaShift));
        }
        //            if (auxPatchSets.size()==0)
        //            {
        //                auxPatchSets.insert(std::shared_ptr<PeriodicPlanePatch<dim>>(nullptr));
        //            }
        // std::cout<<"Trying to create the loop node with loopID "<<loop->sID<<" if loop has GP "<<(loop->glidePlane!=nullptr)<<" networkID is "
        // <<netNode->sID<<std::endl;
        std::cout<<std::endl;
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
        std::cout<<" Inserting loop "<<loop.sID<<std::endl;
        this->insertLoop(this->loops().get(loop.sID),loopNodes);
    }
    updateGeometry();
    
    
    
    for(const auto& inclusion : evl.sphericalInclusions())
    {
        std::cout<<"Creating spherical inclusion "<<inclusion.inclusionID<<std::endl;
        const std::pair<bool,const Simplex<dim,dim>*> searchPair(mesh.search(inclusion.C));
        if(searchPair.first)
        {
            
            const auto& grain(poly.grain(searchPair.second->region->regionID));
            if(inclusion.phaseID<int(grain.singleCrystal->secondPhases().size()))
            {
            const auto secondPhase(grain.singleCrystal->secondPhases()[inclusion.phaseID]);
            EshelbyInclusionBase<dim>::set_count(inclusion.inclusionID);
            
            
            std::shared_ptr<EshelbyInclusionBase<dim>> iptr(new SphericalInclusion<dim>(inclusion.C,inclusion.a,inclusion.eT,poly.nu,poly.mu,inclusion.mobilityReduction,inclusion.phaseID,secondPhase));
            
            eshelbyInclusions().emplace(inclusion.inclusionID,iptr);
            }
            else
            {
                throw std::runtime_error("phaseID does not exist in grain.");
            }
            
//            eshelbyInclusions().emplace(std::piecewise_construct,
//                                        std::make_tuple(inclusion.inclusionID),
//                                        std::make_tuple(inclusion.C,inclusion.a,inclusion.eT,poly.nu,poly.mu,inclusion.mobilityReduction,inclusion.phaseID,secondPhase) );
            
        }
    }
    
    for(const auto& piNode : evl.polyhedronInclusionNodes())
    {
        polyhedronInclusionNodes().emplace(piNode.nodeID,piNode);
    }
    
    std::map<size_t,std::map<size_t,std::vector<size_t>>> faceMap;
    for(const auto& edge : evl.polyhedronInclusionEdges())
    {
        const size_t& iID(edge.inclusionID);
        const size_t& fID(edge.faceID);
        const size_t& sourceID(edge.sourceID);
        faceMap[iID][fID].push_back(sourceID);
    }

    
    for(const auto& inclusion : evl.polyhedronInclusions())
    {
        std::cout<<"Creating polyhedron inclusion "<<inclusion.inclusionID<<std::endl;

        const auto faceIter(faceMap.find(inclusion.inclusionID));
        if(faceIter!=faceMap.end())
        {
            const auto& faces(faceIter->second);
            std::cout<<"    #faces= "<<faces.size()<<std::endl;
            std::set<const PolyhedronInclusionNodeIO<dim>*> uniquePolyNodes;
            for(const auto& pair : faces)
            {
                for(const auto& nodeID : pair.second)
                {
                    uniquePolyNodes.emplace(&polyhedronInclusionNodes().at(nodeID));
                }
            }
            std::cout<<"    #nodes= "<<uniquePolyNodes.size()<<std::endl;
            if(uniquePolyNodes.size()>=dim+1)
            {
                // Find grain
                std::set<size_t> grainIDs;
                for(const auto& nodePtr : uniquePolyNodes)
                {
                    const std::pair<bool,const Simplex<dim,dim>*> searchPair(mesh.search(nodePtr->P));
                    if(searchPair.first)
                    {
                        grainIDs.insert(searchPair.second->region->regionID);
                    }
                    else
                    {
                        throw std::runtime_error("inclusion node outside mesh");
                    }
                }
                
                // Add inclusion
                if(grainIDs.size()==1)
                {
                    const auto& grain(poly.grain(*grainIDs.begin()));
                    if(inclusion.phaseID<int(grain.singleCrystal->secondPhases().size()))
                    {
                    const auto secondPhase(grain.singleCrystal->secondPhases()[inclusion.phaseID]);
                    EshelbyInclusionBase<dim>::set_count(inclusion.inclusionID);
                    std::shared_ptr<EshelbyInclusionBase<dim>> iptr(new PolyhedronInclusion<dim>( polyhedronInclusionNodes(),faces,inclusion.eT,poly.nu,poly.mu,inclusion.mobilityReduction,inclusion.phaseID,secondPhase));
                    eshelbyInclusions().emplace(inclusion.inclusionID,iptr);
                    }
                    else
                    {
                        throw std::runtime_error("phaseID does not exist in grain.");
                    }
                }
                else
                {
                    throw std::runtime_error("inclusion across grain boundaries");
                }
            }
            else
            {
                throw std::runtime_error("inclusion does not have enough nodes");
            }
        }
        else
        {
            throw std::runtime_error("inclusionID not found in faceMap");
        }
    }
        
}

/**********************************************************************/
template <int dim, short unsigned int corder>
const typename DislocationNetwork<dim,corder>::EshelbyInclusionContainerType& DislocationNetwork<dim,corder>::eshelbyInclusions() const
{
    return *this;
}

template <int dim, short unsigned int corder>
typename DislocationNetwork<dim,corder>::EshelbyInclusionContainerType& DislocationNetwork<dim,corder>::eshelbyInclusions()
{
    return *this;
}

template <int dim, short unsigned int corder>
const typename DislocationNetwork<dim,corder>::PolyhedronInclusionNodeContainerType& DislocationNetwork<dim,corder>::polyhedronInclusionNodes() const
{
    return *this;
}

template <int dim, short unsigned int corder>
typename DislocationNetwork<dim,corder>::PolyhedronInclusionNodeContainerType& DislocationNetwork<dim,corder>::polyhedronInclusionNodes()
{
    return *this;
}

/**********************************************************************/
//    template <int dim, short unsigned int corder>
//    void DislocationNetwork<dim,corder>::createEshelbyInclusions()
//    {
//        for(const auto& grain : poly.grains())
//        {
//            EshelbyInclusion<dim>::addSlipSystems(grain.second.slipSystems());
//        }
//
//
//        IDreader<'E',1,14,double> inclusionsReader;
//        inclusionsReader.read(0,true);
//
//        const std::vector<double> inclusionsMobilityReduction(TextFileParser("./inputFiles/initialMicrostructure.txt").readArray<double>("inclusionsMobilityReduction",true));
//        for(const auto& pair : inclusionsReader)
//        {
//
//            const size_t& inclusionID(pair.first);
//            Eigen::Map<const Eigen::Matrix<double,1,14>> row(pair.second.data());
//
//            const VectorDim C(row.template segment<dim>(0));
//            const double a(row(dim+0));
//            MatrixDim eT(MatrixDim::Zero());
//            const int typeID(row(13));
//            int k=dim+1;
//            for(int i=0;i<dim;++i)
//            {
//                for(int j=0;j<dim;++j)
//                {
//                    eT(i,j)=row(k);
//                    k++;
//                }
//            }
//
//
//
//            EshelbyInclusion<dim>::set_count(inclusionID);
//            eshelbyInclusions().emplace(std::piecewise_construct,
//                                        std::make_tuple(inclusionID),
//                                        std::make_tuple(C,a,eT,poly.nu,poly.mu,inclusionsMobilityReduction[typeID],typeID) );
//        }
//    }

template <int dim, short unsigned int corder>
void DislocationNetwork<dim,corder>::updateGeometry()
{
    VerboseDislocationNetwork(2,"DislocationNetwork::updateGeometry"<<std::endl;);
    for(auto& loop : this->loops())
    {// copmute slipped areas and right-handed normal // TODO: PARALLELIZE THIS LOOP
        loop.second.lock()->updateGeometry();
    }
    // updatePlasticDistortionRateFromAreas();
    VerboseDislocationNetwork(3,"DislocationNetwork::updateGeometry DONE"<<std::endl;);
}

template <int dim, short unsigned int corder>
void DislocationNetwork<dim,corder>::updateRates()
{
    VerboseDislocationNetwork(2,"DislocationNetwork::updateRates"<<std::endl;);
    for(auto& loop : this->loops())
    {// copmute slipped areas and right-handed normal // TODO: PARALLELIZE THIS LOOP
        loop.second.lock()->updateRates();
    }
    // updatePlasticDistortionRateFromAreas();
    VerboseDislocationNetwork(3,"DislocationNetwork::updateRates DONE"<<std::endl;);
}


template <int dim, short unsigned int corder>
typename DislocationNetwork<dim,corder>::DislocationNetworkIOType& DislocationNetwork<dim,corder>::io()
{
    //        return DislocationNetworkIOType(*this);
    return networkIO;
}

template <int dim, short unsigned int corder>
const typename DislocationNetwork<dim,corder>::DislocationNetworkIOType& DislocationNetwork<dim,corder>::io() const
{
    //        return DislocationNetworkIOType(*this);
    return networkIO;
}

template <int dim, short unsigned int corder>
typename DislocationNetwork<dim,corder>::MatrixDim DislocationNetwork<dim,corder>::plasticDistortion() const
{
    MatrixDim temp(MatrixDim::Zero());
    for(const auto& loop : this->loops())
    {
        temp+= loop.second.lock()->plasticDistortion();
    }
    return temp;
}

template <int dim, short unsigned int corder>
std::map<std::pair<int,int>,double> DislocationNetwork<dim,corder>::slipSystemPlasticDistortion() const
{
    std::map<std::pair<int,int>,double> temp; // <grainID,slipSystemID>
    for(const auto& weakloop : this->loops())
    {
        const auto loop(weakloop.second.lock());
        if(loop->slipSystem())
        {
            const double val = loop->slipSystem()->unitSlip.transpose()*loop->plasticDistortion()*loop->slipSystem()->unitNormal;
            const std::pair<int,int> key(std::make_pair(loop->grain.grainID,loop->slipSystem()->sID));
            auto iter(temp.find(key));
            if(iter!=temp.end())
            {// iter found
                iter->second += val ;
            }
            else
            {// iter not found
                temp.emplace(key,val);
            }
        }
    }
    return temp;
}

template <int dim, short unsigned int corder>
typename DislocationNetwork<dim,corder>::MatrixDim DislocationNetwork<dim,corder>::plasticDistortionRate() const
{
    MatrixDim temp(MatrixDim::Zero());
    for(const auto& loop : this->loops())
    {
        temp+= loop.second.lock()->plasticDistortionRate();
    }
    return temp;
}

// template <int dim, short unsigned int corder>
// const typename DislocationNetwork<dim,corder>::MatrixDim& DislocationNetwork<dim,corder>::plasticDistortionRate() const
// {
//     return  _plasticDistortionRateFromAreas;
// }

template <int dim, short unsigned int corder>
typename DislocationNetwork<dim,corder>::MatrixDim DislocationNetwork<dim,corder>::plasticStrain() const
{/*!\returns the plastic strain rate tensor generated during the last time step.
  */
    return (plasticDistortion()+plasticDistortion().transpose())*0.5;
}

template <int dim, short unsigned int corder>
typename DislocationNetwork<dim,corder>::MatrixDim DislocationNetwork<dim,corder>::plasticStrainRate() const
{/*!\returns the plastic strain rate tensor generated during the last time step.
  */
    return (plasticDistortionRate()+plasticDistortionRate().transpose())*0.5;
}

template <int dim, short unsigned int corder>
std::tuple<double,double,double,double> DislocationNetwork<dim,corder>::networkLength() const
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


template <int dim, short unsigned int corder>
bool DislocationNetwork<dim,corder>::contract(std::shared_ptr<NetworkNodeType> nA,
                                              std::shared_ptr<NetworkNodeType> nB)
{
    return nodeContractor.contract(nA,nB);
}





template <int dim, short unsigned int corder>
typename DislocationNetwork<dim,corder>::VectorDim DislocationNetwork<dim,corder>::displacement(const VectorDim& x) const
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
            temp-=loop.second.lock()->solidAngle(x+shift)/4.0/std::numbers::pi*loop.second.lock()->burgers();
        }
    }
    
    for(const auto& link : this->networkLinks())
    {// sum line-integral part of displacement field per segment
        if(   !link.second.lock()->hasZeroBurgers()
           //                   && (!link.second->isBoundarySegment() || simulationParameters.simulationType==DDtraitsIO::FINITE_NO_FEM)
           //                   && !(link.second->isVirtualBoundarySegment() && simulationParameters.simulationType==DDtraitsIO::PERIODIC_IMAGES)
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
template <int dim, short unsigned int corder>
void DislocationNetwork<dim,corder>::displacement(std::vector<FEMnodeEvaluation<ElementType,dim,1>>& fieldPoints) const
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
template <int dim, short unsigned int corder>
typename DislocationNetwork<dim,corder>::MatrixDim DislocationNetwork<dim,corder>::stress(const VectorDim& x) const
{/*!\param[in] P position vector
  * \returns The stress field generated by the DislocationNetwork at P
  *
  * Note:
  */
    MatrixDim temp(MatrixDim::Zero());
    for(const auto& link : this->networkLinks())
    {// sum stress field per segment
        if(   !link.second.lock()->hasZeroBurgers()
           && !(link.second.lock()->isBoundarySegment() && simulationParameters.simulationType==DDtraitsIO::FINITE_NO_FEM) // exclude boundary segments even if they are non-zero Burgers
           //                   && !(link.second->isVirtualBoundarySegment() && simulationParameters.simulationType==DDtraitsIO::PERIODIC_IMAGES)
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
template <int dim, short unsigned int corder>
void DislocationNetwork<dim,corder>::stress(std::deque<FEMfaceEvaluation<ElementType,dim,dim>>& fieldPoints) const
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

template <int dim, short unsigned int corder>
void DislocationNetwork<dim, corder>::solveGlide(const long int &runID)
{
//    const auto t0 = std::chrono::system_clock::now();
//    std::cout <<"Solving glide "<< std::flush;
    DislocationGlideSolver<LoopNetworkType>(*this).solve(runID);
//    std::cout << magentaColor << std::setprecision(3) << std::scientific << " [" << (std::chrono::duration<double>(std::chrono::system_clock::now() - t0)).count() << " sec]." << defaultColor << std::endl;
}

template <int dim, short unsigned int corder>
void DislocationNetwork<dim, corder>::assembleGlide(const long int &runID, const double &maxVelocity)
{ /*! Performs the following operatons:
   */
#ifdef _OPENMP
    const size_t nThreads = omp_get_max_threads();
#else
    const size_t nThreads = 1;
#endif
    
    //! -1 Compute the interaction StressField between dislocation particles
    
    std::map<int, int> velocityBinMap;
    for (const auto &binVal : simulationParameters.subcyclingBins)
    {
        velocityBinMap.emplace(binVal, 0);
    }
    if (corder == 0)
    { // For straight segments use analytical expression of stress field
        const auto t0 = std::chrono::system_clock::now();
        std::cout <<"Updating quadrature points "<< std::flush;
        for (const auto &links : this->networkLinks())
        {
            
            const int velGroup(simulationParameters.useSubCycling ? links.second.lock()->velocityGroup(maxVelocity, simulationParameters.subcyclingBins) : 1);
            auto velocityBinIter(velocityBinMap.find(velGroup));
            assert(velocityBinIter != velocityBinMap.end());
            velocityBinIter->second++;
            
            if ((runID % velGroup) == 0)
            {
                links.second.lock()->updateQuadraturePointsSeg();
            }
        }
        std::cout << magentaColor << std::setprecision(3) << std::scientific << " [" << (std::chrono::duration<double>(std::chrono::system_clock::now() - t0)).count() << " sec]." << defaultColor << std::endl;
        
        
        const auto t1 = std::chrono::system_clock::now();
        std::cout <<"Computing stacking fault forces at quadrature points (" << nThreads << " threads) " << std::flush;
//        const double eps=1.0e-2;
#ifdef _OPENMP
#pragma omp parallel for // THERE COULD BE A WRITING RACE HERE SINCE DIFFERENT LOOPS MAY TRY TO SUM TO THE SAME GAUSS POINT
#endif
        for (size_t k = 0; k < this->loops().size(); ++k)
        {
            auto loopIter(this->loops().begin());
            std::advance(loopIter, k);
            const auto& fieldLoop(loopIter->second.lock());
            fieldLoop->computeStackingFaultForces();
//            if(fieldLoop->slipSystem() && fieldLoop->glidePlane)
//            {
//                for(const auto& loopLink : fieldLoop->loopLinks())
//                {
//                    if(loopLink->networkLink())
//                    {
//
//                        VectorDim outDir((loopLink->sink->get_P() - loopLink->source->get_P()).cross(fieldLoop->rightHandedUnitNormal()));
//                        const double outDirNorm(outDir.norm());
//                        if(outDirNorm>FLT_EPSILON)
//                        {
//                            outDir/=outDirNorm;
//                            std::vector<std::pair<VectorDim,VectorDim>> qPointSlip(loopLink->networkLink()->quadraturePoints().size(),std::make_pair(VectorDim::Zero(),VectorDim::Zero())); // accumulated slip vectors (outside,inside) for each qPoint
//
//                            for(const auto& weakSourceLoop : this->loops())
//                            {
//                                const auto sourceLoop(weakSourceLoop.second.lock());
//                                if(sourceLoop->slipSystem())
//                                {
//                                    if(fieldLoop->slipSystem()->n==sourceLoop->slipSystem()->n)
//                                    {// same glide plane family
//
//                                        for(size_t q=0;q<loopLink->networkLink()->quadraturePoints().size();++q)
//                                        {
//                                            const auto& qPoint(loopLink->networkLink()->quadraturePoints()[q]);
//                                            qPointSlip[q].first -=sourceLoop->windingNumber(qPoint.r + eps*outDir)*sourceLoop->burgers(); // slip vector is negative the burgers vector
//                                            qPointSlip[q].second-=sourceLoop->windingNumber(qPoint.r - eps*outDir)*sourceLoop->burgers(); // slip vector is negative the burgers vector
//
//                                        }
//                                    }
//                                }
//                            }
//
//                            for(size_t q=0;q<loopLink->networkLink()->quadraturePoints().size();++q)
//                            {
//                                if((qPointSlip[q].first-qPointSlip[q].second).squaredNorm()>FLT_EPSILON)
//                                {
//                                    auto& qPoint(loopLink->networkLink()->quadraturePoints()[q]);
//                                    if(qPoint.inclusionID<0)
//                                    {// qPoint is not inside an inclusion, we use the matrix gamma-surface
//
//                                        const double gamma1(fieldLoop->slipSystem()->misfitEnergy(qPointSlip[q].first));  // outer point
//                                        const double gamma2(fieldLoop->slipSystem()->misfitEnergy(qPointSlip[q].second)); // inner point
//
//                                        double gammaNoise(0.0);
//                                        if(fieldLoop->slipSystem()->planeNoise)
//                                        {
//                                            if(loopLink->networkLink()->glidePlanes().size()==1)
//                                            {
//                                                const auto& glidePlane(**loopLink->networkLink()->glidePlanes().begin());
//                                                gammaNoise=std::get<2>(fieldLoop->slipSystem()->gridInterp(qPoint.r-glidePlane.P));
//                                            }
//                                        }
//                                        qPoint.stackingFaultForce+= -(gamma2-gamma1+gammaNoise)*outDir;
//                                    }
//                                    else
//                                    {// qPoint is inside an inclusion, we use the inclusion gamma-surface
//
//                                        const auto& secondPhase(eshelbyInclusions().at(qPoint.inclusionID)->secondPhase);
//                                        if(secondPhase)
//                                        {
//                                                const double gamma1(secondPhase->misfitEnergy(qPointSlip[q].first ,fieldLoop->slipSystem()));  // outer point
//                                                const double gamma2(secondPhase->misfitEnergy(qPointSlip[q].second,fieldLoop->slipSystem())); // inner point
//                                            qPoint.stackingFaultForce+= -(gamma2-gamma1)*outDir;
//                                        }
//                                    }
//                                }
//                            }
//                        }
//                    }
//                }
//            }
        }
        std::cout << magentaColor << std::setprecision(3) << std::scientific << " [" << (std::chrono::duration<double>(std::chrono::system_clock::now() - t1)).count() << " sec]." << defaultColor << std::endl;
        
        
        const auto t2 = std::chrono::system_clock::now();
        std::cout <<"Computing stress field at quadrature points (" << nThreads << " threads) " << std::flush;
#ifdef _OPENMP
#pragma omp parallel for
        for (size_t k = 0; k < this->networkLinks().size(); ++k)
        {
            auto linkIter(this->networkLinks().begin());
            std::advance(linkIter, k);
            const int velGroup(simulationParameters.useSubCycling ? linkIter->second.lock()->velocityGroup(maxVelocity, simulationParameters.subcyclingBins) : 1);
            
            if ((runID % velGroup) == 0)
            {
                linkIter->second.lock()->assembleGlide(true);
            }
            else
            {
                linkIter->second.lock()->assembleGlide(false);
            }
        }
#else
        for (auto &linkIter : this->networkLinks())
        {
            const int velGroup(simulationParameters.useSubCycling ? linkIter.second.lock()->velocityGroup(maxVelocity, simulationParameters.subcyclingBins) : 1);
            
            if ((runID % velGroup) == 0)
            {
                linkIter.second.lock()->assembleGlide(true);
            }
            else
            {
                linkIter.second.lock()->assembleGlide(false);
            }
        }
#endif
        std::cout << magentaColor << std::setprecision(3) << std::scientific << " [" << (std::chrono::duration<double>(std::chrono::system_clock::now() - t2)).count() << " sec]." << defaultColor << std::endl;
    }
    else
    { // For curved segments use quandrature integration of stress field
        assert(0 && "ALL THIS MUST BE RE-IMPLEMENTED FOR CURVED SEGMENTS");
    }
    
    //! -3 Loop over DislocationSubNetworks, assemble subnetwork stiffness matrix and force vector, and solve
    //        std::cout <<"Assembling and solving " << std::flush;
//    DislocationGlideSolver<LoopNetworkType>(*this).solve(runID);
    if(simulationParameters.useSubCycling)
    {
        std::cout <<"Velocity bins for segments " << velocityBinMap.size() << std::endl;
        for (const auto &vBins : velocityBinMap)
        {
            std::cout << magentaColor << vBins.first << " " << vBins.second << ", " << defaultColor << std::flush;
        }
        std::cout << std::endl;
    }
}

template <int dim, short unsigned int corder>
void DislocationNetwork<dim,corder>::updateBoundaryNodes()
{
    
    /*!Step 1. Before removing populate the junction information
     *       This populates the network node where the junctions are needed to be preserved
     *       To pupulate junction information we create a map with key=pair(sourceNetNode,sinkNetNode) and val = set(loopIDs passing through nodes)
     */
    std::map<std::pair<std::shared_ptr<NetworkNodeType>,std::shared_ptr<NetworkNodeType>>,std::set<size_t>> networkNodeLoopMap; //Size_t corresponds to the loopID
    for (const auto& ln : this->loopNodes())
    {
        const auto sharedLNptr(ln.second.lock());
        if (sharedLNptr->periodicPlaneEdge.first)
        {//Node on a boundary: possible junction is between sharedLNptr->periodicPrev() and sharedLNptr->periodicNext()
            if (sharedLNptr->networkNode->loopNodes().size()>1)
            {//junction node
                const auto loopsThis (sharedLNptr->networkNode->loopIDs());
                
                const LoopNodeType *pPrev(sharedLNptr->periodicPrev());
                const LoopNodeType *pNext(sharedLNptr->periodicNext());
                
                
                const auto pPrevNetwork (pPrev->networkNode);
                const auto pNextNetwork (pNext->networkNode);
                
                assert(pPrevNetwork!=nullptr);
                assert(pNextNetwork!=nullptr);
                
                const auto loopspPrev(pPrevNetwork->loopIDs());
                const auto loopspNext(pNextNetwork->loopIDs());
                
                std::set<size_t> tempPrev;
                std::set<size_t> tempNext;
                std::set_intersection(loopspPrev.begin(), loopspPrev.end(), loopsThis.begin(), loopsThis.end(), std::inserter(tempPrev, tempPrev.begin()));
                std::set_intersection(loopsThis.begin(), loopsThis.end(), loopspNext.begin(), loopspNext.end(), std::inserter(tempNext, tempNext.begin()));
                
                if (tempPrev!=tempNext)
                {
                    std::cout<<"For bnd network node"<<sharedLNptr->networkNode->sID<<" loops are "<<std::flush;
                    for (const auto& loop : loopsThis)
                    {
                        std::cout<<loop<<", ";
                    }
                    std::cout<<std::endl;
                    
                    std::cout<<"For prev network node"<<pPrevNetwork->sID<<" loops are "<<std::flush;
                    for (const auto& loop : loopspPrev)
                    {
                        std::cout<<loop<<", ";
                    }
                    std::cout<<std::endl;
                    
                    std::cout<<"For next network node"<<pNextNetwork->sID<<" loops are "<<std::flush;
                    for (const auto& loop : loopspNext)
                    {
                        std::cout<<loop<<", ";
                    }
                    std::cout<<std::endl;
                    throw std::runtime_error("BND node must have the same loops as the common loops between the internal nodes");
                }
                else
                {
                    if (pPrevNetwork->sID < pNextNetwork->sID)
                    {
                        networkNodeLoopMap.emplace(std::make_pair(pPrevNetwork, pNextNetwork), tempPrev);
                    }
                    else
                    {
                        networkNodeLoopMap.emplace(std::make_pair(pNextNetwork, pPrevNetwork), tempPrev);
                    }
                }
            }
        }
        else
        {//Node not on a boundary: possible junction is between sharedLNptr and sharedLNptr->periodicNext()
            if (sharedLNptr->networkNode->loopNodes().size()>1)
            {// a junction node
                if (sharedLNptr->boundaryNext().size()==0 && sharedLNptr->periodicNext()->networkNode->loopNodes().size()>1)
                {
                    if (sharedLNptr->periodicPlanePatch()!=sharedLNptr->periodicNext()->periodicPlanePatch())
                    {// Junction is across boundary
                        const auto netLink (sharedLNptr->next.second->networkLink());
                        if(netLink)
                        {
                            std::set<size_t> netLinkLoopIDs (netLink->loopIDs());
                            if (netLink->loopLinks().size()>=2)
                            {
                                //a junction node moving out
                                if (sharedLNptr->networkNode->sID < sharedLNptr->periodicNext()->networkNode->sID)
                                {
                                    networkNodeLoopMap.emplace(std::make_pair(sharedLNptr->networkNode, sharedLNptr->periodicNext()->networkNode), netLinkLoopIDs);
                                }
                                else
                                {
                                    networkNodeLoopMap.emplace(std::make_pair(sharedLNptr->periodicNext()->networkNode,sharedLNptr->networkNode), netLinkLoopIDs);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
    /*!Step 2. Remove nodes in danglingBoundaryLoopNodes
     */
    if (danglingBoundaryLoopNodes.size())
    {
 //       std::cout << "Removing bnd Nodes" << std::endl;
        VerboseDislocationNetwork(1, "Removing bnd Nodes"<< std::endl;);

        for (const auto &node : danglingBoundaryLoopNodes)
        {
            this->removeLoopNode(node->sID);
        }
    }
    
//    std::cout << "Inserting new boundary nodes" << std::endl;
    VerboseDislocationNetwork(1, "Inserting new boundary nodes"<< std::endl;);

    /*!Step 3. Create a map of new boundary NetworkNode to be inserted
     *       key = tuple<souceNetNode,
     *                sinkNetNode,
     *                set<all mesh faces that current link is crossing>,
     *                set<current mesh faces contining new boundary NetworkNode>,
     *                number of mesh intersection crossed by link beyond  new boundary NetworkNode>
     *       value= new  boundary NetworkNode
     */
    std::map<std::tuple<const std::shared_ptr<NetworkNodeType>, const std::shared_ptr<NetworkNodeType>, std::set<const PlanarMeshFace<dim> *>, std::set<const PlanarMeshFace<dim> *>,size_t>, std::shared_ptr<NetworkNodeType>> newNetworkNodesMap;
    for (const auto &weakLoop : this->loops())
    {
        
        const auto loop(weakLoop.second.lock());
        if (loop->periodicGlidePlane)
        {
            VerboseDislocationNetwork(1, " DisloationNetwork::DislocationLoop " << loop->sID << " Updating boundary nodes " << std::endl;);
            
            std::vector<std::pair<VectorLowerDim, const LoopNodeType *const>> loopNodesPos;
            std::map<std::tuple<const LoopNodeType *, const LoopNodeType *, const std::pair<const PeriodicPlaneEdge<dim> *, const PeriodicPlaneEdge<dim> *>>, const LoopNodeType *> bndNodesMap;
            for (const auto &loopLink : loop->linkSequence())
            {
                if (!loopLink->source->periodicPlaneEdge.first)
                {
                    loopNodesPos.emplace_back(loop->periodicGlidePlane->referencePlane->localPosition(loopLink->source->get_P()), loopLink->source.get());
                }
                else
                {
                    
                    if (loopLink->source->periodicPlaneEdge.second)
                    {// loopLink->source is a corner boundary node. For nodes intersecting the corner, storage should be from smaller edgeId and larger edgeID
                        if (loopLink->source->periodicPlaneEdge.first->edgeID<loopLink->source->periodicPlaneEdge.second->edgeID)
                        {
                            bndNodesMap.emplace(std::make_tuple(loopLink->source->periodicPrev(), loopLink->source->periodicNext(), std::make_pair(loopLink->source->periodicPlaneEdge.first.get(), loopLink->source->periodicPlaneEdge.second.get())), loopLink->source.get());
                        }
                        else
                        {
                            bndNodesMap.emplace(std::make_tuple(loopLink->source->periodicPrev(), loopLink->source->periodicNext(), std::make_pair(loopLink->source->periodicPlaneEdge.second.get(),loopLink->source->periodicPlaneEdge.first.get())), loopLink->source.get());
                        }
                    }
                    else
                    {
                        bndNodesMap.emplace(std::make_tuple(loopLink->source->periodicPrev(), loopLink->source->periodicNext(), std::make_pair(loopLink->source->periodicPlaneEdge.first.get(), loopLink->source->periodicPlaneEdge.second.get())), loopLink->source.get());
                    }
                }
            }
            
            if (loopNodesPos.size())
            {
                const auto polyInt(loop->periodicGlidePlane->polygonPatchIntersection(loopNodesPos)); // Note: last element of polyInt is always an internal node
                // 2dPos,shift,edgeIDs ,set of edgeiD(all edges crossed),size_t(number of boundary nodes on same loop link after current node),LoopNodeType* (from loopNodesPos if internal, nullptr otherwise)
                //edges still needed to be traversed will be reverse of the value if junction formation includes link in the opposite direction
                std::map<const LoopNodeType *, size_t> polyIntMap; // ID of internal loopNodes into polyInt
                for (size_t p = 0; p < polyInt.size(); ++p)
                {
                    if (std::get<5>(polyInt[p]))
                    {// an internal node
                        polyIntMap.emplace(std::get<5>(polyInt[p]), p);
                    }
                }
                
                for (size_t k = 0; k < loopNodesPos.size(); ++k)
                {
                    const auto periodicPrev(loopNodesPos[k].second);
                    const auto periodicPrevNetwork(periodicPrev->networkNode);
                    const auto polyIter(polyIntMap.find(periodicPrev));
                    assert(polyIter != polyIntMap.end());
                    const size_t p(polyIter->second);
                    
                    const size_t k1(k < loopNodesPos.size() - 1 ? k + 1 : 0);
                    const auto periodicNext(loopNodesPos[k1].second);
                    const auto periodicNextNetwork(periodicNext->networkNode);
                    const auto polyIter1(polyIntMap.find(periodicNext));
                    assert(polyIter1 != polyIntMap.end());
                    const size_t p1(polyIter1->second);
                    
                    const auto periodicNetworkSource(periodicPrevNetwork->sID < periodicNextNetwork->sID ? periodicPrevNetwork : periodicNextNetwork);
                    const auto periodicNetworkSink(periodicPrevNetwork->sID < periodicNextNetwork->sID ? periodicNextNetwork : periodicPrevNetwork);
                    
                    VerboseDislocationNetwork(2, " PeriodicNetworkSource " << periodicNetworkSource->sID << std::endl;);
                    VerboseDislocationNetwork(2, " PeriodicNetworkSink " << periodicNetworkSink->sID << std::endl;);
                    // std::cout<<"periodicPrev->periodicNext"<<periodicPrev->tag()<<"==>"<<periodicNext->tag()<<std::endl;
                    
                    const LoopNodeType *currentSource(periodicPrev);
                    size_t p2 = (p + 1) % polyInt.size();
                    while (p2 < p1)
                    {
                        const auto periodicPatch(loop->periodicGlidePlane->getPatch(std::get<1>(polyInt[p2])));
                        const auto periodicPatchEdge(std::get<2>(polyInt[p2]).second < 0 ? std::make_pair(periodicPatch->edges()[std::get<2>(polyInt[p2]).first], nullptr) : std::make_pair(periodicPatch->edges()[std::get<2>(polyInt[p2]).first], periodicPatch->edges()[std::get<2>(polyInt[p2]).second]));
                        
                        VerboseDislocationNetwork(2, " First EdgeID on periodicPatchEdge " << periodicPatchEdge.first->edgeID << std::endl;);
                        if (periodicPatchEdge.second)
                        {
                            VerboseDislocationNetwork(2, " Second EdgeID on periodicPatchEdge " << periodicPatchEdge.second->edgeID << std::endl;);
                        }
                        
                        //                            typename std::map<std::tuple<const LoopNodeType *, const LoopNodeType *, const std::pair<const PeriodicPlaneEdge<dim> *, const PeriodicPlaneEdge<dim> *>>, const LoopNodeType *>::iterator bndIter(bndNodesMap.end());
                        auto bndIter(bndNodesMap.end());
                        
                        //For the case where the second edge exists, the pair of edge should be from minimum to maximum.
                        if (periodicPatchEdge.second)
                        {
                            if (periodicPatchEdge.first->edgeID < periodicPatchEdge.second->edgeID)
                            {
                                bndIter = bndNodesMap.find(std::make_tuple(periodicPrev, periodicNext, std::make_pair(periodicPatchEdge.first.get(), periodicPatchEdge.second.get())));
                            }
                            else
                            {
                                bndIter = bndNodesMap.find(std::make_tuple(periodicPrev, periodicNext, std::make_pair(periodicPatchEdge.second.get(), periodicPatchEdge.first.get())));
                            }
                        }
                        else
                        {
                            bndIter = bndNodesMap.find(std::make_tuple(periodicPrev, periodicNext, std::make_pair(periodicPatchEdge.first.get(), periodicPatchEdge.second.get())));
                        }
                        
                        if (bndIter != bndNodesMap.end())
                        { // exising bnd node found
                            // std::cout<<" Using boundary node "<<bndIter->second->tag()<<std::endl;
                            currentSource = bndIter->second;
                        }
                        else
                        {
                            // Original
                            
                            // const VectorDim loopNodePos(loop->periodicGlidePlane->referencePlane->globalPosition(std::get<0>(polyInt[p2])));
                            // const VectorDim networkNodePos(loopNodePos + std::get<1>(polyInt[p2]));
                            
                            const auto periodicPatchEdgesAll(std::get<3>(polyInt[p2])); // Get all the edges which the nodes are intesecting
                            // std::cout<<"periodicPatchEdgesAll size is "<<periodicPatchEdgesAll.size()<<std::endl;
                            //get the mesh faces corresponding to the edges from periodicPatchEdgesAll
                            std::set<const PlanarMeshFace<dim> *> allMeshFaces;
                            size_t totalNumEdgesCrossed(0);
                            for (const auto& allPatchMap : periodicPatchEdgesAll)
                            {
                                // std::cout<<"periodicPatchEdgesAll second size is "<<allPatchMap.second.size()<<std::endl;
                                totalNumEdgesCrossed+=allPatchMap.second.size();
                                const auto periodicPatchTemp(loop->periodicGlidePlane->getPatch(allPatchMap.first)); //From shift grab the patch
                                
                                for (const auto& allPatchEdges : allPatchMap.second)
                                {
                                    const auto periodicPatchEdgeTemp(allPatchEdges.second < 0 ? std::make_pair(periodicPatchTemp->edges()[allPatchEdges.first], nullptr)
                                                                     : std::make_pair(periodicPatchTemp->edges()[allPatchEdges.first], periodicPatchTemp->edges()[allPatchEdges.second]));
                                    
                                    for (const auto &pmface : periodicPatchEdgeTemp.first->meshIntersection->faces)
                                    {
                                        allMeshFaces.emplace(pmface);
                                    }
                                    if (periodicPatchEdgeTemp.second)
                                    {
                                        for (const auto &pmface : periodicPatchEdgeTemp.second->meshIntersection->faces)
                                        {
                                            allMeshFaces.emplace(pmface);
                                        }
                                    }
                                }
                            }
                            
                            
                            
                            
                            const auto networkLoopMapIter(networkNodeLoopMap.find(std::make_pair(periodicNetworkSource, periodicNetworkSink)));
                            // std::set<LoopType *> commonLoops;
                            bool loopBelongtoCommonLoop(false);
                            //                                bool firstLoopInJunction(false); //only for the first loop we need to insert a new network node,otherwise we should be able to find an already existent network node
                            if (networkLoopMapIter != networkNodeLoopMap.end())
                            {
                                // Insert the loops as determined from the map to preserve the junction information
                                const auto networkLoopSetIter (networkLoopMapIter->second.find(loop->sID));
                                loopBelongtoCommonLoop = (networkLoopSetIter != networkLoopMapIter->second.end());
                                //                                    firstLoopInJunction = (networkLoopSetIter == networkLoopMapIter->second.begin()); //If this is true then we need to create a new network node
                            }
                            else
                            {
                                // Insert the common loop
                                loopBelongtoCommonLoop = false;
                            }
                            
                            const VectorDim loopNodePostemp(loop->periodicGlidePlane->referencePlane->globalPosition(std::get<0>(polyInt[p2])));
                            
                            VectorDim networkNodePos(periodicPatchEdge.first->meshIntersection->snap(loopNodePostemp + std::get<1>(polyInt[p2])));
                            //                                std::set<std::shared_ptr<PeriodicPlanePatch<dim>>> auxiliaryPatches; //Aux patches with only be populated if the second patch edge exists
                            // i.e. The interseection is taking place diagonally
                            if (periodicPatchEdge.second)
                            {
                                SegmentSegmentDistance<dim> ssd(periodicPatchEdge.first->meshIntersection->P0, periodicPatchEdge.first->meshIntersection->P1,
                                                                periodicPatchEdge.second->meshIntersection->P0, periodicPatchEdge.second->meshIntersection->P1);
                                assert(ssd.dMin < FLT_EPSILON && "Two edges must intersect");
                                networkNodePos = 0.5 * (ssd.x0 + ssd.x1);
                                
                                //                                    auxiliaryPatches.insert(loop->periodicGlidePlane->getPatch(periodicPatch->shift+periodicPatchEdge.first->deltaShift));
                                //                                    auxiliaryPatches.insert(loop->periodicGlidePlane->getPatch(periodicPatch->shift+periodicPatchEdge.second->deltaShift));
                                
                                //                                    if ((networkNodePos - (loopNodePostemp + std::get<1>(polyInt[p2]))).norm() > FLT_EPSILON)
                                //                                    {
                                //                                        std::cout<<" PeriodicNetworkSource "<<periodicNetworkSource->sID<<std::endl;
                                //                                        std::cout<<" PeriodicNetworkSink "<<periodicNetworkSink->sID<<std::endl;
                                //                                        std::cout<<" First edge mesh intersection "<<std::endl;
                                //                                        for (const auto& mf : periodicPatchEdge.first->meshIntersection->faces)
                                //                                        {
                                //                                            std::cout<<mf->sID<<" with outnormal"<<mf->outNormal().transpose()<< "\t "<<std::flush;
                                //                                        }
                                //                                        std::cout<<std::endl;
                                //                                        std::cout<<" Second edge mesh intersection "<<std::endl;
                                //                                        for (const auto& mf : periodicPatchEdge.second->meshIntersection->faces)
                                //                                        {
                                //                                            std::cout<<mf->sID<<" with outnormal"<<mf->outNormal().transpose()<< "\t "<<std::flush;
                                //                                        }
                                //                                        std::cout<<std::endl;
                                //                                        std::cout<<" Total edges crossed "<<totalNumEdgesCrossed<<" numEdges yet to be crossed "<<(std::get<4>(polyInt[p2]))<<std::endl;
                                //                                        std::cout<<" Trying to set network Node position to "<<networkNodePos.transpose()<<std::endl;
                                //                                        std::cout<<" Actual network node position  "<<(loopNodePostemp + std::get<1>(polyInt[p2])).transpose()<<std::endl;
                                //                                    }
                                assert((networkNodePos - (loopNodePostemp + std::get<1>(polyInt[p2]))).norm() < FLT_EPSILON && "Position mismatch");
                            }
                            //                                else
                            //                                {
                            //                                    auxiliaryPatches.insert(std::shared_ptr<PeriodicPlanePatch<dim>>{nullptr});
                            //                                    networkNodePos = periodicPatchEdge.first->meshIntersection->snap(loopNodePostemp + std::get<1>(polyInt[p2]));
                            //                                }
                            const VectorDim loopNodePos(networkNodePos - std::get<1>(polyInt[p2]));
                            
                            const auto currentLoopLink(currentSource->next.second);
                            
                            const auto currentNetworkLink(currentLoopLink->networkLink());
                            
                            std::set<const PlanarMeshFace<dim> *> tmpMeshFaces;
                            for (const auto &pmface : periodicPatchEdge.first->meshIntersection->faces)
                            {
                                tmpMeshFaces.emplace(pmface);
                            }
                            if (periodicPatchEdge.second)
                            {
                                for (const auto &pmface : periodicPatchEdge.second->meshIntersection->faces)
                                {
                                    tmpMeshFaces.emplace(pmface);
                                }
                            }
                            
                            const size_t edgesStillRemainingtoCross(std::get<4>(polyInt[p2]));
                            VerboseDislocationNetwork(2, " edgesStillRemainingtoCross " << edgesStillRemainingtoCross << std::endl;);
                            VerboseDislocationNetwork(2, " total edges crossed " << totalNumEdgesCrossed << std::endl;);
                            
                            const int u(std::min(edgesStillRemainingtoCross,totalNumEdgesCrossed - edgesStillRemainingtoCross - 1));
                            //                                size_t u(0);
                            //
                            //                                if (loopBelongtoCommonLoop)
                            //                                {
                            //                                    if (periodicNetworkSource == periodicNetworkSink)
                            //                                    {
                            //                                        if (firstLoopInJunction)
                            //                                        {
                            //                                            //We can go as it is...i.e. from periodicPrev to periodicNext
                            //                                            u = edgesStillRemainingtoCross;
                            //                                        }
                            //                                        else
                            //                                        {
                            //                                            //need to check the alignment with the first loop junction nodes
                            //                                            const VectorDim currentLoopLinkChord(periodicNext->get_P()-periodicPrev->get_P());
                            //                                            const double currentLoopLinkChordLength(currentLoopLinkChord.norm());
                            //                                            VectorDim firstLoopLinkChord (VectorDim::Zero());
                            //                                            for (const auto& ln : periodicNetworkSource->loopNodes())
                            //                                            {
                            //                                                if (ln->loop()->sID==*(networkLoopMapIter->second.begin()))
                            //                                                {
                            //                                                    //Can grab any oriented direction
                            //                                                    if (ln->periodicNext()->networkNode==periodicNetworkSource)
                            //                                                    {
                            //                                                        //We are at the source of the loop link
                            //                                                        // std::cout<<" Case A "<<ln->tag()<<ln->get_P().transpose()<<std::endl;
                            //                                                        // std::cout<<" Case A "<<ln->periodicNext()->tag()<<ln->periodicNext()->get_P().transpose()<<std::endl;
                            //                                                        firstLoopLinkChord=ln->periodicNext()->get_P()-ln->get_P();
                            //                                                    }
                            //                                                    else if (ln->periodicPrev()->networkNode==periodicNetworkSource)
                            //                                                    {
                            //                                                        //We are at the sink of the loop link
                            //                                                        // std::cout<<" Case B "<<ln->tag()<<ln->get_P().transpose()<<std::endl;
                            //                                                        // std::cout<<" Case B "<<ln->periodicPrev()->tag()<<ln->periodicPrev()->get_P().transpose()<<std::endl;
                            //                                                        firstLoopLinkChord=ln->get_P()-ln->periodicPrev()->get_P();
                            //                                                    }
                            //                                                    // if (ln->periodicNext()->networkNode==periodicNetworkSource)
                            //                                                    // assert(false && "FINISH HERE");
                            //                                                }
                            //                                            }
                            //                                            const double firstLoopLinkChordLength(firstLoopLinkChord.norm());
                            //                                            assert(firstLoopLinkChordLength>FLT_EPSILON && "First looplink chord must be finite length");
                            //                                            assert(fabs(firstLoopLinkChordLength - currentLoopLinkChordLength)<FLT_EPSILON && "Chord length must be same");
                            //                                            // Determine the alignment
                            //                                            const VectorDim currentLoopLinkDir(currentLoopLinkChord / currentLoopLinkChordLength);
                            //                                            const VectorDim firstLoopLinkDir(firstLoopLinkChord / firstLoopLinkChordLength);
                            //                                            const double dirRef(currentLoopLinkDir.dot(firstLoopLinkDir));
                            //                                            if (fabs(dirRef + 1.0) < FLT_EPSILON)
                            //                                            {
                            //                                                //anti aligned
                            //                                                u = totalNumEdgesCrossed - edgesStillRemainingtoCross - 1; // 1 is to make it to go to 0 (at the end 0 edges should be crossed)
                            //                                            }
                            //                                            else
                            //                                            {
                            //                                                //aligned
                            //                                                u = edgesStillRemainingtoCross;
                            //                                            }
                            //                                        }
                            ////                                        assert(false && "First case occuring... check how the implementation works here");
                            //                                    }
                            //                                    else
                            //                                    {
                            //                                        if (periodicPrevNetwork == periodicNetworkSource)
                            //                                        {
                            //                                            // This condition will give problem if perioicNetworkSource and Sink are same
                            //                                            u = edgesStillRemainingtoCross;
                            //                                        }
                            //                                        else if (periodicNextNetwork == periodicNetworkSource)
                            //                                        {
                            //                                            // This condition will give problem if perioicNetworkSource and Sink are same
                            //                                            u = totalNumEdgesCrossed - edgesStillRemainingtoCross - 1; // 1 is to make it to go to 0 (at the end 0 edges should be crossed)
                            //                                        }
                            //                                    }
                            //                                }
                            
                            // std::cout<<"Total number of edges crossed "<<totalNumEdgesCrossed<<std::endl;
                            const auto key(std::make_tuple(periodicNetworkSource, periodicNetworkSink, allMeshFaces, tmpMeshFaces,u));
                            // std::cout<<" Key is "<<periodicNetworkSource->tag()<<" "<<periodicNetworkSink->tag()<<" "<<u<<" "<<std::endl;
                            // std::cout<<"All Meshfaces "<<std::endl;
                            // for (const auto& mf : allMeshFaces)
                            // {
                            //     std::cout<<mf->sID<<" with outnormal"<<mf->outNormal().transpose()<< "\t "<<std::flush;
                            // }
                            // std::cout<<std::endl<< " Temp mesh faces "<<std::endl;
                            // for (const auto& mf : tmpMeshFaces)
                            // {
                            //     std::cout<<mf->sID<<" with outnormal"<<mf->outNormal().transpose()<< "\t "<<std::flush;
                            // }
                            // std::cout<<std::endl;
                            if (loopBelongtoCommonLoop)
                            {
                                //Either we need to insert a new node or we grab an already existent node
                                const auto networkNodeIter(newNetworkNodesMap.find(key));
                                if (networkNodeIter == newNetworkNodesMap.end())
                                {
                                    // Insert a new node
                                    const auto newNetNode(this->networkNodes().create(networkNodePos, VectorDim::Zero(), 1.0)); // TODO compute velocity and velocityReduction by interpolation
                                    // std::cout << "emplacing " << currentNetworkLink->tag() << "@" << std::setprecision(15) << std::scientific  << ", newNetNode=" << newNetNode->tag() << std::endl;
                                    VerboseDislocationNetwork(2, " Junction case....Inserting a new node " << newNetNode->tag() << std::endl;);
                                    
                                    newNetworkNodesMap.emplace(key, newNetNode);
                                    const auto newLoopNode(this->loopNodes().create(loop, newNetNode, loopNodePos, periodicPatch, periodicPatchEdge));
                                    currentSource = this->expandLoopLink(*currentLoopLink, newLoopNode).get();
                                }
                                else
                                {
                                    // A node has already been inserted... use that node
                                    const VectorDim commensurateLoopNodePos(networkNodeIter->second->get_P() - std::get<1>(polyInt[p2]));
                                    const auto newLoopNode(this->loopNodes().create(loop, networkNodeIter->second, commensurateLoopNodePos, periodicPatch, periodicPatchEdge));
                                    currentSource = this->expandLoopLink(*currentLoopLink, newLoopNode).get();
                                }
                                
                                
                                //                                    if (firstLoopInJunction)
                                //                                    {
                                //                                        //May be a self annihilation case...For the self annihilation case use the already existent key
                                //                                        const auto networkNodeIter(newNetworkNodesMap.find(key));
                                //                                        if (networkNodeIter == newNetworkNodesMap.end())
                                //                                        {
                                //                                            // Insert a new node
                                //                                            const auto newNetNode(this->networkNodes().create(networkNodePos, VectorDim::Zero(), 1.0)); // TODO compute velocity and velocityReduction by interpolation
                                //                                            // std::cout << "emplacing " << currentNetworkLink->tag() << "@" << std::setprecision(15) << std::scientific  << ", newNetNode=" << newNetNode->tag() << std::endl;
                                //                                            VerboseDislocationNetwork(2, " Junction case....Inserting a new node " << newNetNode->tag() << std::endl;);
                                //
                                //                                            newNetworkNodesMap.emplace(key, newNetNode);
                                //                                            const auto newLoopNode(this->loopNodes().create(loop, newNetNode, loopNodePos, periodicPatch, periodicPatchEdge));
                                //                                            currentSource = this->expandLoopLink(*currentLoopLink, newLoopNode).get();
                                //                                        }
                                //                                        else
                                //                                        {
                                //                                            // A node has already been inserted... use that node
                                //                                            const VectorDim commensurateLoopNodePos(networkNodeIter->second->get_P() - std::get<1>(polyInt[p2]));
                                //                                            const auto newLoopNode(this->loopNodes().create(loop, networkNodeIter->second, commensurateLoopNodePos, periodicPatch, periodicPatchEdge));
                                //                                            currentSource = this->expandLoopLink(*currentLoopLink, newLoopNode).get();
                                //                                        }
                                //                                    }
                                //                                    else
                                //                                    {
                                //                                        //grab an already existent node
                                //                                        const auto networkNodeIter(newNetworkNodesMap.find(key));
                                //                                        if (networkNodeIter == newNetworkNodesMap.end())
                                //                                        {
                                //                                            std::cout<<"Printing common loops "<<std::endl;
                                //                                            for (const auto& loop : networkLoopMapIter->second)
                                //                                            {
                                //                                                std::cout<<loop<<std::endl;
                                //                                            }
                                //                                            std::cout<<"Printing already existent keys "<< newNetworkNodesMap.size()<<" new network node map size"<<std::endl;
                                //                                            for (const auto& nnMap : newNetworkNodesMap)
                                //                                            {
                                //                                                std::cout<<std::get<0>(nnMap.first)->sID<<" "<<std::get<1>(nnMap.first)->sID<<" All Mesh face "<<std::endl;
                                //                                                std::cout << "All Meshfaces " << std::endl;
                                //                                                for (const auto &mf : std::get<2>(nnMap.first))
                                //                                                {
                                //                                                    std::cout << mf->sID << " with outnormal" << mf->outNormal().transpose() << "\t " << std::flush;
                                //                                                }
                                //                                                std::cout<<std::endl<< " Temp mesh faces "<<std::endl;
                                //                                                for (const auto &mf : std::get<3>(nnMap.first))
                                //                                                {
                                //                                                    std::cout << mf->sID << " with outnormal" << mf->outNormal().transpose() << "\t " << std::flush;
                                //                                                }
                                //                                                std::cout << std::endl;
                                //                                                std::cout<<" NetworkNode is "<<nnMap.second->sID<<std::endl;
                                //                                            }
                                //
                                //                                            std::cout<<"Printing bnd Nodes Map "<< bndNodesMap.size()<<" bndNodemap size"<<std::endl;
                                //                                            // std::map<std::tuple<const LoopNodeType *, const LoopNodeType *, const std::pair<const PeriodicPlaneEdge<dim> *, const PeriodicPlaneEdge<dim> *>>, const LoopNodeType *> bndNodesMap;
                                //                                            for (const auto& bndNode : bndNodesMap )
                                //                                            {
                                //                                                std::cout<<"Prev Node "<<(std::get<0>(bndNode.first))->tag()<<std::endl;
                                //                                                std::cout<<"Next Node "<<(std::get<1>(bndNode.first))->tag()<<std::endl;
                                //                                                std::cout<<"First Edge ID "<<(std::get<2>(bndNode.first)).first->edgeID<<std::endl;
                                //                                                if ((std::get<2>(bndNode.first)).second)
                                //                                                    std::cout<<"Second Edge ID "<<(std::get<2>(bndNode.first)).second->edgeID<<std::endl;
                                //                                                std::cout<<" Stored Node "<<bndNode.second->tag()<<std::endl;
                                //
                                //                                                // std::cout<<(std::get<0>(bndNode.first))->tag()<<" "<<(std::get<1>(bndNode.first))->tag()<<" "<<(std::get<2>(bndNode.first)).first->edgeID<<" "<<(std::get<2>(bndNode.first)).second->edgeID<<" "<<
                                //                                                // bndNode.second->tag()<<std::endl;
                                //                                            }
                                //
                                //                                        }
                                //                                        assert(networkNodeIter != newNetworkNodesMap.end() && "Inserting bnd node corresponding to a junction... bnd node should be present already");
                                //                                        /* The network node position must be commensurate */
                                //                                        if ((networkNodeIter->second->get_P() - networkNodePos).norm() > 10000 * FLT_EPSILON) // THis condition is just to check for widely different positions
                                //                                        {
                                //                                            std::cout << std::scientific << std::setprecision(15) << " PeriodicNetwork Source " << periodicNetworkSource->sID << "P= " << periodicNetworkSource->get_P().transpose() << "\n PeriodicNetwork Sink " << periodicNetworkSink->sID << " P = " << periodicNetworkSink->get_P().transpose() << std::endl;
                                //                                            std::cout << " Position of the network node " << std::scientific << std::setprecision(15) << networkNodeIter->second->get_P().transpose() << std::endl;
                                //                                            std::cout << " Actual Position of the network node should be " << std::scientific << std::setprecision(15) << networkNodePos.transpose() << std::endl;
                                //                                            std::cout << " Position difference " << std::setprecision(15) << (networkNodeIter->second->get_P() - networkNodePos).transpose() << std::endl;
                                //                                            std::cout << " Position difference norm " << std::setprecision(15) << (networkNodeIter->second->get_P() - networkNodePos).norm() << std::endl;
                                //                                            assert(false && "Loop node position and network node position mismatch");
                                //                                        }
                                //                                        // Here we want the loop node position to be commensurate with the network node (This is important for minimizing accumulating error)
                                //                                        VerboseDislocationNetwork(2, " Junction case....Grabbing an already existent node " << networkNodeIter->second->tag() << std::endl;);
                                //                                        const VectorDim commensurateLoopNodePos(networkNodeIter->second->get_P() - std::get<1>(polyInt[p2]));
                                //                                        const auto newLoopNode(this->loopNodes().create(loop, networkNodeIter->second, commensurateLoopNodePos, periodicPatch, periodicPatchEdge));
                                //                                        currentSource = this->expandLoopLink(*currentLoopLink, newLoopNode).get();
                                //                                    }
                            }
                            else
                            {
                                //create a new node
                                const auto newNetNode(this->networkNodes().create(networkNodePos, VectorDim::Zero(), 1.0)); // TODO compute velocity and velocityReduction by interpolation
                                VerboseDislocationNetwork(2, " non-Junction case....Inserting a new node " << newNetNode->tag() << std::endl;);
                                
                                // std::cout << "emplacing " << currentNetworkLink->tag() << "@" << std::setprecision(15) << std::scientific  << ", newNetNode=" << newNetNode->tag() << std::endl;
                                newNetworkNodesMap.emplace(key, newNetNode);
                                const auto newLoopNode(this->loopNodes().create(loop, newNetNode, loopNodePos, periodicPatch, periodicPatchEdge));
                                currentSource = this->expandLoopLink(*currentLoopLink, newLoopNode).get();
                            }
                        }
                        p2 = (p2 + 1) % polyInt.size();
                    }
                }
            }
        }
    }
    
    danglingBoundaryLoopNodes.clear();
    
}


template <int dim, short unsigned int corder>
void DislocationNetwork<dim,corder>::moveGlide(const double & dt_in)
{/*! Moves all nodes in the DislocationNetwork using the stored velocity and current dt
  */
    
    const auto t0= std::chrono::system_clock::now();
    std::cout<<"Moving DislocationNodes by glide (dt="<<dt_in<< ")... "<<std::flush;
    
    if(simulationParameters.isPeriodicSimulation())
    {
        
        
        danglingBoundaryLoopNodes.clear();
        
        for(auto& node : this->networkNodes())
        {// Expansion
            // std::cout<<"Trying set p for "<<node.second.lock()->sID<<" LoopNdoe size ==>"<<node.second.lock()->loopNodes().size()<<" GlidePlane size ==>"
            // <<node.second.lock()->glidePlanes().size()<<std::endl;
            // std::cout<<std::setprecision(15)<<" Old Position "<<node.second.lock()->get_P().transpose()<<std::endl;
            // std::cout<<std::setprecision(15)<<" New Position "<<(node.second.lock()->get_P()+node.second.lock()->get_V()*dt_in).transpose()<<std::endl;
            // std::cout<<"Checking if the loop contains the position "<<std::endl;
            //                for (const auto& loopN : node.second.lock()->loopNodes())
            //                {
            //                    const VectorDim patchShift(loopN->periodicPlanePatch()? loopN->periodicPlanePatch()->shift : VectorDim::Zero());
            //                    const VectorDim oldPosition(loopN->get_P());
            //                    const VectorDim oldPositionN(node.second.lock()->get_P()-patchShift);
            //                    const VectorDim newPosition(node.second.lock()->get_P()+node.second.lock()->get_V()*dt_in-patchShift);
            //
            //                    // std::cout<<" Loop "<<loopN->loop()->sID<<" contains position OLD Position"<<loopN->loop()->glidePlane->contains(oldPosition)<<std::endl;
            //                    // std::cout<<" Loop "<<loopN->loop()->sID<<" contains position New Position"<<loopN->loop()->glidePlane->contains(newPosition)<<std::endl;
            //                }
            node.second.lock()->trySet_P(node.second.lock()->get_P()+node.second.lock()->get_V()*dt_in);
        }
        std::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<defaultColor<<std::endl;
        
        updateBoundaryNodes();
        
        
    }
    else
    {
        std::cout<<"        Moving DislocationNodes by glide (dt="<<dt_in<< ")... "<<std::flush;
        const auto t0= std::chrono::system_clock::now();
        for (auto& nodeIter : this->networkNodes())
        {
            //                nodeIter.second.lock()->moveGlide(dt_in);
            nodeIter.second.lock()->set_P(nodeIter.second.lock()->get_P()+nodeIter.second.lock()->get_V()*dt_in);
        }
        std::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<defaultColor<<std::endl;
    }
    
    
}

/**********************************************************************/
template <int dim, short unsigned int corder>
void DislocationNetwork<dim,corder>::storeSingleGlideStepDiscreteEvents(const long int&)
{
    
    crossSlipMaker.findCrossSlipSegments();
    crossSlipMaker.execute(); // this is now performed before moving, since internally it computes forces and velocities on the new segments
    
}

/**********************************************************************/
template <int dim, short unsigned int corder>
void DislocationNetwork<dim,corder>::executeSingleGlideStepDiscreteEvents(const long int& runID)
{
    
//    crossSlipMaker.execute();
    //        //! 13- Node redistribution
    networkRemesher.remesh(runID);
    //        //! 12- Form Junctions
    junctionsMaker.formJunctions(3.0*networkRemesher.Lmin);
    //Calling remesh again so that any other topological changes created by junctions whihc are otherwise removable can be removed
    //        //! 13- Node redistribution
    //        this->io().output(runID);
    
    networkRemesher.remesh(runID);
    //        updateVirtualBoundaryLoops();
    
}


template <int dim, short unsigned int corder>
int DislocationNetwork<dim,corder>::verboseDislocationNetwork=0;

template class DislocationNetwork<3,0>;
// template class DislocationNetwork<3,1>;

}
#endif
