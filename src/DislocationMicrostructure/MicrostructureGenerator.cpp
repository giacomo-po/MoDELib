/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_MicrostructureGenerator_cpp_
#define model_MicrostructureGenerator_cpp_

#include <fstream>
#include <filesystem>


#include <MicrostructureGenerator.h>
#include <PeriodicDipoleGenerator.h>
#include <PeriodicLoopGenerator.h>
#include <PrismaticLoopGenerator.h>
#include <SphericalInclusionsGenerator.h>
#include <PolyhedronInclusionsGenerator.h>
#include <IrradiationDefectsGenerator.h>

#include <VTKGenerator.h>

namespace model
{

    /**********************************************************************/
    MicrostructureGenerator::MicrostructureGenerator(DislocationDynamicsBase<3>& ddBase_in) :
    /* init*/ configIO(ddBase_in.simulationParameters.traitsIO.evlFolder)
    /* init*/,auxIO(ddBase_in.simulationParameters.traitsIO.auxFolder)
    /* init*/,ddBase(ddBase_in)
    /* init*/,outputBinary(TextFileParser(ddBase.simulationParameters.traitsIO.ddFile).readScalar<int>("outputBinary",true))
    /* init*/,minSize(0.1*std::min(ddBase.mesh.xMax(0)-ddBase.mesh.xMin(0),std::min(ddBase.mesh.xMax(1)-ddBase.mesh.xMin(1),ddBase.mesh.xMax(2)-ddBase.mesh.xMin(2))))
    /* init*/,maxSize(std::max(ddBase.mesh.xMax(0)-ddBase.mesh.xMin(0),std::max(ddBase.mesh.xMax(1)-ddBase.mesh.xMin(1),ddBase.mesh.xMax(2)-ddBase.mesh.xMin(2))))
    {
        
        std::cout<<greenBoldColor<<"Generating microstructure for folder "<<ddBase.simulationParameters.traitsIO.simulationFolder<<defaultColor<<std::endl;
        
        
        // Some sanity checks
        if(ddBase.mesh.volume()<FLT_EPSILON)
        {
            throw std::runtime_error("mesh "+ddBase.simulationParameters.traitsIO.meshFile+" is empty.");
        }
        
        
//        std::ifstream initialMicrostructureFile(ddBase.simulationParameters.traitsIO.microstructureFile);
       
        
        

        
//        writeConfigFiles(0);
    }

void MicrostructureGenerator::readMicrostructureFile()
{
    this->clear();
    configIO.clear();
    auxIO.clear();
//    configIO.reset(new DDconfigIO<3>(ddBase.simulationParameters.traitsIO.evlFolder));
//    auxIO.reset(new DDauxIO<3>(ddBase.simulationParameters.traitsIO.auxFolder));


    const auto microstructureFiles(TextFileParser(ddBase.simulationParameters.traitsIO.microstructureFile).readStringVector("microstructureFile"));
    for(const auto& pair : microstructureFiles)
    {
        const std::string microstructureFileName(std::filesystem::path(ddBase.simulationParameters.traitsIO.microstructureFile).parent_path().string()+"/"+TextFileParser::removeSpaces(pair.first));
        const std::string microstructureType(TextFileParser(microstructureFileName).readString("type",false));
        const std::string tag(TextFileParser(microstructureFileName).readString("tag",false));
        bool success(false);
        if(microstructureType=="PeriodicDipole")
        {
            success=this->emplace(tag,new PeriodicDipoleGenerator(microstructureFileName)).second;
        }
        else if(microstructureType=="PeriodicLoop")
        {
            success=this->emplace(tag,new PeriodicLoopGenerator(microstructureFileName)).second;
        }
        else if(microstructureType=="PrismaticLoop")
        {
            success=this->emplace(tag,new PrismaticLoopGenerator(microstructureFileName)).second;
        }
        else if(microstructureType=="SphericalInclusions")
        {
            success=this->emplace(tag,new SphericalInclusionsGenerator(microstructureFileName)).second;
        }
        else if(microstructureType=="PolyhedronInclusions")
        {
            success=this->emplace(tag,new PolyhedronInclusionsGenerator(microstructureFileName)).second;
        }
        else if(microstructureType=="Irradiation")
        {
            success=this->emplace(tag,new IrradiationDefectsGenerator(microstructureFileName)).second;
        }
        else if(microstructureType=="VTK")
        {
            success=this->emplace(tag,new VTKGenerator(microstructureFileName)).second;
        }
        else
        {
            std::cout<<"unkown microstructure type "<<microstructureType<<std::endl;
        }
        if(!success)
        {
            throw std::runtime_error("Duplicate microstructure tag "+tag+".");
        }
    }
    
    for(auto& gen : *this)
    {
        if(gen.second->style=="individual")
        {
            gen.second->generateIndividual(*this);
        }
        else if(gen.second->style=="density")
        {
            gen.second->generateDensity(*this);
        }
        else
        {
            throw std::runtime_error("Uknown style for generator "+gen.second->tag);
        }
    }
    
    configIO.finalize();

}

const DDtraitsIO& MicrostructureGenerator::traits() const
{
    return ddBase.simulationParameters.traitsIO;
}


    const DDconfigIO<3>& MicrostructureGenerator::config() const
    {
        return configIO;
    }

    const DDauxIO<3>& MicrostructureGenerator::aux() const
    {
        return auxIO;
    }

DDconfigIO<3>& MicrostructureGenerator::config()
{
    return configIO;
}

DDauxIO<3>& MicrostructureGenerator::aux()
{
    return auxIO;
}

    void MicrostructureGenerator::insertJunctionLoop(const std::vector<VectorDimD>& loopNodePos,
                                                 const std::shared_ptr<PeriodicGlidePlane<3>>& periodicPlane,
                                                 const VectorDimD& b,
                                                 const VectorDimD& unitNormal,
                                                 const VectorDimD& P0,
                                                 const size_t& grainID,
                                                 const DislocationLoopIO<dim>::DislocationLoopType& loopType)
{
    std::vector<PolyPoint> dummyPolyPoints;
    std::vector<std::pair<VectorDimD, const PolyPoint *const>> loopNodePosTemp;
    for(const auto& pos : loopNodePos)
    {
        dummyPolyPoints.push_back(PolyPoint());
        loopNodePosTemp.emplace_back(pos, &dummyPolyPoints.back());
    }
    
    const auto ppi(periodicPlane->polygonPatchIntersection(loopNodePosTemp,true));
    const size_t loopID(insertLoop(b,unitNormal,P0,grainID,loopType));
    std::vector<size_t> loopNodeIDs;
    for(const auto &tup : ppi)
    {
        const VectorDimD loopNodePos(periodicPlane->referencePlane->globalPosition(std::get<0>(tup)));
        const VectorDimD networkNodePos(loopNodePos+std::get<1>(tup));
        const auto networkNodeIter(uniqueNetworkNodeMap.find(networkNodePos));
        if(networkNodeIter==uniqueNetworkNodeMap.end())
        {// no NetworkNode found at current position
            uniqueNetworkNodeMap.emplace(networkNodePos,insertNetworkNode(networkNodePos)); // insert NetworkNode and store its ID
        }
        loopNodeIDs.push_back(insertLoopNode(loopID,loopNodePos,uniqueNetworkNodeMap.at(networkNodePos),std::get<1>(tup),std::get<2>(tup))); // insert LoopNode and store its ID
    }
    insertLoopLinks(loopID,loopNodeIDs);
}

    size_t MicrostructureGenerator::insertLoop(const VectorDimD& b,const VectorDimD& unitNormal,const VectorDimD& P0,const size_t& grainID,const DislocationLoopType& loopType)
    {
        const size_t loopID(configIO.loops().size());
        configIO.loops().emplace_back(loopID, b,unitNormal,P0,grainID,loopType);
        return loopID;
    }

    size_t MicrostructureGenerator::insertLoopNode(const size_t& loopID,const VectorDimD& loopNodePos,const size_t& networkNodeID,const VectorDimD& loopNodeShift,const std::pair<short int,short int>& periodicEdgeIDs)
    {
        const size_t loopNodeID(configIO.loopNodes().size());
        configIO.loopNodes().emplace_back(loopNodeID,loopID,loopNodePos,networkNodeID,loopNodeShift,periodicEdgeIDs);
        return loopNodeID;
    }

    std::vector<size_t> MicrostructureGenerator::insertLoopLinks(const size_t& loopID,const std::vector<size_t>& loopNodeIDs)
    {
        std::vector<size_t> temp;
        for(size_t k=0;k<loopNodeIDs.size();++k)
        {
            const size_t k1=(k+1)<loopNodeIDs.size()? k+1 : 0;
            temp.push_back(configIO.loopLinks().size());
            const size_t sourceNodeID(loopNodeIDs[k ]);
            const size_t   sinkNodeID(loopNodeIDs[k1]);
            const auto sourceNode(configIO.loopNodes()[sourceNodeID]);
            const auto   sinkNode(configIO.loopNodes()[sinkNodeID]);
            configIO.loopLinks().emplace_back(loopID,sourceNodeID,sinkNodeID,(sourceNode.P-sinkNode.P).norm()>FLT_EPSILON,0);
        }
        return temp;
    }


    size_t MicrostructureGenerator::insertNetworkNode(const VectorDimD& networkNodePos)
    {
        const size_t networkNodeID(configIO.nodes().size());
        configIO.nodes().emplace_back(networkNodeID,networkNodePos,Eigen::Matrix<double,1,3>::Zero(),1.0,0);
        return networkNodeID;
    }

    size_t MicrostructureGenerator::insertInclusion(const VectorDimD& pos,const double& R, const Eigen::Matrix<double,dim,dim>& eT, const double& vrc,const int&type)
    {
        const size_t inclusionID(configIO.sphericalInclusions().size()+configIO.polyhedronInclusions().size());
        configIO.sphericalInclusions().emplace_back(inclusionID,pos,R,eT,vrc,type);
        return inclusionID;
    }

size_t MicrostructureGenerator::insertInclusion(const std::map<size_t,Eigen::Vector3d>& polyNodes,const std::map<size_t,std::vector<size_t>>& faceMap, const Eigen::Matrix<double,dim,dim>& eT, const double& vrc,const int&type)
{
    const size_t inclusionID(configIO.sphericalInclusions().size()+configIO.polyhedronInclusions().size());
    configIO.polyhedronInclusions().emplace_back(inclusionID,eT,vrc,type);

    const size_t startNodeID(configIO.polyhedronInclusionNodes().size());
    
    size_t nodeCounter(0);
    for(const auto& node : polyNodes)
    {
        configIO.polyhedronInclusionNodes().emplace_back(startNodeID+nodeCounter,node.second);
        nodeCounter++;
    }
    
    for(const auto& pair : faceMap)
    {
        const size_t& faceID(pair.first);
        for(size_t k=0;k<pair.second.size();++k)
        {
            const size_t k1(k<pair.second.size()-1? k+1 : 0);
            

            const auto sourceIter(polyNodes.find(pair.second[k]));
            const auto sinkIter(polyNodes.find(pair.second[k1]));
            if(sourceIter!=polyNodes.end() && sinkIter!=polyNodes.end())
            {                
                const size_t sourceID(startNodeID +std::distance(polyNodes.begin(),sourceIter));
                const size_t sinkID(startNodeID +std::distance(polyNodes.begin(),sinkIter));
                configIO.polyhedronInclusionEdges().emplace_back(inclusionID,faceID,sourceID,sinkID);
            }
        }
    }
    
    return inclusionID;
}



    /**********************************************************************/
    void MicrostructureGenerator::writeConfigFiles(const size_t& fileID)
    {
        
//        const int outputGlidePlanes(TextFileParser(ddBase.simulationParameters.traitsIO.ddFile).readScalar<int>("outputGlidePlanes",true));

        
//        if(outputGlidePlanes)
//        {
//            for(const auto& loop : configIO.loops())
//            {
//                GlidePlaneKey<dim> loopPlaneKey(loop.P, poly.grain(loop.grainID).singleCrystal->reciprocalLatticeDirection(loop.N));
//                auxIO.glidePlanes().emplace_back(loopPlaneKey);
//            }
//        }
                

            if(outputBinary)
            {
                std::cout<<greenBoldColor<<"Writing configuration to "<<configIO.getBinFilename(fileID)<<defaultColor<<std::endl;
                configIO.writeBin(fileID);
                auxIO.writeBin(fileID);
            }
            else
            {
                std::cout<<greenBoldColor<<"Writing configuration to "<<configIO.getTxtFilename(fileID)<<defaultColor<<std::endl;
                configIO.writeTxt(fileID);
                auxIO.writeTxt(fileID);
            }
        
    }

bool MicrostructureGenerator::allPointsInGrain(const std::vector<VectorDimD>& points,const int& grainID)
{
    bool temp=true;
    for(const auto& point : points)
    {
        temp*=ddBase.mesh.searchRegion(grainID,point).first;
        if(!temp)
        {
            break;
        }
    }
    return temp;
}

}
#endif
