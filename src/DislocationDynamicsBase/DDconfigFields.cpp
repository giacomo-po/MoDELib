/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2020 by Danny Perez <danny_perez@lanl.gov>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DDconfigFields_cpp_
#define model_DDconfigFields_cpp_

#include <DDconfigFields.h>
#include <StressStraight.h>


namespace model
{

    template <int dim>
    DDconfigFields<dim>::DDconfigFields(DislocationDynamicsBase<dim>& ddBase_in,const DDconfigIO<dim>& configIO_in):
    /* init */ ddBase(ddBase_in)
    /* init */,periodicShifts(ddBase.mesh.periodicShifts(ddBase.simulationParameters.periodicImageSize))
    /* init */,configIO(configIO_in)
    {
        
    }


    template <int dim>
    void DDconfigFields<dim>::updateConfiguration()
    {
        
        // clean up
        loopPatches().clear();
        segments()=configIO.segments();
        polyhedronInclusionNodes().clear();
        eshelbyInclusions().clear();
        EshelbyInclusionBase<dim>::force_count(0);
        
        
        // Update loop patches
        for(const auto& pair : configIO.loopNodeSequence())
        {
            const auto& loopID(pair.first);
            const auto& loopIO(configIO.loops()[configIO.loopMap().at(loopID)]);
            const auto& grain(ddBase.poly.grain(loopIO.grainID));
            GlidePlaneKey<3> glidePlaneKey(loopIO.P, grain.singleCrystal->reciprocalLatticeDirection(loopIO.N));
            std::shared_ptr<PeriodicGlidePlane<3>> periodicGlidePlane(ddBase.periodicGlidePlaneFactory.get(glidePlaneKey));
            
            std::vector<Eigen::Matrix<double,3,1>> nodeShifts;
            std::vector<Eigen::Matrix<double,3,1>> nodePos;
            
            for(const auto& loopNodeID : pair.second)
            {
                const auto& loopNodeIO(configIO.loopNodes()[configIO.loopNodeMap().at(loopNodeID)]);
                nodeShifts.push_back(loopNodeIO.periodicShift);
                if(loopNodeIO.edgeIDs.first<0 && loopNodeIO.edgeIDs.second<0)
                {
                    nodePos.push_back(loopNodeIO.P);
                }
            }
            
            DislocationLoopPatches<3> currentPatches(periodicGlidePlane);
            currentPatches.update(nodeShifts,nodePos);
            loopPatches().emplace(loopID,currentPatches);
        }
        
        // update inclusions
        for(const auto& inclusion : configIO.sphericalInclusions())
        {
            //        std::cout<<"Creating spherical inclusion "<<inclusion.inclusionID<<std::endl;
            const std::pair<bool,const Simplex<dim,dim>*> searchPair(ddBase.mesh.search(inclusion.C));
            if(searchPair.first)
            {
                
                const auto& grain(ddBase.poly.grain(searchPair.second->region->regionID));
                if(inclusion.phaseID<int(grain.singleCrystal->secondPhases().size()))
                {
                    const auto secondPhase(grain.singleCrystal->secondPhases()[inclusion.phaseID]);
                    EshelbyInclusionBase<dim>::set_count(inclusion.inclusionID);
                    
                    
                    std::shared_ptr<EshelbyInclusionBase<dim>> iptr(new SphericalInclusion<dim>(inclusion.C,inclusion.a,inclusion.eT,ddBase.poly.nu,ddBase.poly.mu,inclusion.mobilityReduction,inclusion.phaseID,secondPhase));
                    
                    eshelbyInclusions().emplace(inclusion.inclusionID,iptr);
                }
                else
                {
                    throw std::runtime_error("phaseID does not exist in grain.");
                }
            }
        }
        
        for(const auto& piNode : configIO.polyhedronInclusionNodes())
        {
            polyhedronInclusionNodes().emplace(piNode.nodeID,piNode);
        }
        
        std::map<size_t,std::map<size_t,std::vector<size_t>>> faceMap;
        for(const auto& edge : configIO.polyhedronInclusionEdges())
        {
            const size_t& iID(edge.inclusionID);
            const size_t& fID(edge.faceID);
            const size_t& sourceID(edge.sourceID);
            faceMap[iID][fID].push_back(sourceID);
        }
        
        
        for(const auto& inclusion : configIO.polyhedronInclusions())
        {
            //        std::cout<<"Creating polyhedron inclusion "<<inclusion.inclusionID<<std::endl;
            
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
                        const std::pair<bool,const Simplex<dim,dim>*> searchPair(ddBase.mesh.search(nodePtr->P));
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
                        const auto& grain(ddBase.poly.grain(*grainIDs.begin()));
                        if(inclusion.phaseID<int(grain.singleCrystal->secondPhases().size()))
                        {
                            const auto secondPhase(grain.singleCrystal->secondPhases()[inclusion.phaseID]);
                            EshelbyInclusionBase<dim>::set_count(inclusion.inclusionID);
                            std::shared_ptr<EshelbyInclusionBase<dim>> iptr(new PolyhedronInclusion<dim>( polyhedronInclusionNodes(),faces,inclusion.eT,ddBase.poly.nu,ddBase.poly.mu,inclusion.mobilityReduction,inclusion.phaseID,secondPhase));
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

    template <int dim>
    const std::map<size_t,DislocationLoopPatches<dim>>& DDconfigFields<dim>::loopPatches() const
    {
        return *this;
    }

    template <int dim>
    std::map<size_t,DislocationLoopPatches<dim>>& DDconfigFields<dim>::loopPatches()
    {
        return *this;
    }

    template <int dim>
    const std::map<std::pair<size_t,size_t>,DislocationSegmentIO<dim>>& DDconfigFields<dim>::segments() const
    {
        return *this;
    }

    template <int dim>
    std::map<std::pair<size_t,size_t>,DislocationSegmentIO<dim>>& DDconfigFields<dim>::segments()
    {
        return *this;
    }

    template <int dim>
    const typename DDconfigFields<dim>::PolyhedronInclusionNodeContainerType& DDconfigFields<dim>::polyhedronInclusionNodes() const
    {
        return *this;
    }

    template <int dim>
    typename DDconfigFields<dim>::PolyhedronInclusionNodeContainerType& DDconfigFields<dim>::polyhedronInclusionNodes()
    {
        return *this;
    }

    template <int dim>
    const typename DDconfigFields<dim>::EshelbyInclusionContainerType& DDconfigFields<dim>::eshelbyInclusions() const
    {
        return *this;
    }

    template <int dim>
    typename DDconfigFields<dim>::EshelbyInclusionContainerType& DDconfigFields<dim>::eshelbyInclusions()
    {
        return *this;
    }

    template <int dim>
    double DDconfigFields<dim>::solidAngle(const VectorDim& x) const
    {
        double temp(0.0);
        for(const auto& patch : loopPatches())
        {
            for(const auto& shift : periodicShifts)
            {
                temp+=patch.second.solidAngle(x+shift);
            }
        }
        return temp;
    }

template <int dim>
typename  DDconfigFields<dim>::MatrixDim DDconfigFields<dim>::dislocationStress(const VectorDim& x) const
{
    MatrixDim temp(MatrixDim::Zero());
    for (const auto& segment : segments())
    {
        if(segment.second.meshLocation==0 || segment.second.meshLocation==2)
        {// segment inside mesh or on grain boundary
            auto itSource(configIO.nodeMap().find(segment.second.sourceID)); //source
            auto   itSink(configIO.nodeMap().find(segment.second.sinkID)); //sink
            if(itSource!=configIO.nodeMap().end() && itSink!=configIO.nodeMap().end())
            {
                const auto& sourceNode(configIO.nodes()[itSource->second]);
                const auto&   sinkNode(configIO.nodes()[itSink->second]);
                StressStraight<3> ss(ddBase.poly,sourceNode.P,sinkNode.P,segment.second.b);
                for(const auto& shift : periodicShifts)
                {
                    temp+=ss.stress(x+shift);
                }
            }
            else
            {

            }
        }
    }
    return temp;
}

template <int dim>
typename  DDconfigFields<dim>::MatrixDim DDconfigFields<dim>::inclusionStress(const VectorDim& x) const
{
    MatrixDim temp(MatrixDim::Zero());
    for(const auto& inclusion : eshelbyInclusions() )
    {
        for(const auto& shift : periodicShifts)
        {
            temp+=inclusion.second->stress(x+shift);
        }
    }
    return temp;
}

    template class DDconfigFields<3>;

}
#endif
