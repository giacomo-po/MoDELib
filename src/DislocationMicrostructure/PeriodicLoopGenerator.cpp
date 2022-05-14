/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PeriodicDipoleGenerator_cpp_
#define model_PeriodicDipoleGenerator_cpp_


#include <chrono>
#include <random>
#include <cmath>
#include <list>
#include <assert.h>
#include <Eigen/LU>
#include <Eigen/Cholesky>
#include <limits>

//#include <Simplex.h>
#include <SimplicialMesh.h>
#include <Polycrystal.h>
#include <PolycrystallineMaterial.h>
#include <LatticeModule.h>
//#include <PlaneMeshIntersection.h>
#include <DislocationNodeIO.h>
#include <DislocationLoopIO.h>
#include <DislocationLoopLinkIO.h>
#include <DislocationLoopNodeIO.h>
#include <DDconfigIO.h>
#include <DDauxIO.h>

#include <DislocationLinkingNumber.h>
#include <TextFileParser.h>
#include <DislocationInjector.h>
#include <MeshBoundarySegment.h>
//#include <ConfinedDislocationObject.h>
#include <GlidePlaneModule.h>
#include <MeshModule.h>
#include <Plane.h>
#include <MicrostructureGenerator.h>
#include <PeriodicLoopGenerator.h>
#include <PlaneLineIntersection.h>

namespace model
{

    PeriodicLoopGenerator::PeriodicLoopGenerator(const std::string& fileName) :
    /* init */ MicrostructureGeneratorBase(fileName)
    {
    
    }

    void PeriodicLoopGenerator::generate(MicrostructureGenerator& mg)
    {
        generateSingle(mg);
    }

    void PeriodicLoopGenerator::generateDensity(MicrostructureGenerator& mg)
    {
        
    }

    void PeriodicLoopGenerator::generateSingle(MicrostructureGenerator& mg)
    {
        const std::vector<int> periodicDipoleSlipSystemIDs(this->parser.readArray<int>("periodicLoopSlipSystemIDs",true));
        
        if(periodicDipoleSlipSystemIDs.size())
        {
            std::cout<<magentaBoldColor<<"Generating individual periodic straight dislocations"<<defaultColor<<std::endl;
            const std::vector<int> periodicDipoleExitFaceIDs(this->parser.readArray<int>("periodicDipoleExitFaceIDs",true));
            const Eigen::Matrix<double,Eigen::Dynamic,dim> periodicDipolePoints(this->parser.readMatrix<double>("periodicDipolePoints",periodicDipoleSlipSystemIDs.size(),dim,true));
            const std::vector<double> periodicDipoleHeights(this->parser.readArray<double>("periodicDipoleHeights",true));
            
            if(periodicDipoleSlipSystemIDs.size()!=periodicDipoleExitFaceIDs.size())
            {
                throw std::runtime_error("periodicDipoleSlipSystemIDs.size()="+std::to_string(periodicDipoleSlipSystemIDs.size())+" NOT EQUAL TO periodicDipoleExitFaceIDs.size()="+std::to_string(periodicDipoleExitFaceIDs.size()));
            }
            if(int(periodicDipoleSlipSystemIDs.size())!=periodicDipolePoints.rows())
            {
                throw std::runtime_error("periodicDipoleSlipSystemIDs.size()="+std::to_string(periodicDipoleSlipSystemIDs.size())+" NOT EQUAL TO periodicDipolePoints.size()="+std::to_string(periodicDipolePoints.size()));
            }
            if(periodicDipoleSlipSystemIDs.size()!=periodicDipoleHeights.size())
            {
                throw std::runtime_error("periodicDipoleSlipSystemIDs.size()="+std::to_string(periodicDipoleSlipSystemIDs.size())+" NOT EQUAL TO periodicDipoleHeights.size()="+std::to_string(periodicDipoleHeights.size()));
            }

            for(size_t k=0;k<periodicDipoleSlipSystemIDs.size();++k)
            {
                const int& rSS(periodicDipoleSlipSystemIDs[k]);


                if(rSS>=0)
                {
                    std::pair<bool,const Simplex<dim,dim>*> found(mg.mesh.search(periodicDipolePoints.row(k)));
                    if(!found.first)
                    {
                        std::cout<<"Point "<<periodicDipolePoints.row(k)<<" is outside mesh. EXITING."<<std::endl;
                        exit(EXIT_FAILURE);
                    }

                    const int grainID(found.second->region->regionID);
                    assert(mg.poly.grains().size()==1 && "Periodic dislocations only supported for single crystals");
                    const auto& grain(mg.poly.grain(grainID));

                    const auto periodicFaceIter(grain.region.faces().find(periodicDipoleExitFaceIDs[k]));
                    if(periodicFaceIter!=grain.region.faces().end())
                    {
                        const auto periodicFaceA(periodicFaceIter->second);
                        const auto periodicFaceB(periodicFaceA->periodicFacePair.second);

                        if(periodicFaceB!=nullptr)
                        {
                            const auto& faceAshift(periodicFaceA->periodicFacePair.first);
                            const auto faceAlatticeShift(mg.poly.grain(grainID).latticeVector(faceAshift));
                            
                            const auto& slipSystem(*mg.poly.grain(grainID).slipSystems()[rSS]);
                            
                            if(slipSystem.n.dot(faceAlatticeShift)==0)
                            {
                                const std::pair<bool,long int> heightPair=LatticePlane::computeHeight(slipSystem.n,periodicDipolePoints.row(k));

                                GlidePlaneKey<3> glidePlaneKey(heightPair.second, slipSystem.n);
                                std::shared_ptr<GlidePlane<3>> glidePlane(mg.glidePlaneFactory.getFromKey(glidePlaneKey));
                                const VectorDimD P0(glidePlane->snapToPlane(periodicDipolePoints.row(k)));
                                PlaneLineIntersection<3> pliA(periodicFaceA->center(),periodicFaceA->outNormal(),P0,faceAshift);
                                PlaneLineIntersection<3> pliB(periodicFaceB->center(),periodicFaceB->outNormal(),P0,faceAshift);

                                if(pliA.type==PlaneLineIntersection<3>::INCIDENT && pliB.type==PlaneLineIntersection<3>::INCIDENT)
                                {
                                    const VectorDimD AB(pliB.P-pliA.P);
                                    if((AB-faceAshift).norm()<FLT_EPSILON)
                                    {
                                        
                                        GlidePlaneKey<3> parallelGlidePlaneKey(heightPair.second+periodicDipoleHeights[k], slipSystem.n);
                                        std::shared_ptr<GlidePlane<3>> parallelglidePlane(mg.glidePlaneFactory.getFromKey(parallelGlidePlaneKey));
                                        
                                        std::vector<PolyPoint> dummyPolyPoints;
                                        std::vector<std::pair<VectorDimD, const PolyPoint *const>> loopNodePosTemp;
                                        
                                        dummyPolyPoints.push_back(PolyPoint());
                                        loopNodePosTemp.push_back(std::make_pair(0.5*(pliB.P+pliA.P), &dummyPolyPoints.back()));
                                        dummyPolyPoints.push_back(PolyPoint());
                                        loopNodePosTemp.push_back(std::make_pair(0.5*(pliB.P+pliA.P)+AB, &dummyPolyPoints.back()));
                                        dummyPolyPoints.push_back(PolyPoint());
                                        loopNodePosTemp.push_back(std::make_pair(parallelglidePlane->snapToPlane(0.5*(pliB.P+pliA.P)+AB), &dummyPolyPoints.back()));
                                        dummyPolyPoints.push_back(PolyPoint());
                                        loopNodePosTemp.push_back(std::make_pair(parallelglidePlane->snapToPlane(0.5*(pliB.P+pliA.P)), &dummyPolyPoints.back()));

                                        GlidePlaneKey<3> prismaticPlaneKey(P0, grain.reciprocalLatticeDirection(glidePlane->unitNormal.cross(AB)));
                                        std::shared_ptr<PeriodicGlidePlane<3>> prismaticGlidePlane(mg.periodicGlidePlaneFactory.get(prismaticPlaneKey));

                                        
                                        const auto ppi(prismaticGlidePlane->polygonPatchIntersection(loopNodePosTemp));

                                        const size_t prismaticLoopID(mg.insertLoop(slipSystem.s.cartesian(),prismaticGlidePlane->referencePlane->unitNormal,P0,grainID,DislocationLoopIO<dim>::SESSILELOOP));

                                        std::vector<size_t> loopNodeIDs;
                                        std::map<VectorDimD,size_t,CompareVectorsByComponent<double,dim,float>> uniqueNetworkNodeMap; // networkNodePosition->networkNodeID
                                        for(const auto &tup : ppi)
                                        {
                                            const VectorDimD loopNodePos(prismaticGlidePlane->referencePlane->globalPosition(std::get<0>(tup)));
                                            const VectorDimD networkNodePos(loopNodePos+std::get<1>(tup));
                                            const auto networkNodeIter(uniqueNetworkNodeMap.find(networkNodePos));
                                            if(networkNodeIter==uniqueNetworkNodeMap.end())
                                            {// no NetworkNode found at current position
                                                uniqueNetworkNodeMap.emplace(networkNodePos,mg.insertNetworkNode(networkNodePos)); // insert NetworkNode and store its ID
                                            }
                                            loopNodeIDs.push_back(mg.insertLoopNode(prismaticLoopID,loopNodePos,uniqueNetworkNodeMap.at(networkNodePos),std::get<1>(tup),std::get<2>(tup))); // insert LoopNode and store its ID
                                        }
                                        mg.insertLoopLinks(prismaticLoopID,loopNodeIDs);
                                    }
                                    else
                                    {
                                        throw std::runtime_error("periodic line intersection with faces is not the face shift vector");
                                    }
                                }
                                else
                                {
                                    throw std::runtime_error("periodic line direction does not form an incident intersecitn with periodic faces");
                                }

                            }
                            else
                            {
                                std::cout<<"Normal to slipSystem "<<rSS<<" is not othogonal to face "<<periodicDipoleExitFaceIDs[k]<<std::endl;
                            }
                        }
                        else
                        {
                            throw std::runtime_error("Mesh face "+std::to_string(periodicDipoleExitFaceIDs[k])+" is not a periodic face.");
                        }
                    }
                    else
                    {
                        throw std::runtime_error("Mesh face "+std::to_string(periodicDipoleExitFaceIDs[k])+" not found.");
                    }
                }
                else
                {
                    std::cout<<"negative slip system ID. Skipping entries."<<std::endl;
                }
            }
        }
        
    }

}
#endif
