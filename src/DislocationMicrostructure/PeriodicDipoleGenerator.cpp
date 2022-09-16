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
#include <random>
#include <iomanip>

//#include <Simplex.h>
#include <SimplicialMesh.h>
#include <Polycrystal.h>
#include <PolycrystallineMaterialBase.h>
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
#include <PeriodicDipoleGenerator.h>
#include <PlaneLineIntersection.h>

namespace model
{

    PeriodicDipoleGenerator::PeriodicDipoleGenerator(const std::string& fileName) :
    /* init */ MicrostructureGeneratorBase(fileName)
    {
        
    }

    void PeriodicDipoleGenerator::generateDensity(MicrostructureGenerator& mg)
    {
        std::cout<<magentaBoldColor<<"Generating periodic dipole density"<<defaultColor<<std::endl;
        const double targetPeriodicDipoleDensity(this->parser.readScalar<double>("targetPeriodicDipoleDensity",true));
        if(targetPeriodicDipoleDensity>0.0)
        {
            std::mt19937 generator;
            double density=0.0;
            while(density<targetPeriodicDipoleDensity)
            {
                const std::pair<LatticeVector<dim>, int> rp(mg.poly.randomLatticePointInMesh());
                const LatticeVector<dim> L0=rp.first;
                const size_t grainID=rp.second;
                std::uniform_int_distribution<> ssDist(0,mg.poly.grain(grainID).singleCrystal->slipSystems().size()-1);
                const int rSS(ssDist(generator)); // a random SlipSystem
//                const auto& slipSystem(*poly.grain(grainID).singleCrystal->slipSystems()[rSS]);
                std::uniform_int_distribution<> fDist(0,mg.poly.grain(grainID).region.faces().size()-1);
                const int rF(fDist(generator)); // a random face
                auto faceIter(mg.poly.grain(grainID).region.faces().begin());
                std::advance(faceIter,rF);

                try
                {
                    generateSingle(mg,rSS,L0.cartesian(),faceIter->first,100);
                    density+=2.0*faceIter->second->periodicFacePair.first.norm()/mg.mesh.volume()/std::pow(mg.poly.b_SI,2);
                    std::cout<<"periodic dipole density="<<density<<std::endl;
                }
                catch(const std::exception& e)
                {
                    
                }
            }
        }
    }

    void PeriodicDipoleGenerator::generateIndividual(MicrostructureGenerator& mg)
    {
        const std::vector<int> periodicDipoleSlipSystemIDs(this->parser.readArray<int>("periodicDipoleSlipSystemIDs",true));
        
        if(periodicDipoleSlipSystemIDs.size())
        {
            std::cout<<magentaBoldColor<<"Generating individual periodic dipole"<<defaultColor<<std::endl;
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
                generateSingle(mg,periodicDipoleSlipSystemIDs[k],periodicDipolePoints.row(k),periodicDipoleExitFaceIDs[k],periodicDipoleHeights[k]);
            }
        }
        
    }

//    void PeriodicDipoleGenerator::insertJunctionLoop(MicrostructureGenerator& mg,
//                                                     std::map<VectorDimD,size_t,CompareVectorsByComponent<double,dim,float>>& uniqueNetworkNodeMap,
//                                                     const std::vector<VectorDimD>& loopNodePos,
//                                                     const std::shared_ptr<PeriodicGlidePlane<3>>& periodicPlane,
//                                                     const VectorDimD& b,
//                                                     const VectorDimD& unitNormal,
//                                                     const VectorDimD& P0,
//                                                     const size_t& grainID,
//                                                     const DislocationLoopIO<dim>::DislocationLoopType& loopType)
//    {
//        std::vector<PolyPoint> dummyPolyPoints;
//        std::vector<std::pair<VectorDimD, const PolyPoint *const>> loopNodePosTemp;
//        for(const auto& pos : loopNodePos)
//        {
//            dummyPolyPoints.push_back(PolyPoint());
//            loopNodePosTemp.emplace_back(pos, &dummyPolyPoints.back());
//        }
//
//        const auto ppi(periodicPlane->polygonPatchIntersection(loopNodePosTemp));
//        const size_t loopID(mg.insertLoop(b,unitNormal,P0,grainID,loopType));
//        std::vector<size_t> loopNodeIDs;
//        for(const auto &tup : ppi)
//        {
//            const VectorDimD loopNodePos(periodicPlane->referencePlane->globalPosition(std::get<0>(tup)));
//            const VectorDimD networkNodePos(loopNodePos+std::get<1>(tup));
//            const auto networkNodeIter(uniqueNetworkNodeMap.find(networkNodePos));
//            if(networkNodeIter==uniqueNetworkNodeMap.end())
//            {// no NetworkNode found at current position
//                uniqueNetworkNodeMap.emplace(networkNodePos,mg.insertNetworkNode(networkNodePos)); // insert NetworkNode and store its ID
//            }
//            loopNodeIDs.push_back(mg.insertLoopNode(loopID,loopNodePos,uniqueNetworkNodeMap.at(networkNodePos),std::get<1>(tup),std::get<2>(tup))); // insert LoopNode and store its ID
//        }
//        mg.insertLoopLinks(loopID,loopNodeIDs);
//    }


    void PeriodicDipoleGenerator::generateSingle(MicrostructureGenerator& mg,const int& rSS,const VectorDimD& dipolePoint,const int& exitFaceID,const int dipoleHeight)
    {
        std::pair<bool,const Simplex<dim,dim>*> found(mg.mesh.search(dipolePoint));
        if(!found.first)
        {
            std::cout<<"Point "<<dipolePoint.transpose()<<" is outside mesh. EXITING."<<std::endl;
            exit(EXIT_FAILURE);
        }
        
        const int grainID(found.second->region->regionID);
        assert(mg.poly.grains.size()==1 && "Periodic dislocations only supported for single crystals");
        const auto& grain(mg.poly.grain(grainID));
        
        if(rSS>=0 && rSS<int(grain.singleCrystal->slipSystems().size()))
        {
            const auto periodicFaceIter(grain.region.faces().find(exitFaceID));
            if(periodicFaceIter!=grain.region.faces().end())
            {
                const auto periodicFaceA(periodicFaceIter->second);
                const auto periodicFaceB(periodicFaceA->periodicFacePair.second);

                if(periodicFaceB!=nullptr)
                {

                    const auto& faceAshift(periodicFaceA->periodicFacePair.first);
                    const auto faceAlatticeShift(grain.singleCrystal->latticeVector(faceAshift));
                    
                    const auto& slipSystem(*grain.singleCrystal->slipSystems()[rSS]);
                    
                    if(slipSystem.n.dot(faceAlatticeShift)==0)
                    {

                        //                        const std::pair<bool,long int> heightPair=LatticePlane::computeHeight(slipSystem.n,dipolePoint);
                        
                        const long int planeIndex(slipSystem.n.closestPlaneIndexOfPoint(dipolePoint));
                        GlidePlaneKey<3> glidePlaneKey(planeIndex, slipSystem.n);
                        std::shared_ptr<PeriodicGlidePlane<3>> glidePlane(mg.periodicGlidePlaneFactory.get(glidePlaneKey));
                        //                        const VectorDimD P0(glidePlane->snapToPlane(dipolePoint));
                        const VectorDimD P0(grain.singleCrystal->snapToLattice(dipolePoint).cartesian());
                        
                        PlaneLineIntersection<3> pliA(periodicFaceA->center(),periodicFaceA->outNormal(),P0,faceAshift);
                        PlaneLineIntersection<3> pliB(periodicFaceB->center(),periodicFaceB->outNormal(),P0,faceAshift);
                        
                        if(pliA.type==PlaneLineIntersection<3>::INCIDENT && pliB.type==PlaneLineIntersection<3>::INCIDENT)
                        {
                            const VectorDimD AB(pliB.P-pliA.P);
                            if((AB-faceAshift).norm()<FLT_EPSILON)
                            {

                                GlidePlaneKey<3> parallelGlidePlaneKey(planeIndex+dipoleHeight, slipSystem.n);
                                std::shared_ptr<PeriodicGlidePlane<3>> parallelglidePlane(mg.periodicGlidePlaneFactory.get(parallelGlidePlaneKey));

                                GlidePlaneKey<3> prismaticPlaneKey(P0, grain.singleCrystal->reciprocalLatticeDirection(glidePlane->referencePlane->unitNormal.cross(AB)));
                                std::shared_ptr<PeriodicGlidePlane<3>> prismaticGlidePlane(mg.periodicGlidePlaneFactory.get(prismaticPlaneKey));

                                if(parallelglidePlane && prismaticGlidePlane)
                                {

                                    
                                    std::map<VectorDimD,size_t,CompareVectorsByComponent<double,dim,float>> uniqueNetworkNodeMap; // networkNodePosition->networkNodeID
                                    // The prismatic loop
                                    const int nShift(2);
                                    const VectorDimD lShift(nShift*AB);
                                    const VectorDimD rShift(lShift+AB);
                                    std::vector<VectorDimD> prismaticNodePos;
                                    prismaticNodePos.push_back(0.5*(pliB.P+pliA.P)+lShift); //HERE ADD N*AB
                                    prismaticNodePos.push_back(0.5*(pliB.P+pliA.P)+rShift);
                                    prismaticNodePos.push_back(parallelglidePlane->referencePlane->snapToPlane(0.5*(pliB.P+pliA.P)+rShift));
                                    prismaticNodePos.push_back(parallelglidePlane->referencePlane->snapToPlane(0.5*(pliB.P+pliA.P)+lShift));
                                    

                                    
//                                    insertJunctionLoop(mg,uniqueNetworkNodeMap,prismaticNodePos,prismaticGlidePlane,
//                                                       slipSystem.s.cartesian(),prismaticGlidePlane->referencePlane->unitNormal,
//                                                       P0,grainID,DislocationLoopIO<dim>::SESSILELOOP);
                                    
                                    mg.insertJunctionLoop(prismaticNodePos,prismaticGlidePlane,
                                                       slipSystem.s.cartesian(),prismaticGlidePlane->referencePlane->unitNormal,
                                                       P0,grainID,DislocationLoopIO<dim>::SESSILELOOP);
                                    

                                    
                                    // First glide loop
                                    const double glideStep=1.0;
                                    std::vector<VectorDimD> firstNodePos;
                                    firstNodePos.push_back(0.5*(pliB.P+pliA.P)+lShift);
                                    firstNodePos.push_back(0.5*(pliB.P+pliA.P)+rShift);
                                    firstNodePos.push_back(0.5*(pliB.P+pliA.P)+rShift+glideStep*prismaticGlidePlane->referencePlane->unitNormal);
                                    firstNodePos.push_back(0.5*(pliB.P+pliA.P)+lShift+glideStep*prismaticGlidePlane->referencePlane->unitNormal);
                                    
//                                    insertJunctionLoop(mg,uniqueNetworkNodeMap,firstNodePos,glidePlane,
//                                                       -slipSystem.s.cartesian(),glidePlane->referencePlane->unitNormal,
//                                                       P0,grainID,DislocationLoopIO<dim>::GLISSILELOOP);
                                    
                                    mg.insertJunctionLoop(firstNodePos,glidePlane,
                                                       -slipSystem.s.cartesian(),glidePlane->referencePlane->unitNormal,
                                                       P0,grainID,DislocationLoopIO<dim>::GLISSILELOOP);


                                    
                                    // Second glide loop
                                    std::vector<VectorDimD> secondNodePos;
                                    for(const auto& pos : firstNodePos)
                                    {
                                        secondNodePos.push_back(parallelglidePlane->referencePlane->snapToPlane(pos));
                                    }
                                    

//                                    insertJunctionLoop(mg,uniqueNetworkNodeMap,secondNodePos,parallelglidePlane,
//                                                       slipSystem.s.cartesian(),parallelglidePlane->referencePlane->unitNormal,
//                                                       parallelglidePlane->referencePlane->snapToPlane(P0),grainID,DislocationLoopIO<dim>::GLISSILELOOP);
                                    
                                    mg.insertJunctionLoop(secondNodePos,parallelglidePlane,
                                                       slipSystem.s.cartesian(),parallelglidePlane->referencePlane->unitNormal,
                                                       parallelglidePlane->referencePlane->snapToPlane(P0),grainID,DislocationLoopIO<dim>::GLISSILELOOP);
                                }
                                else
                                {
                                    std::cout<<"Cannot create glide planes"<<std::endl;
                                }
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
                        throw std::runtime_error("planeNormal of slipSystem "+std::to_string(rSS)+" is not othogonal to faceNormal "+std::to_string(exitFaceID));
                    }
                }
                else
                {
                    throw std::runtime_error("Mesh face "+std::to_string(exitFaceID)+" is not a periodic face.");
                }
            }
            else
            {
                throw std::runtime_error("Mesh face "+std::to_string(exitFaceID)+" not found.");
            }
        }
        else
        {
            if(rSS<0)
            {
                std::cout<<"Skipping slip system "<<rSS<<std::endl;
            }
            else
            {
                throw std::runtime_error("slipSystem "+std::to_string(rSS)+" not found, skipping.");
            }
        }
    }

}
#endif
