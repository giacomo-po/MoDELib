/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_IrradiationDefectsGenerator_cpp_
#define model_IrradiationDefectsGenerator_cpp_

#include <numbers>
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
#include <GlidePlaneModule.h>
#include <MeshModule.h>
#include <Plane.h>
#include <MicrostructureGenerator.h>
#include <IrradiationDefectsGenerator.h>
#include <PlaneLineIntersection.h>

namespace model
{

    IrradiationDefectsGenerator::IrradiationDefectsGenerator(const std::string& fileName) :
    /* init */ MicrostructureGeneratorBase(fileName)
    {
        
    }


    void IrradiationDefectsGenerator::generateDensity(MicrostructureGenerator& mg)
    {
        if(mg.ddBase.poly.crystalStructure=="FCC")
        {
            generateDensityFCC(mg);
        }
        else if(mg.ddBase.poly.crystalStructure=="BCC")
        {
            
        }
        else if(mg.ddBase.poly.crystalStructure=="HCP")
        {
            
        }
        else
        {
            throw std::runtime_error("Unknown crystal structure"+mg.ddBase.poly.crystalStructure);
        }
    }

    void IrradiationDefectsGenerator::generateIndividual(MicrostructureGenerator& mg)
    {
        if(mg.ddBase.poly.crystalStructure=="FCC")
        {
            generateIndividualFCC(mg);
        }
        else if(mg.ddBase.poly.crystalStructure=="BCC")
        {
            
        }
        else if(mg.ddBase.poly.crystalStructure=="HCP")
        {
            
        }
        else
        {
            throw std::runtime_error("Unknown crystal structure"+mg.ddBase.poly.crystalStructure);
        }

    }

bool IrradiationDefectsGenerator::generateSingleSFT(MicrostructureGenerator& mg,const int& planeID,const VectorDimD& basePoint,const bool& inverted,const int& sftSize)
{
    if(sftSize>1.0)
    {
        std::pair<bool,const Simplex<dim,dim>*> found=mg.ddBase.mesh.search(basePoint);
        if(found.first)
        {
            const int grainID=found.second->region->regionID;
            const auto& grain(mg.ddBase.poly.grain(grainID));
            const auto& singleCrystal(grain.singleCrystal);


            if(planeID>=0 && planeID<4)
            {
//                const auto basePlane(inverted? std::shared_ptr<LatticePlaneBase>(new LatticePlaneBase(singleCrystal->planeNormals()[planeID]->primitiveVectors.second*(-1),singleCrystal->planeNormals()[planeID]->primitiveVectors.first*(-1))) : singleCrystal->planeNormals()[planeID]);
                const auto basePlane(inverted? std::shared_ptr<LatticePlaneBase>(new LatticePlaneBase(singleCrystal->planeNormals()[planeID]->primitiveVectors.second,singleCrystal->planeNormals()[planeID]->primitiveVectors.first)) : singleCrystal->planeNormals()[planeID]);
                const auto unitNormal(basePlane->cartesian().normalized());
                const auto faceBaryShift((basePlane->primitiveVectors.first+basePlane->primitiveVectors.second).cartesian()*sftSize/3.0);
                const VectorDimD corner(singleCrystal->snapToLattice(basePoint-faceBaryShift).cartesian());

                std::vector<VectorDimD> nodePos;
                nodePos.push_back(corner);
                nodePos.push_back(corner+basePlane->primitiveVectors.first.cartesian()*sftSize);
                nodePos.push_back(corner+basePlane->primitiveVectors.second.cartesian()*sftSize);
                nodePos.push_back((nodePos[0]+nodePos[1]+nodePos[2])/3.0-unitNormal*sftSize*sqrt(2.0/3.0));
//                nodePos.push_back((nodePos[0]+nodePos[1]+nodePos[2])/3.0+unitNormal*sftSize*sqrt(2.0/3.0));

                
                if(mg.allPointsInGrain(nodePos,grainID))
                {
                    //3 ___________2____________ 3
                    //  \          /\          /
                    //   \  230   /  \  132   /
                    //    \      /    \      /
                    //     \    /      \    /
                    //      \  /  012   \  /
                    //       \/__________\/
                    //      0 \          / 1
                    //         \  031   /
                    //          \      /
                    //           \    /
                    //            \  /
                    //             \/
                    //              3

                    if(true)
                    {// SFT loop 0->1->2 (base Frank loop)
                        std::vector<VectorDimD> facePos012;
                        facePos012.push_back(nodePos[0]);
                        facePos012.push_back(nodePos[1]);
                        facePos012.push_back(nodePos[2]);
                        const auto n012(singleCrystal->reciprocalLatticeDirection((nodePos[1]-nodePos[0]).cross(nodePos[2]-nodePos[1])));
                        const VectorDimD b012(n012.interplaneVector());
                        const std::pair<bool,long int> heightPair012=LatticePlane::computeHeight(n012,nodePos[0]);
                        if(!heightPair012.first)
                        {
                            throw std::runtime_error("Cannot determine plane of SFT face 012.");
                        }
                        GlidePlaneKey<3> basePlaneKey012(heightPair012.second, n012);
                        std::shared_ptr<PeriodicGlidePlane<3>> glidePlane012(mg.ddBase.periodicGlidePlaneFactory.get(basePlaneKey012));
                        mg.insertJunctionLoop(facePos012, glidePlane012,
                                              b012,n012.cartesian().normalized(),
                                              nodePos[0],grainID,DislocationLoopIO<dim>::SESSILELOOP);
                    }

                    if(true)
                    {// SFT loop 2->3->0 (side Shockley loop)
                        std::vector<VectorDimD> facePos230;
                        facePos230.push_back(nodePos[2]);
                        facePos230.push_back(nodePos[3]);
                        facePos230.push_back(nodePos[0]);
                        const auto n230(singleCrystal->reciprocalLatticeDirection((nodePos[3]-nodePos[2]).cross(nodePos[0]-nodePos[3])));
                        const VectorDimD b230((0.5*(nodePos[2]+nodePos[0])-nodePos[3]).normalized()*sqrt(3.0)/3.0);
                        const std::pair<bool,long int> heightPair230=LatticePlane::computeHeight(n230,nodePos[2]);
                        if(!heightPair230.first)
                        {
                            throw std::runtime_error("Cannot determine plane of SFT face 230.");
                        }
                        GlidePlaneKey<3> basePlaneKey230(heightPair230.second, n230);
                        std::shared_ptr<PeriodicGlidePlane<3>> glidePlane230(mg.ddBase.periodicGlidePlaneFactory.get(basePlaneKey230));
                        mg.insertJunctionLoop(facePos230, glidePlane230,
                                              b230,n230.cartesian().normalized(),
                                              nodePos[2],grainID,DislocationLoopIO<dim>::GLISSILELOOP);
                    }
                    
                    if(true)
                    {// SFT loop 0->3->1 (side Shockley loop)
                        std::vector<VectorDimD> facePos031;
                        facePos031.push_back(nodePos[0]);
                        facePos031.push_back(nodePos[3]);
                        facePos031.push_back(nodePos[1]);
                        const auto n031(singleCrystal->reciprocalLatticeDirection((nodePos[3]-nodePos[0]).cross(nodePos[1]-nodePos[3])));
                        const VectorDimD b031((0.5*(nodePos[0]+nodePos[1])-nodePos[3]).normalized()*sqrt(3.0)/3.0);
                        const std::pair<bool,long int> heightPair031=LatticePlane::computeHeight(n031,nodePos[0]);
                        if(!heightPair031.first)
                        {
                            throw std::runtime_error("Cannot determine plane of SFT face 031.");
                        }
                        GlidePlaneKey<3> basePlaneKey031(heightPair031.second, n031);
                        std::shared_ptr<PeriodicGlidePlane<3>> glidePlane031(mg.ddBase.periodicGlidePlaneFactory.get(basePlaneKey031));
                        mg.insertJunctionLoop(facePos031, glidePlane031,
                                              b031,n031.cartesian().normalized(),
                                              nodePos[0],grainID,DislocationLoopIO<dim>::GLISSILELOOP);
                    }
                    
                    if(true)
                    {// SFT loop 1->3->2 (side Shockley loop)
                        std::vector<VectorDimD> facePos132;
                        facePos132.push_back(nodePos[1]);
                        facePos132.push_back(nodePos[3]);
                        facePos132.push_back(nodePos[2]);
                        const auto n132(singleCrystal->reciprocalLatticeDirection((nodePos[3]-nodePos[1]).cross(nodePos[2]-nodePos[3])));
                        const VectorDimD b132((0.5*(nodePos[1]+nodePos[2])-nodePos[3]).normalized()*sqrt(3.0)/3.0);
                        const std::pair<bool,long int> heightPair132=LatticePlane::computeHeight(n132,nodePos[1]);
                        if(!heightPair132.first)
                        {
                            throw std::runtime_error("Cannot determine plane of SFT face132.");
                        }
                        GlidePlaneKey<3> basePlaneKey132(heightPair132.second, n132);
                        std::shared_ptr<PeriodicGlidePlane<3>> glidePlane132(mg.ddBase.periodicGlidePlaneFactory.get(basePlaneKey132));
                        mg.insertJunctionLoop(facePos132, glidePlane132,
                                              b132,n132.cartesian().normalized(),
                                              nodePos[1],grainID,DislocationLoopIO<dim>::GLISSILELOOP);
                    }
                    return true;
                }
                else
                {
                    std::cout<<"SFT points outside mesh, skipping generation."<<std::endl;
                    return false;
                }
            }
            else
            {
                std::cout<<"SFT planeID not a valid index into planeNormals, skipping generation."<<std::endl;
                return false;
            }
        }
        else
        {
            std::cout<<"SFT basePoint outside mesh, skipping generation."<<std::endl;
            return false;
        }
    }
    else
    {
        std::cout<<"SFT size <1b, skipping generation."<<std::endl;
        return false;
    }
    
}


void IrradiationDefectsGenerator::generateIndividualFCC(MicrostructureGenerator& mg)
{
    const std::string defect(this->parser.readString("defect",true));
    if(defect=="SFT")
    {
        std::cout<<magentaBoldColor<<"Generating individual SFTs"<<defaultColor<<std::endl;
        const std::vector<int> sftPlaneIDs(this->parser.readArray<int>("sftPlaneIDs",true));
        const std::vector<int> sftIsInverted(this->parser.readArray<int>("sftIsInverted",true));
        const std::vector<double> sftSizes(this->parser.readArray<double>("sftSizes",true));
        const Eigen::Matrix<double,Eigen::Dynamic,dim> sftBasePoints(this->parser.readMatrix<double>("sftBasePoints",sftPlaneIDs.size(),dim,true));

        if(int(sftPlaneIDs.size())!=sftBasePoints.rows())
        {
            std::cout<<"sftPlaneIDs.size()="<<sftPlaneIDs.size()<<std::endl;
            std::cout<<"sftBasePoints.rows()="<<sftBasePoints.rows()<<std::endl;
            throw std::runtime_error("You must provide one point for each SFT. Each point is a row of the matrix sftBasePoints.");
        }
        
        if(sftPlaneIDs.size()!=sftSizes.size())
        {
            std::cout<<"sftPlaneIDs.size()="<<sftPlaneIDs.size()<<std::endl;
            std::cout<<"sftSizes.size()="<<sftSizes.size()<<std::endl;
            throw std::runtime_error("You must provide one size for each SFT. Each size is an element of the vector sftSizes.");
        }
        
        if(sftPlaneIDs.size()!=sftIsInverted.size())
        {
            std::cout<<"sftPlaneIDs.size()="<<sftPlaneIDs.size()<<std::endl;
            std::cout<<"sftIsInverted()="<<sftSizes.size()<<std::endl;
            throw std::runtime_error("You must provide one boolean value of sftIsInverted for each SFT. Each value is an element of the vector sftIsInverted.");
        }
        
        size_t ndefects=0;
        for(size_t k=0;k<sftPlaneIDs.size();++k)
        {
            const int& planeID(sftPlaneIDs[k]);
            const VectorDimD basePoint(sftBasePoints.row(k));
            const bool inverted(sftIsInverted[k]);
            const int sftSize(std::round(sftSizes[k]/mg.ddBase.poly.b_SI));
            ndefects+=generateSingleSFT(mg,planeID,basePoint,inverted,sftSize);
        }
        std::cout<<"Generated "<<ndefects<<" STFs"<<std::endl;
    }
    else
    {
        throw std::runtime_error("Unkonwn FCC irradiation defect "+defect);
    }
    
}

void IrradiationDefectsGenerator::generateDensityFCC(MicrostructureGenerator& mg)
{
    const std::string defect(this->parser.readString("defect",true));
    if(defect=="SFT")
    {
        std::cout<<magentaBoldColor<<"Generating SFT density"<<defaultColor<<std::endl;
        const double targetSFTdensity(this->parser.readScalar<double>("targetSFTdensity",true));
        const double sftSizeMean(this->parser.readScalar<double>("sftSizeMean",true));
        const double sftSizeStd(this->parser.readScalar<double>("sftSizeStd",true));
        
        std::mt19937 generator;
        std::uniform_int_distribution<> planeNormalDistribution(0,3);
        std::uniform_int_distribution<> invertedThompsonTetrahedronDistribution(0,1);
        std::normal_distribution<double> sftSizeDistribution(sftSizeMean,sftSizeStd);
        
        if(targetSFTdensity>0.0)
        {
            size_t ndefects=0;
            double defectsDensity=ndefects/mg.ddBase.mesh.volume()/std::pow(mg.ddBase.poly.b_SI,3);
            while(defectsDensity<targetSFTdensity)
            {
                const int sftSize(std::round(sftSizeDistribution(generator)/mg.ddBase.poly.b_SI));
                ndefects+=generateSingleSFT(mg,planeNormalDistribution(generator),mg.ddBase.poly.randomPoint(),invertedThompsonTetrahedronDistribution(generator),sftSize);
                defectsDensity=ndefects/mg.ddBase.mesh.volume()/std::pow(mg.ddBase.poly.b_SI,3);
                std::cout<<"SFT density="<<defectsDensity<<std::endl;
            }
        }
    }
    else
    {
        throw std::runtime_error("Unkonwn FCC irradiation defect "+defect);
    }
}

}
#endif
