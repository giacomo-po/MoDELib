/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PeriodicLoopGenerator_cpp_
#define model_PeriodicLoopGenerator_cpp_

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

    void PeriodicLoopGenerator::generateSingle(MicrostructureGenerator& mg,const int& rSS,const VectorDimD& center,const double& radius,const size_t& sides)
    {
        std::pair<bool,const Simplex<dim,dim>*> found(mg.ddBase.mesh.search(center));
        if(!found.first)
        {
            std::cout<<"Point "<<center.transpose()<<" is outside mesh. EXITING."<<std::endl;
            exit(EXIT_FAILURE);
        }
        
        const int grainID(found.second->region->regionID);
        assert(mg.ddBase.poly.grains.size()==1 && "Periodic dislocations only supported for single crystals");
        const auto& grain(mg.ddBase.poly.grain(grainID));
        
        if(rSS>=0 && rSS<int(grain.singleCrystal->slipSystems().size()))
        {
            const auto& slipSystem(*grain.singleCrystal->slipSystems()[rSS]);
            
            
            const long int planeIndex(slipSystem.n.closestPlaneIndexOfPoint(center));
            GlidePlaneKey<3> glidePlaneKey(planeIndex, slipSystem.n);
            std::shared_ptr<PeriodicGlidePlane<3>> glidePlane(mg.ddBase.periodicGlidePlaneFactory.get(glidePlaneKey));
            const VectorDimD P0(glidePlane->referencePlane->snapToPlane(center));
            
            std::vector<VectorDimD> loopNodePos;
            for(size_t k=0;k< sides;++k)
            {
                loopNodePos.push_back(P0+Eigen::AngleAxisd(k*2.0*std::numbers::pi/sides,slipSystem.unitNormal)*slipSystem.s.cartesian().normalized()*radius);
            }
            
            
//            std::map<VectorDimD,size_t,CompareVectorsByComponent<double,dim,float>> uniqueNetworkNodeMap; // networkNodePosition->networkNodeID

            mg.insertJunctionLoop(loopNodePos,glidePlane,
                               slipSystem.s.cartesian(),glidePlane->referencePlane->unitNormal,
                                  P0,grainID,DislocationLoopIO<dim>::GLISSILELOOP);
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

    void PeriodicLoopGenerator::generateDensity(MicrostructureGenerator& mg)
    {
        std::cout<<magentaBoldColor<<"Generating periodic loop density"<<defaultColor<<std::endl;
        const double targetDensity(this->parser.readScalar<double>("targetDensity",true));
        if(targetDensity>0.0)
        {
            const int numberOfSides(this->parser.readScalar<int>("numberOfSides",true));
            const double radiusDistributionMean(this->parser.readScalar<double>("radiusDistributionMean",true));
            const double radiusDistributionStd(this->parser.readScalar<double>("radiusDistributionStd",true));
            std::normal_distribution<double> radiusDistribution(radiusDistributionMean/mg.ddBase.poly.b_SI,radiusDistributionStd/mg.ddBase.poly.b_SI);
            std::mt19937 generator;
            double density=0.0;
            while(density<targetDensity)
            {
                const std::pair<LatticeVector<dim>, int> rp(mg.ddBase.poly.randomLatticePointInMesh());
                const LatticeVector<dim> L0=rp.first;
                const size_t grainID=rp.second;
                std::uniform_int_distribution<> ssDist(0,mg.ddBase.poly.grain(grainID).singleCrystal->slipSystems().size()-1);
                const int rSS(ssDist(generator)); // a random SlipSystem
                const double radius(radiusDistribution(generator));
                try
                {
                    generateSingle(mg,rSS,L0.cartesian(),radius,numberOfSides);
                    density+=2.0*std::numbers::pi*radius/mg.ddBase.mesh.volume()/std::pow(mg.ddBase.poly.b_SI,2);
                    std::cout<<"periodic loop density="<<density<<std::endl;
                }
                catch(const std::exception& e)
                {
                    
                }
            }
        }
    }

    void PeriodicLoopGenerator::generateIndividual(MicrostructureGenerator& mg)
    {
        const std::vector<int> periodicLoopSlipSystemIDs(this->parser.readArray<int>("periodicLoopSlipSystemIDs",true));
        
        if(periodicLoopSlipSystemIDs.size())
        {
            std::cout<<magentaBoldColor<<"Generating individual periodic loops"<<defaultColor<<std::endl;
            const std::vector<double> periodicLoopRadii(this->parser.readArray<double>("periodicLoopRadii_SI",true));
            const Eigen::Matrix<double,Eigen::Dynamic,dim> periodicLoopCenters(this->parser.readMatrix<double>("periodicLoopCenters",periodicLoopSlipSystemIDs.size(),dim,true));
            const std::vector<int> periodicLoopSides(this->parser.readArray<int>("periodicLoopSides",true));
            
            if(periodicLoopSlipSystemIDs.size()!=periodicLoopRadii.size())
            {
                throw std::runtime_error("periodicLoopSlipSystemIDs.size()="+std::to_string(periodicLoopSlipSystemIDs.size())+" NOT EQUAL TO periodicLoopRadii.size()="+std::to_string(periodicLoopRadii.size()));
            }
            if(int(periodicLoopSlipSystemIDs.size())!=periodicLoopCenters.rows())
            {
                throw std::runtime_error("periodicLoopSlipSystemIDs.size()="+std::to_string(periodicLoopSlipSystemIDs.size())+" NOT EQUAL TO periodicLoopCenters.rows()="+std::to_string(periodicLoopCenters.rows()));
            }
            if(periodicLoopSlipSystemIDs.size()!=periodicLoopSides.size())
            {
                throw std::runtime_error("periodicLoopSlipSystemIDs.size()="+std::to_string(periodicLoopSlipSystemIDs.size())+" NOT EQUAL TO periodicLoopSides.size()="+std::to_string(periodicLoopSides.size()));
            }
            
            for(size_t k=0;k<periodicLoopSlipSystemIDs.size();++k)
            {
                generateSingle(mg,periodicLoopSlipSystemIDs[k],periodicLoopCenters.row(k),periodicLoopRadii[k]/mg.ddBase.poly.b_SI,periodicLoopSides[k]);
            }
        }
        
    }

}
#endif
