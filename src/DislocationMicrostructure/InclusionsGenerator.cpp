/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_InclusionsGenerator_cpp_
#define model_InclusionsGenerator_cpp_


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
#include <InclusionsGenerator.h>
#include <PlaneLineIntersection.h>

namespace model
{

    InclusionsGenerator::InclusionsGenerator(const std::string& fileName) :
    /* init */ MicrostructureGeneratorBase(fileName)
    /* init */,allowOverlap(false)
    {
        
    }

//    void PeriodicDipoleGenerator::generate(MicrostructureGenerator& mg)
//    {
//        generateIndividual(mg);
//    }


    void InclusionsGenerator::generateDensity(MicrostructureGenerator& mg)
    {
        const double targetInclusionDensity(this->parser.readScalar<double>("targetDensity",true));
        if(targetInclusionDensity>0)
        {
            std::cout<<magentaBoldColor<<"Generating inclusions density"<<defaultColor<<std::endl;
            const double inclusionsDiameterLognormalDistribution_M(this->parser.readScalar<double>("diameterLognormalDistribution_M",true) );
            const double inclusionsDiameterLognormalDistribution_S(this->parser.readScalar<double>("diameterLognormalDistribution_S",true));
            const double inclusionsDiameterLognormalDistribution_A(this->parser.readScalar<double>("diameterLognormalDistribution_A",true));
            const Eigen::Matrix<double,1,dim*dim> inclusionsTransformationEigenDistortion(this->parser.readMatrix<double>("transformationEigenDistortion",1,dim*dim,true));
            const Eigen::Matrix<double,1,dim> patternVector(this->parser.readMatrix<double>("patternVector",1,dim,true)/mg.poly.b_SI);
            const double velocityReductionFactor(this->parser.readScalar<double>("velocityReductionFactor",true));
            const int phaseIDs(this->parser.readScalar<int>("phaseID",true));

            
            std::lognormal_distribution<double> distribution(log(inclusionsDiameterLognormalDistribution_M/inclusionsDiameterLognormalDistribution_A),inclusionsDiameterLognormalDistribution_S);

            const double patternHeigth(patternVector.norm());
            const bool applyPattern(patternHeigth>0.0);
            const VectorDimD patternDir(applyPattern? (patternVector/patternHeigth).eval() : VectorDimD::Zero());

            std::mt19937 generator;
            double density(0.0);
            while(density<targetInclusionDensity)
            {
                
                const double diameter = distribution(generator)*inclusionsDiameterLognormalDistribution_A/mg.poly.b_SI;
                const double radius(0.5*diameter);
                std::pair<LatticeVector<dim>,int> pointPair=mg.poly.randomLatticePointInMesh();
                VectorDimD P=pointPair.first.cartesian();
                const int& grainID(pointPair.second);

                if(applyPattern)
                {
                    const VectorDimD globalVector(mg.poly.grain(grainID).singleCrystal->C2G*patternVector.transpose());
                    const VectorDimD globalDir(mg.poly.grain(grainID).singleCrystal->C2G*patternDir);
                    const long long pointHeigth=std::round(P.dot(globalDir)/patternHeigth);
                    const VectorDimD O(pointHeigth*globalVector);
                    P-=(P-O).dot(globalDir)*globalDir;
                }
                if(generateSingle(mg,P,radius,inclusionsTransformationEigenDistortion,velocityReductionFactor,phaseIDs))
                {
                    density+=1.0/mg.mesh.volume()/std::pow(mg.poly.b_SI,3);
                    std::cout<<"inclusion density="<<density<<std::endl;
                }

                
            }
            
            
            //                            const double patternHeigth(currentPattern.norm());
            //                            const bool applyPattern(patternHeigth>0.0);

            
//            const double inclusionsDiameterLognormalDistribution_S(isInclusionsUsed()? TextFileParser("./inputFiles/initialMicrostructure.txt").readArray<double>("inclusionsDiameterLognormalDistribution_S",true) : std::vector<double>())
//            const double inclusionsDiameterLognormalDistribution_A(isInclusionsUsed()? TextFileParser("./inputFiles/initialMicrostructure.txt").readArray<double>("inclusionsDiameterLognormalDistribution_A",true) : std::vector<double>())
//            const double inclusionsTransformationStrains(isInclusionsUsed()? TextFileParser("./inputFiles/initialMicrostructure.txt").readMatrix<double>("inclusionsTransformationStrains",targetInclusionDensities.size(),MicrostructureGenerator::dim*MicrostructureGenerator::dim,true) : Eigen::Matrix<double,Eigen::Dynamic,MicrostructureGenerator::dim*MicrostructureGenerator::dim>::Zero(1,MicrostructureGenerator::dim*MicrostructureGenerator::dim))
//            const double patternVectors(isInclusionsUsed()? TextFileParser("./inputFiles/initialMicrostructure.txt").readMatrix<double>("patternVectors",targetInclusionDensities.size(),MicrostructureGenerator::dim,true) : Eigen::Matrix<double,Eigen::Dynamic,MicrostructureGenerator::dim>::Zero(1,MicrostructureGenerator::dim))

        }
    }

    void InclusionsGenerator::generateIndividual(MicrostructureGenerator& mg)
    {
        const std::vector<double> inclusionRadii_SI(this->parser.readArray<double>("inclusionRadii_SI",true));
        
        if(inclusionRadii_SI.size())
        {
            std::cout<<magentaBoldColor<<"Generating individual inclusions"<<defaultColor<<std::endl;
//            const std::vector<int> periodicDipoleExitFaceIDs(this->parser.readArray<int>("periodicDipoleExitFaceIDs",true));
            const Eigen::Matrix<double,Eigen::Dynamic,dim> inclusionsCenters(this->parser.readMatrix<double>("inclusionsCenters",inclusionRadii_SI.size(),dim,true));
            const Eigen::Matrix<double,Eigen::Dynamic,dim*dim> inclusionsEigenDistortions(this->parser.readMatrix<double>("inclusionsEigenDistortions",inclusionRadii_SI.size(),dim*dim,true));
            const std::vector<double> inclusionVRF(this->parser.readArray<double>("inclusionVelocityReductionFactors",true));
            const std::vector<int> phaseIDs(this->parser.readArray<int>("phaseIDs",true));

            if(int(inclusionRadii_SI.size())!=inclusionsCenters.rows())
            {
                throw std::runtime_error("inclusionRadii_SI.size()="+std::to_string(inclusionRadii_SI.size())+" NOT EQUAL TO inclusionsCenters.rows()="+std::to_string(inclusionsCenters.rows()));
            }
            if(int(inclusionRadii_SI.size())!=inclusionsEigenDistortions.rows())
            {
                throw std::runtime_error("inclusionRadii_SI.size()="+std::to_string(inclusionRadii_SI.size())+" NOT EQUAL TO inclusionsEigenDistortions.rows()="+std::to_string(inclusionsEigenDistortions.rows()));
            }
            if(inclusionRadii_SI.size()!=inclusionVRF.size())
            {
                throw std::runtime_error("inclusionRadii_SI.size()="+std::to_string(inclusionRadii_SI.size())+" NOT EQUAL TO inclusionVRF.size()="+std::to_string(inclusionVRF.size()));
            }
            if(inclusionRadii_SI.size()!=phaseIDs.size())
            {
                throw std::runtime_error("inclusionRadii_SI.size()="+std::to_string(inclusionRadii_SI.size())+" NOT EQUAL TO phaseIDs.size()="+std::to_string(phaseIDs.size()));
            }
            
            for(size_t k=0;k<inclusionRadii_SI.size();++k)
            {
                if(generateSingle(mg,inclusionsCenters.row(k),inclusionRadii_SI[k]/mg.poly.b_SI,inclusionsEigenDistortions.row(k),inclusionVRF[k],phaseIDs[k]))
                {
                    std::cout<<"generated inclusion"<<std::endl;
                }
            }
        }
        
    }

    bool InclusionsGenerator::generateSingle(MicrostructureGenerator& mg,const VectorDimD& C,const double& R, const Eigen::Matrix<double,1,dim*dim>& eTrow, const double& vrc,const int&type)
    {
        
        Eigen::Matrix<double,dim,dim> eT(Eigen::Map<const Eigen::Matrix<double,dim,dim>>(eTrow.data(),dim,dim).transpose());

        if(   mg.mesh.search(C+VectorDimD::UnitX()*R).first
           && mg.mesh.search(C-VectorDimD::UnitX()*R).first
           && mg.mesh.search(C+VectorDimD::UnitY()*R).first
           && mg.mesh.search(C-VectorDimD::UnitY()*R).first
           && mg.mesh.search(C+VectorDimD::UnitZ()*R).first
           && mg.mesh.search(C-VectorDimD::UnitZ()*R).first)
        {
            if(allowOverlap)
            {
                const size_t inclusionID(mg.insertInclusion(C,R,eT,vrc,type));
                return true;
            }
            else
            {
                bool isOutside(true);
                for(const auto& incl : mg.config().eshelbyInclusions())
                {
                    isOutside=(isOutside && ((C-incl.C).norm()>R+incl.a));
                }
                if(isOutside)
                {
                    const size_t inclusionID(mg.insertInclusion(C,R,eT,vrc,type));
                    return true;
                }
                else
                {
                    return false;
                }
            }
        }
        else
        {
            return false;
        }
        
        
        
        
    }

}
#endif
