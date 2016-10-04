/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DipolarMicrostructureGenerator_H_
#define model_DipolarMicrostructureGenerator_H_

#include <math.h>       /* round, floor, ceil, trunc */
#include <random>
#include <vector>
#include <model/DislocationDynamics/MicrostructureGeneration/MicrostructureGenerator.h>
#include <model/LatticeMath/LatticeVector.h>
#include <model/DislocationDynamics/Materials/CrystalOrientation.h>
#include <model/Utilities/SequentialOutputFile.h>


namespace model
{
    
    class DipolarMicrostructureGenerator : public MicrostructureGenerator
    {
        constexpr static int dim=3;
        typedef Eigen::Matrix<double,dim,1> VectorDimD;
        
        std::random_device rd;
        //std::default_random_engine generator;
        std::mt19937 generator;
        std::uniform_int_distribution<> distribution;
        
        SequentialOutputFile<'E',1> edgeFile;
        SequentialOutputFile<'V',1> vertexFile;
        
    public:
        DipolarMicrostructureGenerator() :
        /* init list */ generator(rd()),
        /* init list */ distribution(0,CrystalOrientation<dim>::slipSystems().size()-1)
        {
            
            EigenDataReader EDR;
            
            double targetDensity=0.0;
            EDR.readScalarInFile("./microstructureInput.txt","targetDensity",targetDensity);
            
            double fractionSessile=0.0;
            EDR.readScalarInFile("./microstructureInput.txt","fractionSessile",fractionSessile);
            
            std::cout<<"Generating dipoles..."<<std::endl;
            double density=0.0;
            size_t nodeID=0;
            
            std::vector<int> distVector(CrystalOrientation<dim>::slipSystems().size(),0);
            
            while(density<targetDensity)
            {
                
                LatticeVector<dim> L0=this->randomPointInMesh();
                
                const int rSS=distribution(generator); // a random SlipSystem
                
                const auto& slipSystem=CrystalOrientation<dim>::slipSystems()[rSS];
                
                std::set<int> planeIDs;
                for (unsigned int k=0;k<CrystalOrientation<dim>::planeNormals().size();++k)
                {
                    if(slipSystem.s.dot(CrystalOrientation<dim>::planeNormals()[k])==0)
                    {
                        planeIDs.insert(k);
                    }
                }
                assert(planeIDs.size()==2 && "ONLY FCC IS SUPPORTED AT THE MOMENT.");
                
                LatticeVector<3> v1(LatticeVector<dim>(slipSystem.s.cross(CrystalOrientation<dim>::planeNormals()[*planeIDs.begin()])));
                LatticeVector<3> v2(LatticeVector<dim>(slipSystem.s.cross(CrystalOrientation<dim>::planeNormals()[*planeIDs.rbegin()])));
                if(density<fractionSessile*targetDensity)
                { // overwrite d2
                    for(const auto& slip : CrystalOrientation<dim>::slipSystems())
                    {
                        if(fabs(slip.s.cartesian().dot(slipSystem.s.cartesian()))<FLT_EPSILON)
                        {
                            v2=(LatticeDirection<3>(slip.s));
                            break;
                        }
                    }
                    //                        assert(0 && "SESSILE LOOPS NOT SUPPORTED YET.");
                }
                
                const LatticeDirection<3> d1(v1);
                const LatticeDirection<3> d2(v2);
                
                const double d1cNorm(d1.cartesian().norm());
                const double d2cNorm(d2.cartesian().norm());
                
                int a1=this->randomSize()/d1cNorm;
                int a2=this->randomSize()/d2cNorm;
                
                LatticeVector<dim> L1=L0+d1*a1;
                LatticeVector<dim> L2=L1+d2*a2;
                LatticeVector<dim> L3=L2-d1*a1;
                
                if(   mesh.search(L1.cartesian()).first
                   && mesh.search(L2.cartesian()).first
                   && mesh.search(L3.cartesian()).first)
                {
                    distVector[rSS]++;
                    
                    
                    density += 2.0*(d1cNorm*a1 + d2cNorm*a2)/this->mesh.volume()/pow(Material<Isotropic>::b_real,2);
                    std::cout<<"density="<<density<<std::endl;
                    
                    vertexFile << nodeID+0<<"\t" << std::setprecision(15)<<std::scientific<<L0.cartesian().transpose()<<"\t" <<VectorDimD::Zero().transpose()<<"\t"<< 0 <<"\t"<< 0<<"\n";
                    vertexFile << nodeID+1<<"\t" << std::setprecision(15)<<std::scientific<<L1.cartesian().transpose()<<"\t" <<VectorDimD::Zero().transpose()<<"\t"<< 0 <<"\t"<< 0<<"\n";
                    vertexFile << nodeID+2<<"\t" << std::setprecision(15)<<std::scientific<<L2.cartesian().transpose()<<"\t" <<VectorDimD::Zero().transpose()<<"\t"<< 0 <<"\t"<< 0<<"\n";
                    vertexFile << nodeID+3<<"\t" << std::setprecision(15)<<std::scientific<<L3.cartesian().transpose()<<"\t" <<VectorDimD::Zero().transpose()<<"\t"<< 0 <<"\t"<< 0<<"\n";
                    
                    
                    edgeFile << nodeID+0<<"\t"<< nodeID+1<<"\t"<< std::setprecision(15)<<std::scientific<<slipSystem.s.cartesian().transpose()<<"\t" <<VectorDimD::Zero().transpose()<<"\t"<< 1.0<<"\t"<< 1.0<<"\t"<< 0<<"\n";
                    edgeFile << nodeID+1<<"\t"<< nodeID+2<<"\t"<< std::setprecision(15)<<std::scientific<<slipSystem.s.cartesian().transpose()<<"\t" <<VectorDimD::Zero().transpose()<<"\t"<< 1.0<<"\t"<< 1.0<<"\t"<< 0<<"\n";
                    edgeFile << nodeID+2<<"\t"<< nodeID+3<<"\t"<< std::setprecision(15)<<std::scientific<<slipSystem.s.cartesian().transpose()<<"\t" <<VectorDimD::Zero().transpose()<<"\t"<< 1.0<<"\t"<< 1.0<<"\t"<< 0<<"\n";
                    edgeFile << nodeID+3<<"\t"<< nodeID+0<<"\t"<< std::setprecision(15)<<std::scientific<<slipSystem.s.cartesian().transpose()<<"\t" <<VectorDimD::Zero().transpose()<<"\t"<< 1.0<<"\t"<< 1.0<<"\t"<< 0<<"\n";
                    
                    nodeID+=4;
                    
                }
            }
            
            // Add  "targetDensitySessile"
            if(false)
            {
                std::cout<<"Adding sessile segments"<<std::endl;
                const double targetDensitySessile=10*targetDensity; // maybe readFromFile
                density=0.0;
                while(density<targetDensitySessile)
                {
                    LatticeVector<dim> L0=this->randomPointInMesh();
                    const int rSS=distribution(generator); // a random SlipSystem
                    const auto& slipSystem=CrystalOrientation<dim>::slipSystems()[rSS];
                    
                    LatticeVector<3> v1;
                    for(const auto& slip : CrystalOrientation<dim>::slipSystems())
                    {
                        if(fabs(slip.s.cartesian().dot(slipSystem.s.cartesian()))<FLT_EPSILON)
                        {
                            v1=(LatticeDirection<3>(slip.s));
                            break;
                        }
                    }
                    const LatticeDirection<3> d1(v1);
                    const double d1cNorm(d1.cartesian().norm());
                    int a1=this->randomSize()/d1cNorm;
                    LatticeVector<dim> L1=L0+d1*a1;
                    
                    if(   mesh.search(L1.cartesian()).first)
                    {
                        //distVector[rSS]++;
                        
                        
                        density += d1cNorm*a1/this->mesh.volume()/pow(Material<Isotropic>::b_real,2);
                        std::cout<<"density="<<density<<std::endl;
                        
                        vertexFile << nodeID+0<<"\t" << std::setprecision(15)<<std::scientific<<L0.cartesian().transpose()<<"\t" <<VectorDimD::Zero().transpose()<<"\t"<< 0 <<"\t"<< 0<<"\n";
                        vertexFile << nodeID+1<<"\t" << std::setprecision(15)<<std::scientific<<L1.cartesian().transpose()<<"\t" <<VectorDimD::Zero().transpose()<<"\t"<< 0 <<"\t"<< 0<<"\n";
                        
                        
                        edgeFile << nodeID+0<<"\t"<< nodeID+1<<"\t"<< std::setprecision(15)<<std::scientific<<slipSystem.s.cartesian().transpose()<<"\t" <<VectorDimD::Zero().transpose()<<"\t"<< 1.0<<"\t"<< 1.0<<"\t"<< 0<<"\n";
                        
                        nodeID+=2;
                        
                    }
                    
                }
            }
            
            
            std::cout<<"Slip systems distributions:"<<std::endl;
            for(unsigned int k=0;k<distVector.size();++k)
            {
                std::cout<<"    slip system "<<k<<": "<<distVector[k]<<std::endl;
                
            }
            
        }
        
        
    };
    
}
#endif
