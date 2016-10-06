/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_FrankReadMicrostructureGenerator_H_
#define model_FrankReadMicrostructureGenerator_H_

#include <math.h>       /* round, floor, ceil, trunc */
#include <random>
#include <model/DislocationDynamics/MicrostructureGeneration/MicrostructureGenerator.h>
#include <model/LatticeMath/LatticeVector.h>
//#include <model/DislocationDynamics/Materials/CrystalOrientation.h>
#include <model/Utilities/SequentialOutputFile.h>


namespace model
{
    
    class FrankReadMicrostructureGenerator : public MicrostructureGenerator
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
        FrankReadMicrostructureGenerator() :
        /* init list */ generator(rd())
        //        /* init list */ distribution(0,CrystalOrientation<dim>::slipSystems().size()-1)
        {
            //        this->poly.grains().
            EigenDataReader EDR;
            
            double targetDensity=0.0;
            EDR.readScalarInFile("./microstructureInput.txt","targetDensity",targetDensity);
            
            double fractionEdge=0.0;
            EDR.readScalarInFile("./microstructureInput.txt","fractionEdge",fractionEdge);
            
            std::cout<<"Generating Frank-Read sources..."<<std::endl;
            double density=0.0;
            double edgeDensity=0.0;
            
            size_t nodeID=0;
            while(density<targetDensity)
            {
                //std::cout<<"I'm here 1"<<std::endl;
                const std::pair<LatticeVector<dim>,int> L0=this->randomPointInMesh();
                
                std::uniform_int_distribution<> distribution(0,this->poly.grain(L0.second).slipSystems().size()-1);
                
                //std::cout<<"I'm here 2"<<std::endl;
                const int rSS=distribution(generator); // a random SlipSystem
                
                const auto& slipSystem=this->poly.grain(L0.second).slipSystems()[rSS];
                
                //std::cout<<"I'm here 3"<<std::endl;
                ReciprocalLatticeDirection<3> sr(this->poly.grain(L0.second).reciprocalLatticeDirection(slipSystem.s.cartesian()));
                //std::cout<<"I'm here 3a"<<std::endl;
                
                bool isEdge=true;
                
                
                
                LatticeDirection<3> d1(LatticeVector<dim>(sr.cross(slipSystem.n)));
                double d1cNorm(d1.cartesian().norm());
                int a1=this->randomSize()/d1cNorm;
                LatticeVector<dim> L1=L0.first+d1*a1;
                //std::cout<<"I'm here 3b"<<std::endl;

                
                if(edgeDensity>=fractionEdge*density) // overwrite with screw dislocaiton
                {
                    isEdge=false;
                    d1cNorm=slipSystem.s.cartesian().norm();
                    a1=this->randomSize()/d1cNorm;
                    L1=L0.first+slipSystem.s*a1;
                }
                
                //std::cout<<"I'm here 4"<<std::endl;
                
                const auto search1(mesh.search(L1.cartesian()));
                
                //std::cout<<"I'm here 5"<<std::endl;
                if(search1.first && search1.second->region->regionID==L0.second)
                {
                    //std::cout<<"I'm here 6"<<std::endl;
                    density += d1cNorm*a1/this->mesh.volume()/pow(Material<Isotropic>::b_real,2);
                    if(isEdge)
                    {
                        edgeDensity+=d1cNorm*a1/this->mesh.volume()/pow(Material<Isotropic>::b_real,2);
                        std::cout<<"edgeDensity="<<edgeDensity<<std::endl;
                        
                    }
                    std::cout<<"density="<<density<<std::endl;
                    
                    vertexFile << nodeID+0<<"\t" << std::setprecision(15)<<std::scientific<<L0.first.cartesian().transpose()<<"\t" <<VectorDimD::Zero().transpose()<<"\t"<< 0 <<"\t"<< 0<<"\n";
                    vertexFile << nodeID+1<<"\t" << std::setprecision(15)<<std::scientific<<L1.cartesian().transpose()<<"\t" <<VectorDimD::Zero().transpose()<<"\t"<< 0 <<"\t"<< 0<<"\n";
                    //                    vertexFile << nodeID+2<<"\t" << std::setprecision(15)<<std::scientific<<L2.cartesian().transpose()<<"\t" <<VectorDimD::Zero().transpose()<<"\t"<< 0 <<"\t"<< 0<<"\n";
                    //                    vertexFile << nodeID+3<<"\t" << std::setprecision(15)<<std::scientific<<L3.cartesian().transpose()<<"\t" <<VectorDimD::Zero().transpose()<<"\t"<< 0 <<"\t"<< 0<<"\n";
                    
                    
                    edgeFile << nodeID+0<<"\t"<< nodeID+1<<"\t"<< std::setprecision(15)<<std::scientific<<slipSystem.s.cartesian().transpose()<<"\t" <<VectorDimD::Zero().transpose()<<"\t"<< 1.0<<"\t"<< 1.0<<"\t"<< 0<<"\n";
                    //                    edgeFile << nodeID+1<<"\t"<< nodeID+2<<"\t"<< std::setprecision(15)<<std::scientific<<slipSystem.s.cartesian().transpose()<<"\t" <<VectorDimD::Zero().transpose()<<"\t"<< 1.0<<"\t"<< 1.0<<"\t"<< 0<<"\n";
                    //                    edgeFile << nodeID+2<<"\t"<< nodeID+3<<"\t"<< std::setprecision(15)<<std::scientific<<slipSystem.s.cartesian().transpose()<<"\t" <<VectorDimD::Zero().transpose()<<"\t"<< 1.0<<"\t"<< 1.0<<"\t"<< 0<<"\n";
                    //                    edgeFile << nodeID+3<<"\t"<< nodeID+0<<"\t"<< std::setprecision(15)<<std::scientific<<slipSystem.s.cartesian().transpose()<<"\t" <<VectorDimD::Zero().transpose()<<"\t"<< 1.0<<"\t"<< 1.0<<"\t"<< 0<<"\n";
                    
                    nodeID+=2;
                    
                    //std::cout<<"I'm here 7"<<std::endl;
                }
            }
            
        }
        
        
    };
    
}
#endif
