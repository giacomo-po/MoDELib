/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_B2BGenerator_H_
#define model_B2BGenerator_H_

#include <math.h>       /* round, floor, ceil, trunc */
#include <random>
#include <model/DislocationDynamics/MicrostructureGeneration/MicrostructureGenerator.h>
#include <model/LatticeMath/LatticeVector.h>
#include <model/DislocationDynamics/Materials/CrystalOrientation.h>
#include <model/Utilities/SequentialOutputFile.h>


namespace model
{
    
    class B2BMicrostructureGenerator : public MicrostructureGenerator
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
        B2BMicrostructureGenerator() :
        /* init list */ generator(rd()),
        /* init list */ distribution(0,CrystalOrientation<dim>::slipSystems().size()-1)
        {
            
            EigenDataReader EDR;
            
            double targetDensity=0.0;
            EDR.readScalarInFile("./microstructureInput.txt","targetDensity",targetDensity);
            
            double fractionEdge=0.0;
            EDR.readScalarInFile("./microstructureInput.txt","fractionEdge",fractionEdge);
            
            std::cout<<"Generating bonudary-t0-boundary dislocatinos..."<<std::endl;
            double density=0.0;
            double edgeDensity=0.0;
            
            size_t nodeID=0;
            while(density<targetDensity)
            {
                
                LatticeVector<dim> L0=this->randomPointInMesh();
                
                const int rSS=distribution(generator); // a random SlipSystem
                
                const auto& slipSystem=CrystalOrientation<dim>::slipSystems()[rSS];
                
                
                
                bool isEdge=true;
                
                LatticeVector<dim> L1=L0;
                LatticeVector<dim> L2=L0;
                
                
                
                if(edgeDensity>=fractionEdge*density)
                {
                    isEdge=false;
                    while(mesh.search(L1.cartesian()).first)
                    {
                        L1+=slipSystem.s; // move in "screw direction"
                    }
                    L1-=slipSystem.s;
                    
                    while(mesh.search(L2.cartesian()).first)
                    {
                        L2-=slipSystem.s;
                    }
                    L2+=slipSystem.s;
                }
                else
                {
                    isEdge=true;
                    // compute "edge direction" d1
                    ReciprocalLatticeDirection<3> sr(slipSystem.s.cartesian());
                    LatticeDirection<3> d1(LatticeVector<dim>(sr.cross(slipSystem.n)));
                    
                    while(mesh.search(L1.cartesian()).first)
                    {
                        L1+=d1;
                    }
                    L1-=d1;
                    
                    while(mesh.search(L2.cartesian()).first)
                    {
                        L2-=d1;
                    }
                    L2+=d1;
                    
                }
                
                
                
                
                assert(mesh.search(L1.cartesian()).first && "L1 outside mesh");
                
                assert(mesh.search(L2.cartesian()).first && "L2 outside mesh");
                
                density += (L1.cartesian()-L2.cartesian()).norm()/this->mesh.volume()/pow(Material<Isotropic>::b_real,2);
                if(isEdge)
                {
                    edgeDensity+=(L1.cartesian()-L2.cartesian()).norm()/this->mesh.volume()/pow(Material<Isotropic>::b_real,2);
                    std::cout<<"edgeDensity="<<edgeDensity<<std::endl;
                    
                }
                std::cout<<"density="<<density<<std::endl;
                
                vertexFile << nodeID+0<<"\t" << std::setprecision(15)<<std::scientific<<L2.cartesian().transpose()<<"\t" <<VectorDimD::Zero().transpose()<<"\t"<< 0 <<"\t"<< 0<<"\n";
                vertexFile << nodeID+1<<"\t" << std::setprecision(15)<<std::scientific<<L0.cartesian().transpose()<<"\t" <<VectorDimD::Zero().transpose()<<"\t"<< 0 <<"\t"<< 0<<"\n";
                vertexFile << nodeID+2<<"\t" << std::setprecision(15)<<std::scientific<<L1.cartesian().transpose()<<"\t" <<VectorDimD::Zero().transpose()<<"\t"<< 0 <<"\t"<< 0<<"\n";

                //                    vertexFile << nodeID+2<<"\t" << std::setprecision(15)<<std::scientific<<L2.cartesian().transpose()<<"\t" <<VectorDimD::Zero().transpose()<<"\t"<< 0 <<"\t"<< 0<<"\n";
                //                    vertexFile << nodeID+3<<"\t" << std::setprecision(15)<<std::scientific<<L3.cartesian().transpose()<<"\t" <<VectorDimD::Zero().transpose()<<"\t"<< 0 <<"\t"<< 0<<"\n";
                
                
                edgeFile << nodeID+0<<"\t"<< nodeID+1<<"\t"<< std::setprecision(15)<<std::scientific<<slipSystem.s.cartesian().transpose()<<"\t" <<VectorDimD::Zero().transpose()<<"\t"<< 1.0<<"\t"<< 1.0<<"\t"<< 0<<"\n";
                edgeFile << nodeID+1<<"\t"<< nodeID+2<<"\t"<< std::setprecision(15)<<std::scientific<<slipSystem.s.cartesian().transpose()<<"\t" <<VectorDimD::Zero().transpose()<<"\t"<< 1.0<<"\t"<< 1.0<<"\t"<< 0<<"\n";

                //                    edgeFile << nodeID+1<<"\t"<< nodeID+2<<"\t"<< std::setprecision(15)<<std::scientific<<slipSystem.s.cartesian().transpose()<<"\t" <<VectorDimD::Zero().transpose()<<"\t"<< 1.0<<"\t"<< 1.0<<"\t"<< 0<<"\n";
                //                    edgeFile << nodeID+2<<"\t"<< nodeID+3<<"\t"<< std::setprecision(15)<<std::scientific<<slipSystem.s.cartesian().transpose()<<"\t" <<VectorDimD::Zero().transpose()<<"\t"<< 1.0<<"\t"<< 1.0<<"\t"<< 0<<"\n";
                //                    edgeFile << nodeID+3<<"\t"<< nodeID+0<<"\t"<< std::setprecision(15)<<std::scientific<<slipSystem.s.cartesian().transpose()<<"\t" <<VectorDimD::Zero().transpose()<<"\t"<< 1.0<<"\t"<< 1.0<<"\t"<< 0<<"\n";
                
                nodeID+=3;
                
            }
            
        }
        
        
    };
    
}
#endif
