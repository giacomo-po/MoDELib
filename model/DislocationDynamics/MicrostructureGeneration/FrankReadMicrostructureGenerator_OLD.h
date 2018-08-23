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
#include <model/IO/SequentialOutputFile.h>


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
        
//        SequentialOutputFile<'E',1> edgeFile;
//        SequentialOutputFile<'V',1> vertexFile;
//        SequentialOutputFile<'L',1> loopFile;

    public:
        FrankReadMicrostructureGenerator(int argc, char** argv) :
        MicrostructureGenerator(argc,argv),
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
            
//            size_t nodeID=0;
//            size_t loopID=0;
//            size_t snID=0;

            while(density<targetDensity)
            {
                const std::pair<LatticeVector<dim>,int> rp=this->randomPointInMesh();
                const LatticeVector<dim> L0=rp.first;
                const int grainID=rp.second;

                std::uniform_int_distribution<> distribution(0,this->poly.grain(grainID).slipSystems().size()-1);
                
                const int rSS=distribution(generator); // a random SlipSystem
                
                const auto& slipSystem=this->poly.grain(grainID).slipSystems()[rSS];
                
                // Compute the ReciprocalLatticeDirection corresponding to s
                ReciprocalLatticeDirection<3> sr(this->poly.grain(grainID).reciprocalLatticeDirection(slipSystem.s.cartesian()));
                
                bool isEdge=true;
                
                
                
                LatticeDirection<3> d1(LatticeVector<dim>(sr.cross(slipSystem.n)*this->randomSign()));
                double d1cNorm(d1.cartesian().norm());
                int a1=this->randomSize()/d1cNorm;
                LatticeVector<dim> L1=L0+d1*a1;

                
                if(edgeDensity>fractionEdge*density) // overwrite with screw dislocaiton
                {
                    isEdge=false;
                    d1cNorm=slipSystem.s.cartesian().norm();
                    a1=this->randomSize()/d1cNorm;
                    L1=L0+slipSystem.s*a1;
                }
                
                // Compute the LatticeDireciton corresponding to -n
                LatticeDirection<3> d2(this->poly.grain(grainID).latticeDirection(-slipSystem.n.cartesian()*this->randomSign()));
                double d2cNorm(d2.cartesian().norm());

                const int a2=2*a1; // aspect ratio of double FR source
                LatticeVector<dim> L2=L1+d2*a2;
                LatticeVector<dim> L3=L0+d2*a2;
                
                
                const auto search1(mesh.search(L1.cartesian()));
                const auto search2(mesh.search(L2.cartesian()));
                const auto search3(mesh.search(L3.cartesian()));
                
                if(   search1.first && search1.second->region->regionID==grainID
                   && search2.first && search2.second->region->regionID==grainID
                   && search3.first && search3.second->region->regionID==grainID)
                {
                    density += 2.0*(d1cNorm*a1 + d2cNorm*a2)/this->mesh.volume()/std::pow(Material<Isotropic>::b_real,2);
                    if(isEdge)
                    {
                        edgeDensity+=2.0*(d1cNorm*a1 + d2cNorm*a2)/this->mesh.volume()/std::pow(Material<Isotropic>::b_real,2);
                        std::cout<<"edgeDensity="<<edgeDensity<<std::endl;
                        
                    }
                    std::cout<<"density="<<density<<std::endl;
                    
                    const VectorDimD P0=L0.cartesian();
                    const VectorDimD P1=L1.cartesian();
                    const VectorDimD P2=L2.cartesian();
                    const VectorDimD P3=L3.cartesian();
                    const VectorDimD P4=0.5*(P0+P1);
                    const VectorDimD P5=0.5*(P2+P3);

                    const VectorDimD n1=slipSystem.n.cartesian().normalized();
                    const VectorDimD n2=d2.cross(d1).cartesian().normalized();

                    this->nodesIO.emplace_back(this->nodeID+0,P0,Eigen::Matrix<double,1,3>::Zero(),1.0,this->snID,0);
                    this->nodesIO.emplace_back(this->nodeID+1,P1,Eigen::Matrix<double,1,3>::Zero(),1.0,this->snID,0);
                    this->nodesIO.emplace_back(this->nodeID+2,P2,Eigen::Matrix<double,1,3>::Zero(),1.0,this->snID,0);
                    this->nodesIO.emplace_back(this->nodeID+3,P3,Eigen::Matrix<double,1,3>::Zero(),1.0,this->snID,0);
                    this->nodesIO.emplace_back(this->nodeID+4,P4,Eigen::Matrix<double,1,3>::Zero(),1.0,this->snID,0);
                    this->nodesIO.emplace_back(this->nodeID+5,P5,Eigen::Matrix<double,1,3>::Zero(),1.0,this->snID,0);

                    
//                    vertexFile << nodeID+0<<"\t" << std::setprecision(15)<<std::scientific<<P0.transpose()<<"\t"<<Eigen::Matrix<double,1,3>::Zero()<<"\t"<<1.0<<"\t"<< snID <<"\t"<< 0<<"\n";
//                    vertexFile << nodeID+1<<"\t" << std::setprecision(15)<<std::scientific<<P1.transpose()<<"\t"<<Eigen::Matrix<double,1,3>::Zero()<<"\t"<<1.0<<"\t"<< snID <<"\t"<< 0<<"\n";
//                    vertexFile << nodeID+2<<"\t" << std::setprecision(15)<<std::scientific<<P2.transpose()<<"\t"<<Eigen::Matrix<double,1,3>::Zero()<<"\t"<<1.0<<"\t"<< snID <<"\t"<< 0<<"\n";
//                    vertexFile << nodeID+3<<"\t" << std::setprecision(15)<<std::scientific<<P3.transpose()<<"\t"<<Eigen::Matrix<double,1,3>::Zero()<<"\t"<<1.0<<"\t"<< snID <<"\t"<< 0<<"\n";
//                    vertexFile << nodeID+4<<"\t" << std::setprecision(15)<<std::scientific<<P4.transpose()<<"\t"<<Eigen::Matrix<double,1,3>::Zero()<<"\t"<<1.0<<"\t"<< snID <<"\t"<< 0<<"\n";
//                    vertexFile << nodeID+5<<"\t" << std::setprecision(15)<<std::scientific<<P5.transpose()<<"\t"<<Eigen::Matrix<double,1,3>::Zero()<<"\t"<<1.0<<"\t"<< snID <<"\t"<< 0<<"\n";
                    
                    this->loopsIO.emplace_back(this->loopID+0,slipSystem.s.cartesian(),n1,P0,grainID);
                    this->loopsIO.emplace_back(this->loopID+1,slipSystem.s.cartesian(),n2,P0,grainID);
                    this->loopsIO.emplace_back(this->loopID+2,slipSystem.s.cartesian(),n1,P3,grainID);

                    this->edgesIO.emplace_back(this->loopID+0,this->nodeID+0,this->nodeID+1,0);
                    this->edgesIO.emplace_back(this->loopID+0,this->nodeID+1,this->nodeID+4,0);
                    this->edgesIO.emplace_back(this->loopID+0,this->nodeID+4,this->nodeID+0,0);

                    this->edgesIO.emplace_back(this->loopID+1,this->nodeID+0,this->nodeID+3,0);
                    this->edgesIO.emplace_back(this->loopID+1,this->nodeID+3,this->nodeID+2,0);
                    this->edgesIO.emplace_back(this->loopID+1,this->nodeID+2,this->nodeID+1,0);
                    this->edgesIO.emplace_back(this->loopID+1,this->nodeID+1,this->nodeID+0,0);
                    
                    this->edgesIO.emplace_back(this->loopID+2,this->nodeID+3,this->nodeID+5,0);
                    this->edgesIO.emplace_back(this->loopID+2,this->nodeID+5,this->nodeID+2,0);
                    this->edgesIO.emplace_back(this->loopID+2,this->nodeID+2,this->nodeID+3,0);

                    
////                    loopFile <<loopID+0<<"\t"<< std::setprecision(15)<<std::scientific<<slipSystem.s.cartesian().transpose()<<"\t"<< n1.transpose()<<"\t"<<P0.transpose()<<"\t"<<grainID<<"\n";
//                    edgeFile << loopID+0<<"\t" << nodeID+0<<"\t"<< nodeID+1<<"\n";
//                    edgeFile << loopID+0<<"\t" << nodeID+1<<"\t"<< nodeID+4<<"\n";
//                    edgeFile << loopID+0<<"\t" << nodeID+4<<"\t"<< nodeID+0<<"\n";
//
////                    loopFile <<loopID+1<<"\t"<< std::setprecision(15)<<std::scientific<<slipSystem.s.cartesian().transpose()<<"\t"<< n2.transpose()<<"\t"<<P0.transpose()<<"\t"<<grainID<<"\n";
//                    edgeFile << loopID+1<<"\t" << nodeID+0<<"\t"<< nodeID+3<<"\n";
//                    edgeFile << loopID+1<<"\t" << nodeID+3<<"\t"<< nodeID+2<<"\n";
//                    edgeFile << loopID+1<<"\t" << nodeID+2<<"\t"<< nodeID+1<<"\n";
//                    edgeFile << loopID+1<<"\t" << nodeID+1<<"\t"<< nodeID+0<<"\n";
//
//                    loopFile <<loopID+2<<"\t"<< std::setprecision(15)<<std::scientific<<slipSystem.s.cartesian().transpose()<<"\t"<< n1.transpose()<<"\t"<<P3.transpose()<<"\t"<<grainID<<"\n";
//                    edgeFile << loopID+2<<"\t" << nodeID+3<<"\t"<< nodeID+5<<"\n";
//                    edgeFile << loopID+2<<"\t" << nodeID+5<<"\t"<< nodeID+2<<"\n";
//                    edgeFile << loopID+2<<"\t" << nodeID+2<<"\t"<< nodeID+3<<"\n";

                    
                    this->nodeID+=6;
                    this->loopID+=3;
                    this->snID++;
                    
                }
            }
            
            this->write();

            
//            // Testing GB dislocations
//            int grainID=1;
//            LatticeVector<dim> L0=this->poly.grainBoundaries().begin()->second.latticePlane(grainID).P;
//            LatticeVector<dim> L1=L0+this->poly.grainBoundaries().begin()->second.latticePlane(grainID).n.primitiveVectors.first*50;
//            LatticeVector<dim> b=this->poly.grainBoundaries().begin()->second.latticePlane(1).n.primitiveVectors.second;
//            vertexFile << nodeID+0<<"\t" << std::setprecision(15)<<std::scientific<<L0.cartesian().transpose()<<"\t" <<VectorDimD::Zero().transpose()<<"\t"<< 0 <<"\t"<< 0<<"\t"<< grainID<<"\n";
//            vertexFile << nodeID+1<<"\t" << std::setprecision(15)<<std::scientific<<L1.cartesian().transpose()<<"\t" <<VectorDimD::Zero().transpose()<<"\t"<< 0 <<"\t"<< 0<<"\n";
//            edgeFile << nodeID+0<<"\t"<< nodeID+1<<"\t"<< std::setprecision(15)<<std::scientific<<b.cartesian().transpose()<<"\t" <<VectorDimD::Zero().transpose()<<"\t"<< 1.0<<"\t"<< 1.0<<"\t"<< 0<<"\n";

        }
        
        
    };
    
}
#endif
