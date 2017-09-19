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
#include <model/DislocationDynamics/MicrostructureGeneration/MicrostructureGenerator.h>
#include <model/LatticeMath/LatticeVector.h>
#include <model/Mesh/PlaneMeshIntersection.h>
//#include <model/DislocationDynamics/Materials/CrystalOrientation.h>
#include <model/IO/SequentialOutputFile.h>


namespace model
{
    
    class DipolarMicrostructureGenerator : public MicrostructureGenerator
    {
        constexpr static int dim=3;
        typedef Eigen::Matrix<double,dim,1> VectorDimD;
        typedef typename PlaneMeshIntersection<dim>::PlaneMeshIntersectionContainerType PlaneMeshIntersectionContainerType;
        
        std::random_device rd;
        //std::default_random_engine generator;
        std::mt19937 generator;
        
        SequentialOutputFile<'E',1> edgeFile;
        SequentialOutputFile<'V',1> vertexFile;
        SequentialOutputFile<'L',1> loopFile;
        
    public:
        DipolarMicrostructureGenerator() :
        /* init list */ generator(rd())
        //        /* init list */ distribution(0,CrystalOrientation<dim>::slipSystems().size()-1)
        {
            
            EigenDataReader EDR;
            
            double targetDensity=0.0;
            EDR.readScalarInFile("./microstructureInput.txt","targetDensity",targetDensity);
            
            double fractionSessile=0.0;
            EDR.readScalarInFile("./microstructureInput.txt","fractionSessile",fractionSessile);
            
            std::cout<<"Generating dipoles..."<<std::endl;
            double density=0.0;
            size_t nodeID=0;
            size_t loopID=0;
            size_t snID=0;
            
            while(density<targetDensity)
            {
                std::pair<LatticeVector<dim>,int> rp=this->randomPointInMesh();
                const int grainID=rp.second;
                
                const LatticeVector<dim> L0=rp.first;
                //                const LatticeVector<dim> L0(this->poly.grain(grainID));
                //                std::cout<<"CHANGE HERE!!!!!!!"<<std::endl;
                
                
                std::uniform_int_distribution<> distribution(0,this->poly.grain(grainID).slipSystems().size()-1);
                
                
                const int rSS=distribution(generator); // a random SlipSystem
                //                const int rSS=23; // a random SlipSystem
                //                std::cout<<"AND HERE!!!!!!!"<<std::endl;
                
                const SlipSystem& slipSystem=this->poly.grain(grainID).slipSystems()[rSS];
                
                std::set<int> planeIDs;
                for (unsigned int k=0;k<this->poly.grain(grainID).planeNormals().size();++k)
                {
                    if(slipSystem.s.dot(this->poly.grain(grainID).planeNormals()[k])==0)
                    {
                        planeIDs.insert(k);
                    }
                }
                assert(planeIDs.size()==2 && "ONLY FCC IS SUPPORTED AT THE MOMENT.");
                
                
                
                ReciprocalLatticeDirection<3> sr(this->poly.grain(grainID).reciprocalLatticeDirection(slipSystem.s.cartesian()));
                //std::cout<<sr.transpose()<<std::endl;
                
                
                LatticeDirection<3> d1(LatticeVector<dim>(sr.cross(this->poly.grain(grainID).planeNormals()[*planeIDs.begin()])));
                LatticeDirection<3> d2(LatticeVector<dim>(sr.cross(this->poly.grain(grainID).planeNormals()[*planeIDs.rbegin()])));
                LatticeDirection<3> d3(slipSystem.s);
                
                if(density/targetDensity<fractionSessile)
                { // overwrite d2
                    assert(0 && "SESSILE LOOPS NOT SUPPORTED YET.");
                }
                
                
                const double d1cNorm(d1.cartesian().norm());
                const double d2cNorm(d2.cartesian().norm());
                const double d3cNorm(d3.cartesian().norm());
                
                assert(d1cNorm>0.0);
                assert(d2cNorm>0.0);
                assert(d3cNorm>0.0);
                //                std::cout<<d1cNorm<<std::endl;
                //                std::cout<<d2cNorm<<std::endl;
                //
                int a1=this->randomSize()/d1cNorm;
                int a2=this->randomSize()/d2cNorm;
                int a3=this->randomSize()/d3cNorm;
                
                assert(a1!=0);
                assert(a2!=0);
                assert(a3!=0);
                
                LatticeVector<dim> L1=L0+d1*a1;
                LatticeVector<dim> L2=L1+d2*a2;
                LatticeVector<dim> L3=L2-d1*a1;
                
                const VectorDimD P0=L0.cartesian();
                const VectorDimD P1=L1.cartesian();
                const VectorDimD P2=L2.cartesian();
                const VectorDimD P3=L3.cartesian();

                
                //                LatticeVector<dim> L4=L0+d3*a3;
                //                LatticeVector<dim> L5=L1+d3*a3;
                //                LatticeVector<dim> L6=L2+d3*a3;
                //                LatticeVector<dim> L7=L3+d3*a3;
                
                
                if(   mesh.searchRegion(grainID,P1).first
                   && mesh.searchRegion(grainID,P2).first
                   && mesh.searchRegion(grainID,P3).first
                   //                   && mesh.searchRegion(grainID,L4.cartesian()).first
                   //                   && mesh.searchRegion(grainID,L5.cartesian()).first
                   //                   && mesh.searchRegion(grainID,L6.cartesian()).first
                   //                   && mesh.searchRegion(grainID,L7.cartesian()).first
                   )
                {
                    //                    std::cout<<"density="<<density<< "(WARNING: the dislocation length accounts for the part on the boundary)"<<std::endl;
                    
                    const VectorDimD n1=d1.cross(slipSystem.s).cartesian().normalized();
                    const VectorDimD n2=d2.cross(slipSystem.s).cartesian().normalized();
                    
                    
                    
                    PlaneMeshIntersectionContainerType pmi01=PlaneMeshIntersection<dim>(this->mesh).reducedPlaneMeshIntersection(P0,n1,grainID);
                    const VectorDimD P4=this->boundaryProjection(P0,d3.cartesian(),pmi01).second;
                    PlaneMeshIntersectionContainerType pmi12=PlaneMeshIntersection<dim>(this->mesh).reducedPlaneMeshIntersection(P1,n2,grainID);
                    const VectorDimD P5=this->boundaryProjection(P1,d3.cartesian(),pmi12).second;
                    PlaneMeshIntersectionContainerType pmi23=PlaneMeshIntersection<dim>(this->mesh).reducedPlaneMeshIntersection(P2,n1,grainID);
                    const VectorDimD P6=this->boundaryProjection(P2,d3.cartesian(),pmi23).second;
                    PlaneMeshIntersectionContainerType pmi30=PlaneMeshIntersection<dim>(this->mesh).reducedPlaneMeshIntersection(P3,n2,grainID);
                    const VectorDimD P7=this->boundaryProjection(P3,d3.cartesian(),pmi30).second;
                    
                    std::deque<std::pair<int,VectorDimD>> v54=this->boundaryProjection(P1,P0,d3.cartesian(),pmi01);
                    std::deque<std::pair<int,VectorDimD>> v65=this->boundaryProjection(P2,P1,d3.cartesian(),pmi12);
                    std::deque<std::pair<int,VectorDimD>> v76=this->boundaryProjection(P3,P2,d3.cartesian(),pmi23);
                    std::deque<std::pair<int,VectorDimD>> v47=this->boundaryProjection(P0,P3,d3.cartesian(),pmi30);
                    
                    //                    PlaneMeshIntersectionContainerType pmi1=PlaneMeshIntersection<dim>(this->mesh).reducedPlaneMeshIntersection(P1,n1,grainID);
                    //                    const VectorDimD P5=this->boundaryPoint(P1,d3.cartesian(),pmi1);
                    //
                    //                    PlaneMeshIntersectionContainerType pmi2=PlaneMeshIntersection<dim>(this->mesh).reducedPlaneMeshIntersection(P2,n1,grainID);
                    //                    const VectorDimD P6=this->boundaryPoint(P2,d3.cartesian(),pmi2);
                    //                    PlaneMeshIntersectionContainerType pmi3=PlaneMeshIntersection<dim>(this->mesh).reducedPlaneMeshIntersection(P3,n1,grainID);
                    //                    const VectorDimD P7=this->boundaryPoint(P3,d3.cartesian(),pmi3);
                    //
                    //                    std::cout<<"P4="<<P4.transpose()<<std::endl;
                    //                    std::cout<<"P5="<<P5.transpose()<<std::endl;
                    //                                        std::cout<<"P6="<<P6.transpose()<<std::endl;
                    //                                        std::cout<<"P7="<<P7.transpose()<<std::endl;
                    
                    //                    THERE IS A PROBLEM WITH SOME POINTS BEING THE SAME
                    
                    /*! Vertex file format is:
                     * ID Px Py Pz Vx Vy Vz velReducCoeff snID meshLocation grainID
                     */
                    const size_t refNodeID=nodeID;
                    vertexFile << refNodeID+0<<"\t" << std::setprecision(15)<<std::scientific<<P0.transpose()<<"\t"<<Eigen::Matrix<double,1,3>::Zero()<<"\t"<<1.0<<"\t"<< snID <<"\t"<< 0<<"\n";
                    vertexFile << refNodeID+1<<"\t" << std::setprecision(15)<<std::scientific<<P1.transpose()<<"\t"<<Eigen::Matrix<double,1,3>::Zero()<<"\t"<<1.0<<"\t"<< snID <<"\t"<< 0<<"\n";
                    vertexFile << refNodeID+2<<"\t" << std::setprecision(15)<<std::scientific<<P2.transpose()<<"\t"<<Eigen::Matrix<double,1,3>::Zero()<<"\t"<<1.0<<"\t"<< snID <<"\t"<< 0<<"\n";
                    vertexFile << refNodeID+3<<"\t" << std::setprecision(15)<<std::scientific<<P3.transpose()<<"\t"<<Eigen::Matrix<double,1,3>::Zero()<<"\t"<<1.0<<"\t"<< snID <<"\t"<< 0<<"\n";
                    vertexFile << refNodeID+4<<"\t" << std::setprecision(15)<<std::scientific<<P4.transpose()<<"\t"<<Eigen::Matrix<double,1,3>::Zero()<<"\t"<<1.0<<"\t"<< snID <<"\t"<< 0<<"\n";
                    vertexFile << refNodeID+5<<"\t" << std::setprecision(15)<<std::scientific<<P5.transpose()<<"\t"<<Eigen::Matrix<double,1,3>::Zero()<<"\t"<<1.0<<"\t"<< snID <<"\t"<< 0<<"\n";
                    vertexFile << refNodeID+6<<"\t" << std::setprecision(15)<<std::scientific<<P6.transpose()<<"\t"<<Eigen::Matrix<double,1,3>::Zero()<<"\t"<<1.0<<"\t"<< snID <<"\t"<< 0<<"\n";
                    vertexFile << refNodeID+7<<"\t" << std::setprecision(15)<<std::scientific<<P7.transpose()<<"\t"<<Eigen::Matrix<double,1,3>::Zero()<<"\t"<<1.0<<"\t"<< snID <<"\t"<< 0<<"\n";
                    nodeID+=8;
                    
                    
                    writeLoop(P0,0,1,5,v54,4,slipSystem, n1,grainID,refNodeID,nodeID,loopID,snID);
                    writeLoop(P1,1,2,6,v65,5,slipSystem, n2,grainID,refNodeID,nodeID,loopID,snID);
                    writeLoop(P2,2,3,7,v76,6,slipSystem,-n1,grainID,refNodeID,nodeID,loopID,snID);
                    writeLoop(P3,3,0,4,v47,7,slipSystem,-n2,grainID,refNodeID,nodeID,loopID,snID);

                    snID+=1;
                    density += 2.0*(d1cNorm*a1 + d2cNorm*a2)/this->mesh.volume()/std::pow(Material<Isotropic>::b_real,2);
                    
//                    // First loop
//                    loopFile <<loopID+0<<"\t"<< std::setprecision(15)<<std::scientific<<slipSystem.s.cartesian().transpose()<<"\t"<< n1.transpose()<<"\t"<<P0.transpose()<<"\t"<<grainID<<"\n";
//                    edgeFile << loopID+0<<"\t" << refNodeID+0<<"\t"<< refNodeID+1<<"\n";
//                    edgeFile << loopID+0<<"\t" << refNodeID+1<<"\t"<< refNodeID+5<<"\n";
//                    //                    edgeFile << loopID+0<<"\t" << refNodeID+5<<"\t"<< refNodeID+4<<"\n"; // CHANGE THIS
//                    
//                    size_t oldID=refNodeID+5;
//                    for(const auto& pair : v54)
//                    {
//                        vertexFile << nodeID<<"\t" << std::setprecision(15)<<std::scientific<<pair.second.transpose()<<"\t"<<Eigen::Matrix<double,1,3>::Zero()<<"\t"<<1.0<<"\t"<< snID <<"\t"<< 0<<"\n";
//                        edgeFile << loopID+0<<"\t" << oldID<<"\t"<< nodeID<<"\n"; // CHANGE THIS
//                        oldID=nodeID;
//                        nodeID++;
//                    }
//                    edgeFile << loopID+0<<"\t" <<    oldID<<"\t"<< refNodeID+4<<"\n"; // CHANGE THIS
//                    edgeFile << loopID+0<<"\t" << refNodeID+4<<"\t"<< refNodeID+0<<"\n";
//                    
//                    
//                    //                    vertexFile << refNodeID+4<<"\t" << std::setprecision(15)<<std::scientific<<L4.cartesian().transpose()<<"\t"<<Eigen::Matrix<double,1,3>::Zero()<<"\t"<<1.0<<"\t"<< snID <<"\t"<< 0<<"\n";
//                    //                    vertexFile << refNodeID+5<<"\t" << std::setprecision(15)<<std::scientific<<L5.cartesian().transpose()<<"\t"<<Eigen::Matrix<double,1,3>::Zero()<<"\t"<<1.0<<"\t"<< snID <<"\t"<< 0<<"\n";
//                    //                    vertexFile << refNodeID+6<<"\t" << std::setprecision(15)<<std::scientific<<L6.cartesian().transpose()<<"\t"<<Eigen::Matrix<double,1,3>::Zero()<<"\t"<<1.0<<"\t"<< snID <<"\t"<< 0<<"\n";
//                    //                    vertexFile << refNodeID+7<<"\t" << std::setprecision(15)<<std::scientific<<L7.cartesian().transpose()<<"\t"<<Eigen::Matrix<double,1,3>::Zero()<<"\t"<<1.0<<"\t"<< snID <<"\t"<< 0<<"\n";
//                    
//                    /*! Edge file format is:
//                     * loopID sourceID sinkID
//                     */
//                    
//                    edgeFile << loopID+1<<"\t" << refNodeID+1<<"\t"<< refNodeID+2<<"\n";
//                    edgeFile << loopID+1<<"\t" << refNodeID+2<<"\t"<< refNodeID+6<<"\n";
//                    edgeFile << loopID+1<<"\t" << refNodeID+6<<"\t"<< refNodeID+5<<"\n";
//                    edgeFile << loopID+1<<"\t" << refNodeID+5<<"\t"<< refNodeID+1<<"\n";
//                    
//                    edgeFile << loopID+2<<"\t" << refNodeID+2<<"\t"<< refNodeID+3<<"\n";
//                    edgeFile << loopID+2<<"\t" << refNodeID+3<<"\t"<< refNodeID+7<<"\n";
//                    edgeFile << loopID+2<<"\t" << refNodeID+7<<"\t"<< refNodeID+6<<"\n";
//                    edgeFile << loopID+2<<"\t" << refNodeID+6<<"\t"<< refNodeID+2<<"\n";
//                    
//                    edgeFile << loopID+3<<"\t" << refNodeID+3<<"\t"<< refNodeID+0<<"\n";
//                    edgeFile << loopID+3<<"\t" << refNodeID+0<<"\t"<< refNodeID+4<<"\n";
//                    edgeFile << loopID+3<<"\t" << refNodeID+4<<"\t"<< refNodeID+7<<"\n";
//                    edgeFile << loopID+3<<"\t" << refNodeID+7<<"\t"<< refNodeID+3<<"\n";
//                    
//                    
//                    /*! Edge file format is:
//                     * loopID Bx By Bz Nx Ny Nz Lx Ly Lz grainID
//                     * where L is a lattice position in the grain
//                     */
//                    //                    const VectorDimD n1=d1.cross(slipSystem.s).cartesian().normalized();
//                    //                    const VectorDimD n1=d1.cross(slipSystem.s).cartesian().normalized();
//                    
//                    loopFile <<loopID+1<<"\t"<< std::setprecision(15)<<std::scientific<<slipSystem.s.cartesian().transpose()<<"\t"<< n2.transpose()<<"\t"<<P1.transpose()<<"\t"<<grainID<<"\n";
//                    loopFile <<loopID+2<<"\t"<< std::setprecision(15)<<std::scientific<<slipSystem.s.cartesian().transpose()<<"\t"<<-n1.transpose()<<"\t"<<P2.transpose()<<"\t"<<grainID<<"\n";
//                    loopFile <<loopID+3<<"\t"<< std::setprecision(15)<<std::scientific<<slipSystem.s.cartesian().transpose()<<"\t"<<-n2.transpose()<<"\t"<<P3.transpose()<<"\t"<<grainID<<"\n";
//                    
//                    loopID+=4;
//                    snID+=1;
//                    density += 2.0*(d1cNorm*a1 + d2cNorm*a2)/this->mesh.volume()/std::pow(Material<Isotropic>::b_real,2);
                    
                    
                }
            }
            
        }
        

        /**********************************************************************/
        void writeLoop(const VectorDimD& P0,
                       const int& id0,
                       const int& id1,
                       const int& id2,
                       const std::deque<std::pair<int,VectorDimD>>& bndVtx,
                       const int& id3,
                    const SlipSystem& slipSystem,
                       const VectorDimD& n,
                       const int& grainID,
                       const int& refNodeID,
                       size_t& nodeID,
                       size_t& loopID,
                       const size_t& snID)
        {
            // write Loop file
                            loopFile <<loopID+0<<"\t"<< std::setprecision(15)<<std::scientific<<slipSystem.s.cartesian().transpose()<<"\t"<< n.transpose()<<"\t"<<P0.transpose()<<"\t"<<grainID<<"\n";
        
            // write to edge and node files
            edgeFile << loopID<<"\t" << refNodeID+id0<<"\t"<< refNodeID+id1<<"\n";
            edgeFile << loopID<<"\t" << refNodeID+id1<<"\t"<< refNodeID+id2<<"\n";

            size_t oldID=refNodeID+id2;
            for(const auto& pair : bndVtx)
            {
                vertexFile << nodeID<<"\t" << std::setprecision(15)<<std::scientific<<pair.second.transpose()<<"\t"<<Eigen::Matrix<double,1,3>::Zero()<<"\t"<<1.0<<"\t"<< snID <<"\t"<< 0<<"\n";
                edgeFile << loopID<<"\t" << oldID<<"\t"<< nodeID<<"\n"; // CHANGE THIS
                oldID=nodeID;
                nodeID++;
            }
            edgeFile << loopID<<"\t" <<    oldID<<"\t"<< refNodeID+id3<<"\n"; // CHANGE THIS
            edgeFile << loopID<<"\t" << refNodeID+id3<<"\t"<< refNodeID+id0<<"\n";

            loopID+=1;

        }
        
        
    };
    
}
#endif
