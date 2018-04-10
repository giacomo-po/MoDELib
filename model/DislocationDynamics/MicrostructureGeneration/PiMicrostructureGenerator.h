/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PiMicrostructureGenerator_H_
#define model_PiMicrostructureGenerator_H_

#include <math.h>       /* round, floor, ceil, trunc */
#include <random>
#include <model/DislocationDynamics/MicrostructureGeneration/MicrostructureGenerator.h>
#include <model/LatticeMath/LatticeVector.h>
#include <model/Mesh/PlaneMeshIntersection.h>
//#include <model/DislocationDynamics/Materials/CrystalOrientation.h>
#include <model/IO/SequentialOutputFile.h>
#include <model/DislocationDynamics/IO/DislocationLoopIO.h>


namespace model
{
    
    class PiMicrostructureGenerator : public MicrostructureGenerator
    {
        constexpr static int dim=3;
        typedef Eigen::Matrix<double,dim,1> VectorDimD;
        typedef typename PlaneMeshIntersection<dim>::PlaneMeshIntersectionContainerType PlaneMeshIntersectionContainerType;
        
        std::random_device rd;
        std::mt19937 generator;
        
        SequentialOutputFile<'E',1> edgeFile;
        SequentialOutputFile<'V',1> vertexFile;
        SequentialOutputFile<'L',1> loopFile;
        
        
        std::deque<VectorDimD> processLoop(const VectorDimD& b,
                         const VectorDimD& n,
                         const int& grainID,
                         const VectorDimD& P,
                         const VectorDimD& d1,
                         const VectorDimD& d2)
        {
            PlaneMeshIntersectionContainerType pmi=PlaneMeshIntersection<dim>(this->mesh,P,n,grainID);
            std::deque<std::pair<VectorDimD,VectorDimD>> segDeq;
            
            for(int k=0;k<pmi.size();++k)
            {
                const int k1=(k+1)<pmi.size()? k+1 :0;
                segDeq.emplace_back(pmi[k].second,pmi[k1].second);
            }
            
            std::deque<VectorDimD> nodePos;
            int nIntersections=0;
            
            nodePos.push_back(P);
            
            for(const auto& pair : segDeq)
            {
                
                if(nIntersections==0)
                {
                    SegmentSegmentDistance<dim> ssi(P,P+3.0*this->maxSize()*d1,pair.first,pair.second);

                    if(ssi.dMin>FLT_EPSILON) // no intersection
                    {
                        
                    }
                    else //if(ssi.size==1)
                    {
                        nIntersections++;
                        
                        nodePos.push_back((ssi.x0+ssi.x1)*0.5);
                        //                            nodePos.push_back(pair.second);
                    }
                    //                        else
                    //                        {
                    //                            std::cout<<"ssi.size="<<ssi.size<<std::endl;
                    //                            std::cout<<"A="<<A.transpose()<<std::endl;
                    //                                                        std::cout<<"B="<<B.transpose()<<std::endl;
                    //                                                        std::cout<<"pair.first="<<pair.first.transpose()<<std::endl;
                    //                                                        std::cout<<"pair.second="<<pair.second.transpose()<<std::endl;
                    //                            assert(0 && "2 intersections are impossible");
                    //                        }
                }
                else if(nIntersections==1)
                {
                    SegmentSegmentDistance<dim> ssi(P,P+3.0*this->maxSize()*d2,pair.first,pair.second);

                    if(ssi.dMin>FLT_EPSILON) // no intersection
                    {
                        nodePos.push_back(pair.first);
                        
                    }
                    else //if(ssi.size==1)
                    {
                        nIntersections++;
                        nodePos.push_back(pair.first);
                        nodePos.push_back((ssi.x0+ssi.x1)*0.5);
                    }
                    //                        else
                    //                        {
                    //                            std::cout<<"ssi.size="<<ssi.size<<std::endl;
                    //                            std::cout<<"A="<<A.transpose()<<std::endl;
                    //                            std::cout<<"B="<<B.transpose()<<std::endl;
                    //                            std::cout<<"pair.first="<<pair.first.transpose()<<std::endl;
                    //                            std::cout<<"pair.second="<<pair.second.transpose()<<std::endl;
                    //                            assert(0 && "2 intersections are impossible");
                    //                        }
                }
                else
                {
                    
                }
            }
            
            // Write files
//            if(nodePos.size()>=4)
//            {
                // write node and edge file
                for(int k=0;k<nodePos.size();++k)
                {
                    vertexFile << nodeID+k<<"\t" << std::setprecision(15)<<std::scientific<<nodePos[k].transpose()<<"\t"<<Eigen::Matrix<double,1,3>::Zero()<<"\t"<<1.0<<"\t"<< snID <<"\t"<< 0<<"\n";
                    
                    const int nextNodeID=(k+1)<nodePos.size()? nodeID+k+1 : nodeID;
                    edgeFile << loopID<<"\t" <<    nodeID+k<<"\t"<< nextNodeID<<"\n";
                    
                }
                nodeID+=nodePos.size();
                
                // write loop file
                DislocationLoopIO<dim> dlIO(loopID+0, b,n,P,grainID);
                loopFile<< dlIO<<"\n";
                
                loopID+=1;
                snID+=1;
//                density += lineLength/this->mesh.volume()/std::pow(Material<Isotropic>::b_real,2);
//                std::cout<<"theta="<<theta*180.0/M_PI<<", density="<<density<<std::endl;
//            }
            
            return nodePos;
        }
        
        
    public:
        
        double density;
        size_t nodeID;
        size_t loopID;
        size_t snID;
        
        PiMicrostructureGenerator(int argc, char* argv[]) :
        MicrostructureGenerator(argc,argv),
        /* init list */ generator(rd()),
        density(0.0),
        nodeID(0),
        loopID(0),
        snID(0)
        {
            
            EigenDataReader EDR;
            
            double targetDensity=0.0;
            EDR.readScalarInFile("./microstructureInput.txt","targetDensity",targetDensity);
            
            //            double fractionSessile=0.0;
            //            EDR.readScalarInFile("./microstructureInput.txt","fractionSessile",fractionSessile);
            
            std::cout<<"Generating dipoles..."<<std::endl;

            
            while(density<targetDensity)
            {
                std::pair<LatticeVector<dim>,int> rp=this->randomPointInMesh();
                const int grainID=rp.second;
                
                const LatticeVector<dim> L0=rp.first;
                
                std::uniform_int_distribution<> distribution(0,this->poly.grain(grainID).slipSystems().size()-1);
                
                // Random slip system whose Burgers vector define the leg of the pi
                const int rSS=distribution(generator); // a random SlipSystem
                const SlipSystem& slipSystem=this->poly.grain(grainID).slipSystems()[rSS];

                // Find two slip systems which contains slipSystem.s, and sum to a vector orthogonal to slipSystem.s
                std::deque<std::pair<size_t,size_t>> ssIDs;
                for(size_t ss1ID=0;ss1ID<this->poly.grain(grainID).slipSystems().size();++ss1ID)
                {
                    const auto& ss1(this->poly.grain(grainID).slipSystems()[ss1ID]);
                    if(ss1.n.dot(slipSystem.s)==0   // plane of ss1 contains slipSystem.s
                       && ss1.s.cross(slipSystem.s).squaredNorm()!=0) // a different slip system
                    {
                        for(size_t ss2ID=0;ss2ID<this->poly.grain(grainID).slipSystems().size();++ss2ID)
                        {
                            const auto& ss2(this->poly.grain(grainID).slipSystems()[ss2ID]);
                            if(ss2.n.dot(slipSystem.s)==0 // plane of ss2 contains slipSystem.s
                               && ss2.s.cross(slipSystem.s).squaredNorm()!=0 // a different slip system
                               && ss2.s.cross(ss1.s).squaredNorm()!=0 // different from ss1
                               && ss2.n.cross(ss1.n).squaredNorm()!=0) // different from ss1
                            {
                            if(fabs((ss1.s+ss2.s).cartesian().dot(slipSystem.s.cartesian()))<FLT_EPSILON) // ss1.s+ss2.s must sum to a vector orthogonal to slipSystem.s
                                {
                                    ssIDs.emplace_back(ss1ID,ss2ID);
                                }
                            }
                        }
                    }
                }
                if(ssIDs.size()!=16)
                {
                    std::cout<<"ssIDs.size()="<<ssIDs.size()<<std::endl;
                    assert(0 && "ONLY FCC SUPPORTED FOR THIS GENERATOR");
                }
                std::uniform_int_distribution<> distribution1(0,ssIDs.size()-1);
                const int chosenPairID=distribution1(generator); // a random SlipSystem
                const SlipSystem& slipSystem1=this->poly.grain(grainID).slipSystems()[ssIDs[chosenPairID].first];
                const SlipSystem& slipSystem2=this->poly.grain(grainID).slipSystems()[ssIDs[chosenPairID].second];

                
                const VectorDimD P0=L0.cartesian();
                std::deque<VectorDimD> temp0=processLoop(slipSystem1.s.cartesian(),
                            slipSystem1.n.cartesian().normalized(),
                            grainID,
                            P0,
                            slipSystem.s.cartesian().normalized(),
                            slipSystem1.s.cartesian().normalized());
                
                density += (temp0[1]-temp0[0]).norm()/this->mesh.volume()/std::pow(Material<Isotropic>::b_real,2);


                std::deque<VectorDimD> temp1=processLoop(slipSystem2.s.cartesian(),
                                                        slipSystem2.n.cartesian().normalized(),
                                                        grainID,
                                                        P0,
                                                        slipSystem.s.cartesian().normalized(),
                                                        slipSystem2.s.cartesian().normalized());


                
                
                
                
//                // Find slip systems with Burgers vecoctors equal to slipSystem1 and 2, but different plane
//                std::deque<std::pair<size_t,size_t>> ssIDs1;
//                for(size_t ss1ID=0;ss1ID<this->poly.grain(grainID).slipSystems().size();++ss1ID)
//                {
//                    const auto& ss1(this->poly.grain(grainID).slipSystems()[ss1ID]);
//                    if(ss1.n.dot(slipSystem1.s)==0 && ss1.n.cross(slipSystem1.n).squaredNorm()!=0 && (ss1.s-slipSystem1.s).squaredNorm()==0) //
//                    {
////                        std::cout<<"found A"<<std::endl;
//                        for(size_t ss2ID=0;ss2ID<this->poly.grain(grainID).slipSystems().size();++ss2ID)
//                        {
//                            const auto& ss2(this->poly.grain(grainID).slipSystems()[ss2ID]);
//                            if(ss2.n.dot(slipSystem2.s)==0 && ss2.n.cross(slipSystem2.n).squaredNorm()!=0 && (ss2.s-slipSystem2.s).squaredNorm()==0) // a different slip system on same plane of slipSystem
//                            {
////                                                        std::cout<<"found B"<<std::endl;
//                                if(  fabs((ss1.s+ss2.s).cartesian().dot(slipSystem.s.cartesian()))<FLT_EPSILON // ss1.s+ss2.s must sum to a vector orthogonal to slipSystem.s
//                                   && ss2.n.cross(ss1.n).squaredNorm()==0) // ss1 and ss2 must have same normal
//                                {
//                                    ssIDs1.emplace_back(ss1ID,ss2ID);
//                                }
//                            }
//                        }
//                    }
//                }
//                if(ssIDs1.size()!=1)
//                {
//                    std::cout<<"ssIDs1()="<<ssIDs1.size()<<std::endl;
//                    assert(0 && "ONLY FCC SUPPORTED FOR THIS GENERATOR");
//                }
//                const SlipSystem& slipSystem3=this->poly.grain(grainID).slipSystems()[ssIDs1[0].first];
//                const SlipSystem& slipSystem4=this->poly.grain(grainID).slipSystems()[ssIDs1[0].second];
//
//                
//                
//                std::set<int> planeIDs;
//                for (unsigned int k=0;k<this->poly.grain(grainID).planeNormals().size();++k)
//                {
//                    if(slipSystem.s.dot(this->poly.grain(grainID).planeNormals()[k])==0)
//                    {
//                        planeIDs.insert(k);
//                    }
//                }
//                assert(planeIDs.size()==2 && "ONLY FCC IS SUPPORTED AT THE MOMENT.");
//                
//                
//                
//                ReciprocalLatticeDirection<3> sr(this->poly.grain(grainID).reciprocalLatticeDirection(slipSystem.s.cartesian()));
//                
//                LatticeDirection<3> d1(LatticeVector<dim>(sr.cross(this->poly.grain(grainID).planeNormals()[*planeIDs.begin()]))*this->randomSign());
//                LatticeDirection<3> d2(LatticeVector<dim>(sr.cross(this->poly.grain(grainID).planeNormals()[*planeIDs.rbegin()])*this->randomSign()));
//                LatticeDirection<3> d3(slipSystem.s*this->randomSign());
//                
//                //                if(density/targetDensity<fractionSessile)
//                //                { // overwrite d2
//                //                    assert(0 && "SESSILE LOOPS NOT SUPPORTED YET.");
//                //                }
//                
//                
//                const double d1cNorm(d1.cartesian().norm());
//                const double d2cNorm(d2.cartesian().norm());
//                const double d3cNorm(d3.cartesian().norm());
//                
//                assert(d1cNorm>0.0);
//                assert(d2cNorm>0.0);
//                assert(d3cNorm>0.0);
//                
//                int a1=this->randomSize()/d1cNorm;
//                int a2=this->randomSize()/d2cNorm;
//                int a3=this->randomSize()/d3cNorm;
//                
//                assert(a1!=0);
//                assert(a2!=0);
//                assert(a3!=0);
//                
//                LatticeVector<dim> L1=L0+d1*a1;
//                LatticeVector<dim> L2=L1+d2*a2;
//                LatticeVector<dim> L3=L2-d1*a1;
//                
//                const VectorDimD P0=L0.cartesian();
//                const VectorDimD P1=L1.cartesian();
//                const VectorDimD P2=L2.cartesian();
//                const VectorDimD P3=L3.cartesian();
//                
//                
//                if(   mesh.searchRegion(grainID,P1).first
//                   && mesh.searchRegion(grainID,P2).first
//                   && mesh.searchRegion(grainID,P3).first
//                   )
//                {
//                    
//                    const VectorDimD n1=d1.cross(slipSystem.s).cartesian().normalized();
//                    const VectorDimD n2=d2.cross(slipSystem.s).cartesian().normalized();
//                    
//                    
//                    
//                    PlaneMeshIntersectionContainerType pmi01=PlaneMeshIntersection<dim>(this->mesh,P0,n1,grainID);
//                    const VectorDimD P4=this->boundaryProjection(P0,d3.cartesian(),pmi01).second;
//                    PlaneMeshIntersectionContainerType pmi12=PlaneMeshIntersection<dim>(this->mesh,P1,n2,grainID);
//                    const VectorDimD P5=this->boundaryProjection(P1,d3.cartesian(),pmi12).second;
//                    PlaneMeshIntersectionContainerType pmi23=PlaneMeshIntersection<dim>(this->mesh,P2,n1,grainID);
//                    const VectorDimD P6=this->boundaryProjection(P2,d3.cartesian(),pmi23).second;
//                    PlaneMeshIntersectionContainerType pmi30=PlaneMeshIntersection<dim>(this->mesh,P3,n2,grainID);
//                    const VectorDimD P7=this->boundaryProjection(P3,d3.cartesian(),pmi30).second;
//                    
//                    std::deque<std::pair<int,VectorDimD>> v54=this->boundaryProjection(P1,P0,d3.cartesian(),pmi01);
//                    std::deque<std::pair<int,VectorDimD>> v65=this->boundaryProjection(P2,P1,d3.cartesian(),pmi12);
//                    std::deque<std::pair<int,VectorDimD>> v76=this->boundaryProjection(P3,P2,d3.cartesian(),pmi23);
//                    std::deque<std::pair<int,VectorDimD>> v47=this->boundaryProjection(P0,P3,d3.cartesian(),pmi30);
//                    
//                    /*! Vertex file format is:
//                     * ID Px Py Pz Vx Vy Vz velReducCoeff snID meshLocation grainID
//                     */
//                    const size_t refNodeID=nodeID;
//                    vertexFile << refNodeID+0<<"\t" << std::setprecision(15)<<std::scientific<<P0.transpose()<<"\t"<<Eigen::Matrix<double,1,3>::Zero()<<"\t"<<1.0<<"\t"<< snID <<"\t"<< 0<<"\n";
//                    vertexFile << refNodeID+1<<"\t" << std::setprecision(15)<<std::scientific<<P1.transpose()<<"\t"<<Eigen::Matrix<double,1,3>::Zero()<<"\t"<<1.0<<"\t"<< snID <<"\t"<< 0<<"\n";
//                    vertexFile << refNodeID+2<<"\t" << std::setprecision(15)<<std::scientific<<P2.transpose()<<"\t"<<Eigen::Matrix<double,1,3>::Zero()<<"\t"<<1.0<<"\t"<< snID <<"\t"<< 0<<"\n";
//                    vertexFile << refNodeID+3<<"\t" << std::setprecision(15)<<std::scientific<<P3.transpose()<<"\t"<<Eigen::Matrix<double,1,3>::Zero()<<"\t"<<1.0<<"\t"<< snID <<"\t"<< 0<<"\n";
//                    vertexFile << refNodeID+4<<"\t" << std::setprecision(15)<<std::scientific<<P4.transpose()<<"\t"<<Eigen::Matrix<double,1,3>::Zero()<<"\t"<<1.0<<"\t"<< snID <<"\t"<< 0<<"\n";
//                    vertexFile << refNodeID+5<<"\t" << std::setprecision(15)<<std::scientific<<P5.transpose()<<"\t"<<Eigen::Matrix<double,1,3>::Zero()<<"\t"<<1.0<<"\t"<< snID <<"\t"<< 0<<"\n";
//                    vertexFile << refNodeID+6<<"\t" << std::setprecision(15)<<std::scientific<<P6.transpose()<<"\t"<<Eigen::Matrix<double,1,3>::Zero()<<"\t"<<1.0<<"\t"<< snID <<"\t"<< 0<<"\n";
//                    vertexFile << refNodeID+7<<"\t" << std::setprecision(15)<<std::scientific<<P7.transpose()<<"\t"<<Eigen::Matrix<double,1,3>::Zero()<<"\t"<<1.0<<"\t"<< snID <<"\t"<< 0<<"\n";
//                    nodeID+=8;
//                    
//                    
//                    writeLoop(P0,0,1,5,v54,4,slipSystem, n1,grainID,refNodeID,nodeID,loopID,snID);
//                    writeLoop(P1,1,2,6,v65,5,slipSystem, n2,grainID,refNodeID,nodeID,loopID,snID);
//                    writeLoop(P2,2,3,7,v76,6,slipSystem,-n1,grainID,refNodeID,nodeID,loopID,snID);
//                    writeLoop(P3,3,0,4,v47,7,slipSystem,-n2,grainID,refNodeID,nodeID,loopID,snID);
//                    
//                    snID+=1;
//                    density += 2.0*(d1cNorm*a1 + d2cNorm*a2)/this->mesh.volume()/std::pow(Material<Isotropic>::b_real,2);
//                    
//                }
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
            DislocationLoopIO<dim> dlIO(loopID+0, slipSystem.s.cartesian(),n,P0,grainID);
            loopFile<< dlIO<<"\n";
            
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
