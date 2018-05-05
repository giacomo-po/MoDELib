/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_JunctionMicrostructureGenerator_H_
#define model_JunctionMicrostructureGenerator_H_

#include <math.h>       /* round, floor, ceil, trunc */
#include <random>
#include <algorithm>
#include <model/DislocationDynamics/MicrostructureGeneration/MicrostructureGenerator.h>
#include <model/LatticeMath/LatticeVector.h>
#include <model/Mesh/PlaneMeshIntersection.h>
#include <model/IO/SequentialOutputFile.h>
#include <model/DislocationDynamics/IO/DislocationLoopIO.h>
#include <model/Geometry/SegmentSegmentDistance.h>
#include <model/Geometry/PlanePlaneIntersection.h>


namespace model
{
    
    class JunctionMicrostructureGenerator : public MicrostructureGenerator
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
        
        
        static std::deque<VectorDimD> getNodePos(const SimplicialMesh<dim>& mesh,
                                                 const VectorDimD& P0,
                                                 const VectorDimD& n,
                                                 const int& grainID,
                                                 const VectorDimD& A,
                                                 const VectorDimD& B)
        {
            PlaneMeshIntersection<dim> pmi(mesh,P0,n,grainID);
            
            std::deque<std::pair<VectorDimD,VectorDimD>> segDeq;
            
            for(int k=0;k<pmi.size();++k)
            {
                const int k1=(k+1)<pmi.size()? k+1 :0;
                segDeq.emplace_back(pmi[k].second,pmi[k1].second);
            }
            
            
            std::deque<VectorDimD> nodePos;
            int nIntersections=0;
            for(const auto& pair : segDeq)
            {
                SegmentSegmentDistance<dim> ssi(A,B,pair.first,pair.second);
                
                if(nIntersections==0)
                {
                    if(ssi.dMin>FLT_EPSILON) // no intersection
                    {
                        
                    }
                    else
                    {
                        nIntersections++;
                        nodePos.push_back((ssi.x0+ssi.x1)*0.5);
                    }
                }
                else if(nIntersections==1)
                {
                    if(ssi.dMin>FLT_EPSILON) // no intersection
                    {
                        nodePos.push_back(pair.first);
                    }
                    else
                    {
                        nIntersections++;
                        nodePos.push_back(pair.first);
                        nodePos.push_back((ssi.x0+ssi.x1)*0.5);
                    }
                }
                else
                {
                    break;
                }
            }
            return nodePos;
        }
        
    public:
        JunctionMicrostructureGenerator(int argc, char* argv[]) :
        MicrostructureGenerator(argc,argv),
        /* init list */ generator(rd())
        {
            
            // Read target density
            EigenDataReader EDR;
            double targetDensity=0.0;
            EDR.readScalarInFile("./microstructureInput.txt","targetDensity",targetDensity);
            
            //            double fractionSessile=0.0;
            //            EDR.readScalarInFile("./microstructureInput.txt","fractionSessile",fractionSessile);
            
            
            // init counters
            double density=0.0;
            //            double sessileDensity=0.0;
            size_t nodeID=0;
            size_t loopID=0;
            size_t snID=0;
            
            std::cout<<"Generating junctions..."<<std::endl;
            while(density<targetDensity)
            {
                const std::pair<LatticeVector<dim>,int> rp=this->randomPointInMesh();
                const int& grainID=rp.second;   // random grain ID
                const LatticeVector<dim>& L0=rp.first; // random lattice position in the grain
                const VectorDimD P0(L0.cartesian());   // cartesian position of L0
                std::uniform_int_distribution<> distribution(0,this->poly.grain(grainID).slipSystems().size()-1);
                const int rSS1=distribution(generator); // a random SlipSystem ID
                const int rSS2=distribution(generator); // a random SlipSystem ID
                
                const SlipSystem& slipSystem1=this->poly.grain(grainID).slipSystems()[rSS1];
                const SlipSystem& slipSystem2=this->poly.grain(grainID).slipSystems()[rSS2];
                
                const VectorDimD b1=slipSystem1.s.cartesian();    // Burgers vector
                const VectorDimD b2=slipSystem2.s.cartesian();    // Burgers vector

                
                if(   slipSystem1.n.cross(slipSystem2.n).squaredNorm()>0
                   && slipSystem1.s.cross(slipSystem2.s).squaredNorm()>0
                   && b1.dot(b2)<0.0)
                {
                    const VectorDimD n1=slipSystem1.n.cartesian().normalized(); // slip plane normal
                    const VectorDimD n2=slipSystem2.n.cartesian().normalized(); // slip plane normal
                    
                    
                    PlanePlaneIntersection<dim> ppi(P0,n1,P0,n2);
                    
//                    std::cout<<ppi.d.transpose()<<std::endl;
                    
                    const VectorDimD A=P0+3.0*this->maxSize()*ppi.d;
                    const VectorDimD B=P0-3.0*this->maxSize()*ppi.d;
                    
                    std::deque<VectorDimD> nodePos1=getNodePos(this->mesh,
                                                               P0,
                                                               n1,
                                                               grainID,
                                                               A,
                                                               B);
                    
                    std::deque<VectorDimD> nodePos2=getNodePos(this->mesh,
                                                               P0,
                                                               n2,
                                                               grainID,
                                                               A,
                                                               B);
                    
//                    std::cout<<nodePos1.size()<<" "<<nodePos2.size()<<std::endl;
                    
                    if(nodePos1.size()>3 && nodePos2.size()>3)
                    {
                        
                        if((*nodePos2.begin()-*nodePos1.begin()).squaredNorm()>FLT_EPSILON)
                        {
                            std::reverse(std::begin(nodePos2), std::end(nodePos2));
                        }
                        
                        const VectorDimD J1=*nodePos1.begin()+1.0/3.0*(*nodePos1.rbegin()-*nodePos1.begin());
                        const VectorDimD J2=*nodePos1.begin()+2.0/3.0*(*nodePos1.rbegin()-*nodePos1.begin());
                        
                        nodePos1.pop_front();
                        nodePos1.pop_back();
                        
                        nodePos2.pop_front();
                        nodePos2.pop_back();
                        
                        
                        if(   (*nodePos1.begin()-J1).normalized().dot((*nodePos1.rbegin()-J2).normalized())<0.7172
                           && (*nodePos2.begin()-J1).normalized().dot((*nodePos2.rbegin()-J2).normalized())<0.7172)
                        {
                            const double lineLength=(*nodePos1.begin()-J1).norm()+(*nodePos1.rbegin()-J2).norm()
                            /*                   */+(*nodePos2.begin()-J1).norm()+(*nodePos2.rbegin()-J2).norm()
                            /*                   */+(J1-J2).norm();
                            
                            //(nodePos[nodePos.size()-1]-nodePos[0]).norm();
                            //                nodePos.push_back(nodePos[nodePos.size()-1]+1.0/3.0*(nodePos[0]-nodePos[nodePos.size()-1]));
                            //                nodePos.push_back(nodePos[nodePos.size()-1]+2.0/3.0*(nodePos[0]-nodePos[nodePos.size()-1]));
                            
                            // Write files
                            
                            // write node and edge file
                            int J1ID=-1;
                            int J2ID=-1;
                            
                            for(int k=0;k<nodePos1.size();++k)
                            {
                                if(k<nodePos1.size()-1)
                                {
                                    DislocationNodeIO<dim> dlIO(nodeID+k,nodePos1[k],Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
                                    vertexFile<<dlIO<<"\n";
                                    edgeFile << loopID<<"\t" <<    nodeID+k<<"\t"<< nodeID+k+1<<"\n";
                                }
                                else
                                {
                                    J2ID=nodeID+k+1;
                                    J1ID=nodeID+k+2;
                                    
                                    
                                    // Write last node in nodePos1
                                    vertexFile<<DislocationNodeIO<dim>(nodeID+k,nodePos1[k],Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0)<<"\n";
                                    edgeFile << loopID<<"\t" <<    nodeID+k<<"\t"<< J2ID<<"\n";
                                    
                                    vertexFile<<DislocationNodeIO<dim>(J2ID,J2,Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0)<<"\n";
                                    edgeFile << loopID<<"\t" <<    J2ID<<"\t"<< J1ID<<"\n";
                                    
                                    vertexFile<<DislocationNodeIO<dim>(J1ID,J1,Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0)<<"\n";
                                    edgeFile << loopID<<"\t" <<    J1ID<<"\t"<< nodeID<<"\n";
                                    
                                }
                                
                                //                        vertexFile.write(dlIO);
                                //                        vertexFile << nodeID+k<<"\t" << std::setprecision(15)<<std::scientific<<nodePos[k].transpose()<<"\t"<<Eigen::Matrix<double,1,3>::Zero()<<"\t"<<1.0<<"\t"<< snID <<"\t"<< 0<<"\n";
                                
                                //                            const int nextNodeID=(k+1)<nodePos.size()? nodeID+k+1 : nodeID;
                                //                            edgeFile << loopID<<"\t" <<    nodeID+k<<"\t"<< nextNodeID<<"\n";
                                //
                            }
                            nodeID+=nodePos1.size()+2;
                            
                            // write loop file
                            loopFile<< DislocationLoopIO<dim>(loopID+0, b1,n1,P0,grainID)<<"\n";
                            loopID+=1;
                            
                            
                            for(int k=0;k<nodePos2.size();++k)
                            {
                                if(k<nodePos2.size()-1)
                                {
                                    DislocationNodeIO<dim> dlIO(nodeID+k,nodePos2[k],Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
                                    vertexFile<<dlIO<<"\n";
                                    edgeFile << loopID<<"\t" <<    nodeID+k<<"\t"<< nodeID+k+1<<"\n";
                                }
                                else
                                {
                                    //                                J2ID=nodeID+k+1;
                                    //                                J1ID=nodeID+k+2;
                                    
                                    
                                    // Write last node in nodePos1
                                    vertexFile<<DislocationNodeIO<dim>(nodeID+k,nodePos2[k],Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0)<<"\n";
                                    edgeFile << loopID<<"\t" <<    nodeID+k<<"\t"<< J2ID<<"\n";
                                    
                                    //                                vertexFile<<DislocationNodeIO<dim>(J2ID,J2,Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0)<<"\n";
                                    edgeFile << loopID<<"\t" <<    J2ID<<"\t"<< J1ID<<"\n";
                                    
                                    //                                vertexFile<<DislocationNodeIO<dim>(J1ID,J1,Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0)<<"\n";
                                    edgeFile << loopID<<"\t" <<    J1ID<<"\t"<< nodeID<<"\n";
                                    
                                }
                                
                                //                        vertexFile.write(dlIO);
                                //                        vertexFile << nodeID+k<<"\t" << std::setprecision(15)<<std::scientific<<nodePos[k].transpose()<<"\t"<<Eigen::Matrix<double,1,3>::Zero()<<"\t"<<1.0<<"\t"<< snID <<"\t"<< 0<<"\n";
                                
                                //                            const int nextNodeID=(k+1)<nodePos.size()? nodeID+k+1 : nodeID;
                                //                            edgeFile << loopID<<"\t" <<    nodeID+k<<"\t"<< nextNodeID<<"\n";
                                //
                            }
                            nodeID+=nodePos2.size();
                            
                            loopFile<< DislocationLoopIO<dim>(loopID+0, b2,n2,P0,grainID)<<"\n";
                            
                            
                            loopID+=1;
                            snID+=1;
                            density += lineLength/this->mesh.volume()/std::pow(Material<Isotropic>::b_real,2);
                            //                        if(isSessile)
                            //                        {
                            //                            sessileDensity += lineLength/this->mesh.volume()/std::pow(Material<Isotropic>::b_real,2);
                            //                        }
                            std::cout<<"density="<<density<<std::endl;
                        }
                        

                    }
                    
                    
                }
            }
        }
        
    };
    
}
#endif
