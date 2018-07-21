/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SingleEdgeMicrostructureGenerator_H_
#define model_SingleEdgeMicrostructureGenerator_H_

#include <math.h>       /* round, floor, ceil, trunc */
#include <random>
#include <model/DislocationDynamics/MicrostructureGeneration/MicrostructureGenerator.h>
#include <model/LatticeMath/LatticeVector.h>
#include <model/Mesh/PlaneMeshIntersection.h>
//#include <model/IO/SequentialOutputFile.h>
#include <model/DislocationDynamics/IO/DislocationLoopIO.h>
#include <model/Geometry/SegmentSegmentDistance.h>



namespace model
{
    
    class SingleEdgeMicrostructureGenerator : public MicrostructureGenerator
    {
        constexpr static int dim=3;
        typedef Eigen::Matrix<double,dim,1> VectorDimD;
        typedef typename PlaneMeshIntersection<dim>::PlaneMeshIntersectionContainerType PlaneMeshIntersectionContainerType;
        
        std::random_device rd;
        //std::default_random_engine generator;
        std::mt19937 generator;
        
//        SequentialOutputFile<'E',1> edgeFile;
//        SequentialOutputFile<'V',1> vertexFile;
//        SequentialOutputFile<'L',1> loopFile;
        
        
    public:
        SingleEdgeMicrostructureGenerator(int argc, char* argv[]) :
        MicrostructureGenerator(argc,argv),
        /* init list */ generator(rd())
        {
            
//            // Read target density
//            EigenDataReader EDR;
//            double targetDensity=0.0;
//            EDR.readScalarInFile("./microstructureInput.txt","targetDensity",targetDensity);
//
//            double fractionSessile=0.0;
//            EDR.readScalarInFile("./microstructureInput.txt","fractionSessile",fractionSessile);
//
//            
//            // init counters
//            double density=0.0;
//            double sessileDensity=0.0;
            //size_t nodeID=0;
//            size_t loopID=0;
//            size_t snID=0;
            
            
            std::cout<<"SigleEdgeGenerator assuments generates a single edge dislocation, from the first slip system of the crysta structure (see for example FCCcrystal.h).\n";
            std::cout<<"It assumes that single crystal containing the origin of the reference system.\n";
            
            VectorDimD P0(VectorDimD::Zero());
            
            
            
            std::cout<<"Generating straight line..."<<std::endl;
//            while(density<targetDensity)
//            {
//                const std::pair<LatticeVector<dim>,int> rp=this->randomPointInMesh();
                const int grainID=this->poly.grains().begin()->first;   // random grain ID
//                const LatticeVector<dim>& L0=rp.first; // random lattice position in the grain
//                const VectorDimD P0(L0.cartesian());   // cartesian position of L0
//                std::uniform_int_distribution<> distribution(0,this->poly.grain(grainID).slipSystems().size()-1);
                const int rSS=0; // a random SlipSystem ID
                const SlipSystem& slipSystem=this->poly.grain(grainID).slipSystems()[rSS];
                const VectorDimD b=slipSystem.s.cartesian();    // Burgers vector
                
                
                VectorDimD n=slipSystem.n.cartesian().normalized(); // slip plane normal
                //std::uniform_real_distribution<> dis(0.0, 2.0*M_PI);
            const double theta=0.5*M_PI;//
                VectorDimD d=Eigen::AngleAxisd(theta, n)*b.normalized();

//                bool isSessile=false;
//                if(sessileDensity/targetDensity<fractionSessile)
//                {
//                    n=b.normalized();
//                    isSessile=true;
//                    d=Eigen::AngleAxisd(theta, n)*n.cross(VectorDimD::Random()).normalized();
//                }
                
                
                
                // Define line AB containing dislocaiton and piercing the mesh
                const VectorDimD A=P0+3.0*this->maxSize()*d;
                const VectorDimD B=P0-3.0*this->maxSize()*d;
                
                // Compute interseciton between mesh and glide plane
                PlaneMeshIntersection<dim> pmi(this->mesh,P0,n,grainID);
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
                        else //if(ssi.size==1)
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
                        else //if(ssi.size==1)
                        {
                            nIntersections++;
                            nodePos.push_back(pair.first);
                            nodePos.push_back((ssi.x0+ssi.x1)*0.5);
                        }
                    }
                    else
                    {
                        
                    }
                }
                
                const double lineLength=(nodePos[nodePos.size()-1]-nodePos[0]).norm();
//                nodePos.push_back(nodePos[nodePos.size()-1]+1.0/3.0*(nodePos[0]-nodePos[nodePos.size()-1]));
//                nodePos.push_back(nodePos[nodePos.size()-1]+2.0/3.0*(nodePos[0]-nodePos[nodePos.size()-1]));
                
                // Write files
            assert(nodePos.size()>=3);
//                if()
//                {
                    // write node and edge file
                    for(int k=0;k<nodePos.size();++k)
                    {
                        this->nodesIO.emplace_back(this->nodeID+k,nodePos[k],Eigen::Matrix<double,1,3>::Zero(),1.0,this->snID,0);
//                        DislocationNodeIO<dim> dlIO(nodeID+k,nodePos[k],Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
//                        vertexFile<<dlIO<<"\n";
//                        vertexFile.write(dlIO);
//                        vertexFile << nodeID+k<<"\t" << std::setprecision(15)<<std::scientific<<nodePos[k].transpose()<<"\t"<<Eigen::Matrix<double,1,3>::Zero()<<"\t"<<1.0<<"\t"<< snID <<"\t"<< 0<<"\n";
                        
                        const int nextNodeID=(k+1)<nodePos.size()? this->nodeID+k+1 : this->nodeID;
                        this->edgesIO.emplace_back(this->loopID,this->nodeID+k,nextNodeID,0);
//                        edgeFile << loopID<<"\t" <<    nodeID+k<<"\t"<< nextNodeID<<"\n";
                        
                    }
                    this->nodeID+=nodePos.size();
                    
                    // write loop file
                    this->loopsIO.emplace_back(this->loopID+0, b,n,P0,grainID);
//                    DislocationLoopIO<dim> dlIO(loopID+0, b,n,P0,grainID);
//                    loopFile<< dlIO<<"\n";
                    
                    this->loopID+=1;
                    this->snID+=1;
//                    density += lineLength/this->mesh.volume()/std::pow(Material<Isotropic>::b_real,2);
//                    if(isSessile)
//                    {
//                    sessileDensity += lineLength/this->mesh.volume()/std::pow(Material<Isotropic>::b_real,2);
//                    }
//                    std::cout<<"theta="<<theta*180.0/M_PI<<", density="<<density<<" (sessileDensity="<<sessileDensity<<")"<<std::endl;
//                }
//            }
            
            this->write();
        }
        
    };
    
}
#endif
