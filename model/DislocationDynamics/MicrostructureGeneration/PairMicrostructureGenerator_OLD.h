/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PairMicrostructureGenerator_H_
#define model_PairMicrostructureGenerator_H_

#include <math.h>       /* round, floor, ceil, trunc */
#include <random>
#include <model/DislocationDynamics/MicrostructureGeneration/MicrostructureGenerator.h>
#include <model/LatticeMath/LatticeVector.h>
#include <model/Mesh/PlaneMeshIntersection.h>
#include <model/IO/SequentialOutputFile.h>
#include <model/DislocationDynamics/IO/DislocationLoopIO.h>
#include <model/Geometry/SegmentSegmentDistance.h>


namespace model
{
    
    class PairMicrostructureGenerator : public MicrostructureGenerator
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
        
        // init counters
        double density;
        size_t nodeID;
        size_t loopID;
        size_t snID;

        
        /**********************************************************************/
        std::tuple<std::deque<VectorDimD>,VectorDimD,VectorDimD,VectorDimD,int> getLoop()
        {
            const std::pair<LatticeVector<dim>,int> rp=this->randomPointInMesh();
            const int& grainID=rp.second;   // random grain ID
            const LatticeVector<dim>& L0=rp.first; // random lattice position in the grain
            const VectorDimD P0(L0.cartesian());   // cartesian position of L0
            std::uniform_int_distribution<> distribution(0,this->poly.grain(grainID).slipSystems().size()-1);
            const int rSS=distribution(generator); // a random SlipSystem ID
            const SlipSystem& slipSystem=this->poly.grain(grainID).slipSystems()[rSS];
            const VectorDimD b=slipSystem.s.cartesian();    // Burgers vector
            const VectorDimD n=slipSystem.n.cartesian().normalized(); // slip plane normal
            
            std::uniform_real_distribution<> dis(0.0, 2.0*M_PI);
            const double theta=dis(generator); // random angle of the dislocation line in the plane from screw orientation.
            const VectorDimD d=Eigen::AngleAxisd(theta, n)*b.normalized();
            
            // Define line AB containing the dislocaiton and piercing the mesh
            const VectorDimD A=P0+3.0*this->maxSize()*d;
            const VectorDimD B=P0-3.0*this->maxSize()*d;
            
            // Compute interseciton between mesh and glide plane
            PlaneMeshIntersection<dim> pmi(this->mesh,P0,n,grainID);
            
            // Collect sequence of intersection points into boundary segments
            std::deque<std::pair<VectorDimD,VectorDimD>> segDeq;
            for(int k=0;k<pmi.size();++k)
            {
                const int k1=(k+1)<pmi.size()? k+1 :0;
                segDeq.emplace_back(pmi[k].second,pmi[k1].second);
            }
            
            // Intersect line AB with the boundary segments. Stop at the second intersection to close the loop
            std::deque<VectorDimD> nodePos; // this is the container of dislocaiotn nodes
            int nIntersections=0;
            for(const auto& pair : segDeq)
            {
                SegmentSegmentDistance<dim> ssi(A,B,pair.first,pair.second);
                
                if(nIntersections==0)
                {
                    if(ssi.dMin>FLT_EPSILON)
                    {// no intersection
                        
                    }
                    else
                    {// This is the first intersection with the boundary segments, which becomes a dislocaiotn node
                        nodePos.push_back((ssi.x0+ssi.x1)*0.5);
                        nIntersections++;
                    }
                }
                else if(nIntersections==1)
                {
                    if(ssi.dMin>FLT_EPSILON)
                    {// no intersection. All boundary points after the first intersection become dislocation nodes
                        nodePos.push_back(pair.first);
                        
                    }
                    else
                    {// starting point of the boundary segment and intersection point become nodes
                        nodePos.push_back(pair.first);
                        nodePos.push_back((ssi.x0+ssi.x1)*0.5);
                        nIntersections++;
                    }
                }
                else
                {
                    
                }
            }
            
            const double lineLength=(nodePos[nodePos.size()-1]-nodePos[0]).norm(); // compute line length inside the mesh
            // push back two more nodes inside the mesh along AB
            nodePos.push_back(nodePos[nodePos.size()-1]+1.0/3.0*(nodePos[0]-nodePos[nodePos.size()-1]));
            nodePos.push_back(nodePos[nodePos.size()-2]+2.0/3.0*(nodePos[0]-nodePos[nodePos.size()-2]));

            return std::make_tuple(nodePos,b,n,P0,grainID);
        }
        
        /**********************************************************************/
        void writeLoop(const std::deque<VectorDimD>& nodePos,
                       const VectorDimD& b,
                       const VectorDimD& n,
                       const VectorDimD& P0,
                       const int& grainID)
        {
            
            const double lineLength=(nodePos[nodePos.size()-3]-nodePos[0]).norm(); // compute line length inside the mesh

            
            // write node and edge file
            for(int k=0;k<nodePos.size();++k)
            {
                vertexFile << nodeID+k<<"\t" << std::setprecision(15)<<std::scientific<<nodePos[k].transpose()<<"\t"<<Eigen::Matrix<double,1,3>::Zero()<<"\t"<<1.0<<"\t"<< snID <<"\t"<< 0<<"\n";
                
                const int nextNodeID=(k+1)<nodePos.size()? nodeID+k+1 : nodeID;
                edgeFile << loopID<<"\t" <<    nodeID+k<<"\t"<< nextNodeID<<"\n";
                
            }
            nodeID+=nodePos.size();
            
            // write loop file
            DislocationLoopIO<dim> dlIO(loopID+0, b,n,P0,grainID);
            loopFile<< dlIO<<"\n";
            
            loopID+=1;
            snID+=1;
            density += lineLength/this->mesh.volume()/std::pow(Material<Isotropic>::b_real,2);
//            std::cout<<"theta="<<theta*180.0/M_PI<<", density="<<density<<std::endl;

        }
        
    public:
        
        /**********************************************************************/
        PairMicrostructureGenerator(int argc, char* argv[]) :
        MicrostructureGenerator(argc,argv),
        /* init */ generator(rd()),
        /* init */ density(0.0),
        /* init */ nodeID(0),
        /* init */ loopID(0),
        /* init */ snID(0)
        {
            
            // Read target density
            EigenDataReader EDR;
            double targetDensity=0.0;
            EDR.readScalarInFile("./microstructureInput.txt","targetDensity",targetDensity);
            
            bool searchingPair=true;
            while(searchingPair)
            {
                // Get random Loop 1
                const std::tuple<std::deque<VectorDimD>,VectorDimD,VectorDimD,VectorDimD,int> tup1=getLoop();
                const std::deque<VectorDimD>& nodePos1(std::get<0>(tup1)); // node positions for loop 1
                const VectorDimD& b1(std::get<1>(tup1)); // Burgers vector of loop 1
                const VectorDimD& n1(std::get<2>(tup1)); // plane normal for loop 1
                const VectorDimD& P1(std::get<3>(tup1)); // plane origin for loop 1
                const int& grainID1(std::get<4>(tup1)); // grainID of loop 1
                
                // Get random Loop 2
                const std::tuple<std::deque<VectorDimD>,VectorDimD,VectorDimD,VectorDimD,int> tup2=getLoop();
                const std::deque<VectorDimD>& nodePos2(std::get<0>(tup2)); // node positions for loop 1
                const VectorDimD& b2(std::get<1>(tup2)); // Burgers vector of loop 1
                const VectorDimD& n2(std::get<2>(tup2)); // plane normal for loop 1
                const VectorDimD& P2(std::get<3>(tup2)); // plane origin for loop 1
                const int& grainID2(std::get<4>(tup2)); // grainID of loop 1
                
                
                // Check that loops are attractive
                //            const VectorDimD meshCenter(0.5*(this->mesh.xMax()+this->mesh.xMin()));
                const VectorDimD lineDir1=(nodePos1[nodePos1.size()-3]-nodePos1[0]).normalized();
                const VectorDimD lineDir2=(nodePos2[nodePos2.size()-3]-nodePos2[0]).normalized();
                const bool loopsAreAttractive=(b1.dot(b2)*lineDir1.dot(lineDir2)<0.0);
                if(nodePos1.size()>=5 && nodePos2.size()>=5 && loopsAreAttractive)
                {
                    writeLoop(nodePos1,b1,n1,P1,grainID1);
                    writeLoop(nodePos2,b2,n2,P2,grainID2);
                }
                searchingPair=false;
            }

            
//            std::cout<<"Generating straight lines..."<<std::endl;
//            while(density<targetDensity)
//            {
//                
//                // Write files
//                if(nodePos1.size()>=5 && nodePos2.size()>=5)
//                {
//                    
//                }
//            }
        }
        
    };
    
}
#endif
