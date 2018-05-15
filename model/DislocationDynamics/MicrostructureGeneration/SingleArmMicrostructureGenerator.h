/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2018 by Yinan Cui <cuiyinan@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SingleArmMicrostructureGenerator_H_
#define model_SingleArmMicrostructureGenerator_H_

#include <math.h>       /* round, floor, ceil, trunc */
#include <random>
#include <model/DislocationDynamics/MicrostructureGeneration/MicrostructureGenerator.h>
#include <model/LatticeMath/LatticeVector.h>
//#include <model/DislocationDynamics/Materials/CrystalOrientation.h>
#include <model/IO/SequentialOutputFile.h>


namespace model
{
    
    class SingleArmMicrostructureGenerator : public MicrostructureGenerator
    {
        constexpr static int dim=3;
        typedef Eigen::Matrix<double,dim,1> VectorDimD;
        typedef typename PlaneMeshIntersection<dim>::PlaneMeshIntersectionContainerType PlaneMeshIntersectionContainerType;
        
        std::random_device rd;
        //std::default_random_engine generator;
        std::mt19937 generator;
        std::uniform_int_distribution<> distribution;
        
      //  SequentialOutputFile<'E',1> edgeFile;
      //  SequentialOutputFile<'V',1> vertexFile;
      //  SequentialOutputFile<'L',1> loopFile;

    public:
        SingleArmMicrostructureGenerator(int argc, char** argv) :
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
            
            std::cout<<"Generating single arm sources..."<<std::endl;
            double density=0.0;
            double edgeDensity=0.0;
            
            size_t nodeID=0;
            size_t loopID=0;
            size_t snID=0;
            

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
                
                LatticeDirection<3> d1(LatticeVector<dim>(sr.cross(slipSystem.n)*this->randomSign()));  // this is the projection direction
           //     if(edgeDensity>fractionEdge*density) // overwrite with screw dislocaiton
          //      {
          //           d1=slipSystem.s;
          //      }
                // Compute the LatticeDireciton corresponding to -n
                LatticeDirection<3> d2(this->poly.grain(grainID).latticeDirection(-slipSystem.n.cartesian()*this->randomSign()));
                double d2cNorm(d2.cartesian().norm());
                const int a2=this->randomSize()/d2cNorm;  
                LatticeVector<dim> L3=L0+d2*a2;
                
               
                const auto search2(mesh.search(L3.cartesian()));
                
         if(  search2.first && search2.second->region->regionID==grainID)
         {
                const VectorDimD P0=L0.cartesian();  
                const VectorDimD P3=L3.cartesian();
                
 
                const VectorDimD n1=slipSystem.n.cartesian().normalized();
                const VectorDimD n2=d2.cross(d1).cartesian().normalized();

            
                PlaneMeshIntersectionContainerType pmi01=PlaneMeshIntersection<dim>(this->mesh,P0,n2,grainID);
                 const VectorDimD P1=this->boundaryProjection(P0,d1.cartesian(),pmi01).second;
                 const VectorDimD P2=this->boundaryProjection(P3,d1.cartesian(),pmi01).second;
                const std::map<double,VectorDimD> P12=this->boundaryProjection(P0,P3,d1.cartesian(),pmi01);
                
               // const VectorDimD P1=P12.begin().second;
               // const VectorDimD P2=P12.rbegin().second;               
                const VectorDimD P4=(P0+P1)/2.0;  
                const VectorDimD P5=(P3+P2)/2.0;

                if ((P1-P0).norm()>a2*0.5 && (P3-P2).norm()>a2*0.5)  //not too small arm 
                {
					density+=((P1-P0).norm()+(P0-P3).norm()+(P3-P2).norm())/this->mesh.volume()/std::pow(Material<Isotropic>::b_real,2);
						/*! Vertex file format is:
						 * ID Px Py Pz Vx Vy Vz velReducCoeff snID meshLocation grainID
						 */
						const size_t refNodeID=nodeID;
						this->nodesIO.emplace_back(refNodeID+0,P0,Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
						this->nodesIO.emplace_back(refNodeID+1,P3,Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
						this->nodesIO.emplace_back(refNodeID+2,P4,Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
						this->nodesIO.emplace_back(refNodeID+3,P5,Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
						this->nodesIO.emplace_back(refNodeID+4,P1,Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
						this->nodesIO.emplace_back(refNodeID+5,P2,Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
						 nodeID+=6;

					if (P12.size()==0)
					{
						this->edgesIO.emplace_back(loopID,refNodeID+4,refNodeID+5,0);
					}
					else if (P12.size()==1)
					{
						this->nodesIO.emplace_back(nodeID,P12.begin()->second,Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
                        this->edgesIO.emplace_back(loopID,refNodeID+4,nodeID,0);
                        this->edgesIO.emplace_back(loopID,nodeID,refNodeID+5,0);
						nodeID++;
					}
					else
					{
						 for(const auto pair : P12)
						{
							  this->nodesIO.emplace_back(nodeID,pair.second,Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
							  nodeID++;
							   if (pair.first==P12.begin()->first)
							   {
                                    this->edgesIO.emplace_back(loopID,refNodeID+4,nodeID,0);
                                    this->edgesIO.emplace_back(loopID,nodeID,nodeID+1,0);								    
							   }
							   else if (pair.first==P12.rbegin()->first )
							   {
                                    this->edgesIO.emplace_back(loopID,nodeID,refNodeID+5,0);								   
							   }
							   else
							   {
								   this->edgesIO.emplace_back(loopID,nodeID,nodeID+1,0);	
							   }
						 }
					 }
	                    this->loopsIO.emplace_back(loopID,slipSystem.s.cartesian(),n2,P0,grainID);
						this->edgesIO.emplace_back(loopID,refNodeID+5,refNodeID+1,0);	
						this->edgesIO.emplace_back(loopID,refNodeID+1,refNodeID+0,0);
						this->edgesIO.emplace_back(loopID,refNodeID+0,refNodeID+4,0);

	                    this->loopsIO.emplace_back(loopID+1,slipSystem.s.cartesian(),n1,P0,grainID);
						this->edgesIO.emplace_back(loopID+1,refNodeID+0,refNodeID+2,0);	
						this->edgesIO.emplace_back(loopID+1,refNodeID+2,refNodeID+4,0);
						this->edgesIO.emplace_back(loopID+1,refNodeID+4,refNodeID+0,0);					

                        this->loopsIO.emplace_back(loopID+2,slipSystem.s.cartesian(),n1,P3,grainID);
						this->edgesIO.emplace_back(loopID+2,refNodeID+1,refNodeID+5,0);	
						this->edgesIO.emplace_back(loopID+2,refNodeID+5,refNodeID+3,0);
						this->edgesIO.emplace_back(loopID+2,refNodeID+3,refNodeID+1,0);	

						loopID+=3;
						snID++;              


	  
					}
			 }
            }
             this->write();
        
	   }
    };
    
}
#endif
