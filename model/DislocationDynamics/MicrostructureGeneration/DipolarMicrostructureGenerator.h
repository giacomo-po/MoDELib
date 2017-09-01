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
//#include <model/DislocationDynamics/Materials/CrystalOrientation.h>
#include <model/IO/SequentialOutputFile.h>


namespace model
{
    
    class DipolarMicrostructureGenerator : public MicrostructureGenerator
    {
        constexpr static int dim=3;
        typedef Eigen::Matrix<double,dim,1> VectorDimD;
        
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

                const auto& slipSystem=this->poly.grain(grainID).slipSystems()[rSS];
                
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

                LatticeVector<dim> L4=L0+d3*a3;
                LatticeVector<dim> L5=L1+d3*a3;
                LatticeVector<dim> L6=L2+d3*a3;
                LatticeVector<dim> L7=L3+d3*a3;

                
                if(   mesh.searchRegion(grainID,L1.cartesian()).first
                   && mesh.searchRegion(grainID,L2.cartesian()).first
                   && mesh.searchRegion(grainID,L3.cartesian()).first
                   && mesh.searchRegion(grainID,L4.cartesian()).first
                   && mesh.searchRegion(grainID,L5.cartesian()).first
                   && mesh.searchRegion(grainID,L6.cartesian()).first
                   && mesh.searchRegion(grainID,L7.cartesian()).first)
                {
                    density += 4.0*(d1cNorm*a1 + d2cNorm*a2)/this->mesh.volume()/pow(Material<Isotropic>::b_real,2);
                    std::cout<<"density="<<density<< "(WARNING: the dislocation length accounts for the part on the boundary)"<<std::endl;
                    
                    /*! Vertex file format is:
                     * ID Px Py Pz Vx Vy Vz velReducCoeff snID meshLocation grainID
                     */
                    vertexFile << nodeID+0<<"\t" << std::setprecision(15)<<std::scientific<<L0.cartesian().transpose()<<"\t"<<Eigen::Matrix<double,1,3>::Zero()<<"\t"<<1.0<<"\t"<< snID <<"\t"<< 0<<"\n";
                    vertexFile << nodeID+1<<"\t" << std::setprecision(15)<<std::scientific<<L1.cartesian().transpose()<<"\t"<<Eigen::Matrix<double,1,3>::Zero()<<"\t"<<1.0<<"\t"<< snID <<"\t"<< 0<<"\n";
                    vertexFile << nodeID+2<<"\t" << std::setprecision(15)<<std::scientific<<L2.cartesian().transpose()<<"\t"<<Eigen::Matrix<double,1,3>::Zero()<<"\t"<<1.0<<"\t"<< snID <<"\t"<< 0<<"\n";
                    vertexFile << nodeID+3<<"\t" << std::setprecision(15)<<std::scientific<<L3.cartesian().transpose()<<"\t"<<Eigen::Matrix<double,1,3>::Zero()<<"\t"<<1.0<<"\t"<< snID <<"\t"<< 0<<"\n";
                    vertexFile << nodeID+4<<"\t" << std::setprecision(15)<<std::scientific<<L4.cartesian().transpose()<<"\t"<<Eigen::Matrix<double,1,3>::Zero()<<"\t"<<1.0<<"\t"<< snID <<"\t"<< 0<<"\n";
                    vertexFile << nodeID+5<<"\t" << std::setprecision(15)<<std::scientific<<L5.cartesian().transpose()<<"\t"<<Eigen::Matrix<double,1,3>::Zero()<<"\t"<<1.0<<"\t"<< snID <<"\t"<< 0<<"\n";
                    vertexFile << nodeID+6<<"\t" << std::setprecision(15)<<std::scientific<<L6.cartesian().transpose()<<"\t"<<Eigen::Matrix<double,1,3>::Zero()<<"\t"<<1.0<<"\t"<< snID <<"\t"<< 0<<"\n";
                    vertexFile << nodeID+7<<"\t" << std::setprecision(15)<<std::scientific<<L7.cartesian().transpose()<<"\t"<<Eigen::Matrix<double,1,3>::Zero()<<"\t"<<1.0<<"\t"<< snID <<"\t"<< 0<<"\n";

                    /*! Edge file format is:
                     * loopID sourceID sinkID
                     */
                    edgeFile << loopID+0<<"\t" << nodeID+0<<"\t"<< nodeID+1<<"\n";
                    edgeFile << loopID+0<<"\t" << nodeID+1<<"\t"<< nodeID+5<<"\n";
                    edgeFile << loopID+0<<"\t" << nodeID+5<<"\t"<< nodeID+4<<"\n";
                    edgeFile << loopID+0<<"\t" << nodeID+4<<"\t"<< nodeID+0<<"\n";

                    edgeFile << loopID+1<<"\t" << nodeID+1<<"\t"<< nodeID+2<<"\n";
                    edgeFile << loopID+1<<"\t" << nodeID+2<<"\t"<< nodeID+6<<"\n";
                    edgeFile << loopID+1<<"\t" << nodeID+6<<"\t"<< nodeID+5<<"\n";
                    edgeFile << loopID+1<<"\t" << nodeID+5<<"\t"<< nodeID+1<<"\n";

                    edgeFile << loopID+2<<"\t" << nodeID+2<<"\t"<< nodeID+3<<"\n";
                    edgeFile << loopID+2<<"\t" << nodeID+3<<"\t"<< nodeID+7<<"\n";
                    edgeFile << loopID+2<<"\t" << nodeID+7<<"\t"<< nodeID+6<<"\n";
                    edgeFile << loopID+2<<"\t" << nodeID+6<<"\t"<< nodeID+2<<"\n";

                    edgeFile << loopID+3<<"\t" << nodeID+3<<"\t"<< nodeID+0<<"\n";
                    edgeFile << loopID+3<<"\t" << nodeID+0<<"\t"<< nodeID+4<<"\n";
                    edgeFile << loopID+3<<"\t" << nodeID+4<<"\t"<< nodeID+7<<"\n";
                    edgeFile << loopID+3<<"\t" << nodeID+7<<"\t"<< nodeID+3<<"\n";

                    
                    /*! Edge file format is:
                     * loopID Bx By Bz Nx Ny Nz Lx Ly Lz grainID
                     * where L is a lattice position in the grain
                     */
                    const VectorDimD n1=d1.cross(slipSystem.s).cartesian().normalized();
                    const VectorDimD n2=d2.cross(slipSystem.s).cartesian().normalized();
//                    const VectorDimD n1=d1.cross(slipSystem.s).cartesian().normalized();
//                    const VectorDimD n1=d1.cross(slipSystem.s).cartesian().normalized();

                    loopFile <<loopID+0<<"\t"<< std::setprecision(15)<<std::scientific<<slipSystem.s.cartesian().transpose()<<"\t"<< n1.transpose()<<"\t"<<L0.cartesian().transpose()<<"\t"<<grainID<<"\n";
                    loopFile <<loopID+1<<"\t"<< std::setprecision(15)<<std::scientific<<slipSystem.s.cartesian().transpose()<<"\t"<< n2.transpose()<<"\t"<<L1.cartesian().transpose()<<"\t"<<grainID<<"\n";
                    loopFile <<loopID+2<<"\t"<< std::setprecision(15)<<std::scientific<<slipSystem.s.cartesian().transpose()<<"\t"<<-n1.transpose()<<"\t"<<L2.cartesian().transpose()<<"\t"<<grainID<<"\n";
                    loopFile <<loopID+3<<"\t"<< std::setprecision(15)<<std::scientific<<slipSystem.s.cartesian().transpose()<<"\t"<<-n2.transpose()<<"\t"<<L3.cartesian().transpose()<<"\t"<<grainID<<"\n";
                    
                    nodeID+=8;
                    loopID+=4;
                    snID+=1;
                
                }
            }
        
        }
        
    
    };

}
#endif
