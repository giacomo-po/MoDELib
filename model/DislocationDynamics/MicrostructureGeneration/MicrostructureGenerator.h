/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_MicrostructureGenerator_H_
#define model_MicrostructureGenerator_H_

#include <chrono>
#include <random>
#include <cmath>
#include <assert.h>
#include <Eigen/LU>
#include <Eigen/Cholesky>

#include <model/DislocationDynamics/DislocationNetworkTraits.h>
#include <model/DislocationDynamics/Polycrystals/Polycrystal.h> // defines mode::cout
#include <model/IO/EigenDataReader.h>
#include <model/Mesh/SimplicialMesh.h> // defines mode::cout
#include <model/DislocationDynamics/Materials/Material.h>
#include <model/LatticeMath/LatticeVector.h>
#include <model/DislocationDynamics/DislocationNetwork.h>
#include <model/Mesh/PlaneMeshIntersection.h>
#include <model/DislocationDynamics/IO/DislocationNodeIO.h>
#include <model/DislocationDynamics/IO/DislocationLoopIO.h>
#include <model/DislocationDynamics/IO/DislocationEdgeIO.h>
#include <model/DislocationDynamics/IO/EVLio.h>
#include <model/DislocationDynamics/IO/DislocationLinkingNumber.h>


namespace model
{
    
    
    class MicrostructureGenerator
    {
        constexpr static int dim=3;
        
        typedef Eigen::Matrix<double,dim,1> VectorDimD;
        typedef Eigen::Matrix<long int,dim,1> VectorDimI;
        
        typedef Eigen::Matrix<double,dim,dim>	MatrixDimD;
        typedef Eigen::Matrix<long int,dim,dim>	MatrixDimI;
        
        typedef typename PlaneMeshIntersection<dim>::PlaneMeshIntersectionContainerType PlaneMeshIntersectionContainerType;
        
        typedef DislocationNetwork<dim,1,Hermite> DislocationNetworkType;
        
        std::mt19937 generator;
//        std::uniform_real_distribution<double> distribution;
        //std::uniform_real_distribution<double> sizeDistribution;
        double _minSize;
        double _maxSize;
        
        
        
//        /**********************************************************************/
//        VectorDimD randomPoint()
//        {
//            VectorDimD P0;
//            
//            P0 << mesh.xMin(0)+distribution(generator)*(mesh.xMax(0)-mesh.xMin(0)),
//            /* */ mesh.xMin(1)+distribution(generator)*(mesh.xMax(1)-mesh.xMin(1)),
//            /* */ mesh.xMin(2)+distribution(generator)*(mesh.xMax(2)-mesh.xMin(2));
//            
//            return P0;
//            
//        }
        
        /**********************************************************************/
        static double min(const double& a,const double& b)
        {
            return a<b? a : b;
        }
        
        /**********************************************************************/
        static double max(const double& a,const double& b)
        {
            return a>b? a : b;
        }
        
    public:
        
        size_t nodeID;
        size_t snID;
        size_t loopID;
        
        SimplicialMesh<dim> mesh;
        GlidePlaneObserver<dim> gpo;
        Polycrystal<dim> poly;
        bool outputBinary;
        
        std::vector<DislocationNodeIO<dim>> nodesIO;
        std::vector<DislocationLoopIO<dim>> loopsIO;
        std::vector<DislocationEdgeIO<dim>> edgesIO;
        
        
        std::deque<std::deque<VectorDimD>> loopPoints;
        std::deque<VectorDimD> loopBurgers;
        
        double targetIrradiationLoopDensity;
        double averageLoopSize;
        
        
        /**********************************************************************/
        void addIrradiationLoops()
        {
        
            std::mt19937 generator;

            
            size_t ndefects=0;
            double defectsDensity=ndefects/mesh.volume()/std::pow(Material<Isotropic>::b_real,3);
            while(defectsDensity<targetIrradiationLoopDensity)
            {
                
                
                const std::pair<LatticeVector<dim>,int> rp=this->randomPointInMesh();
                const int& grainID=rp.second;   // random grain ID
                const LatticeVector<dim>& L0=rp.first; // random lattice position in the grain
                const VectorDimD P0(L0.cartesian());   // cartesian position of L0
                
                std::uniform_int_distribution<> distribution(0,poly.grain(grainID).slipSystems().size()-1);
                const int rSS=distribution(generator); // a random SlipSystem ID
                const SlipSystem& slipSystem=this->poly.grain(grainID).slipSystems()[rSS];
                const VectorDimD b=slipSystem.s.cartesian();    // Burgers vector
                const VectorDimD a(b.normalized());
                //const VectorDimD c(slipSystem.s.cartesian());
                
                std::vector<VectorDimD> points;
                
                const int NP=6;
                for(int k=0;k<NP;++k)
                {
                    points.push_back(P0+Eigen::AngleAxis<double>(k*2.0*M_PI/NP,a)*slipSystem.unitNormal*0.5*averageLoopSize/Material<Isotropic>::b_real);
                }
                
                bool pointsIncluded=true;
                for(const auto& point : points)
                {
                    pointsIncluded*=mesh.searchRegion(grainID,point).first;
                }
                
                if(pointsIncluded)
                {
                    
                    for(int k=0;k<NP;++k)
                    {
                        nodesIO.emplace_back(nodeID+k,points[k],Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
                        
                        const int nextNodeID=(k+1)<NP? this->nodeID+k+1 : this->nodeID;
                        edgesIO.emplace_back(this->loopID,this->nodeID+k,nextNodeID,0);

                    }
                    
                    loopsIO.emplace_back(this->loopID+0, b,a,P0,grainID);


                    snID++;
                    loopID++;
                    nodeID+=NP;
                    ndefects++;
                    defectsDensity=ndefects/mesh.volume()/std::pow(Material<Isotropic>::b_real,3);
                    std::cout<<"irradiation defects density="<<defectsDensity<<std::endl;
                }
            }
        }
        
        /**********************************************************************/
        void write()
        {
            
            if(targetIrradiationLoopDensity>0.0)
            {
                addIrradiationLoops();
            }
            
            
            if(outputBinary)
            {
                EVLio<dim>::writeBin(0,nodesIO,loopsIO,edgesIO);
            }
            else
            {
                EVLio<dim>::writeTxt(0,nodesIO,loopsIO,edgesIO);
            }
        }
        
        /**********************************************************************/
        MicrostructureGenerator(int argc, char* argv[]) :
        /* init list */ generator(std::chrono::system_clock::now().time_since_epoch().count()),
//        /* init list */ distribution(0.0,1.0),
        //        /* init list */ sizeDistribution(0.1,0.5),
        /* init list */ _minSize(0.0),
        /* init list */ _maxSize(0.0),
        /* init list */ poly(mesh),
        outputBinary(true),
        nodeID(0),
        snID(0),
        loopID(0)
        {
            int meshID(0);
            EigenDataReader EDR;
            bool use_boundary=false;
            EDR.readScalarInFile("./DDinput.txt","use_boundary",use_boundary);
            
            
            if (use_boundary)
            {
                
                EDR.readScalarInFile("./DDinput.txt","meshID",meshID);
                mesh.readMesh(meshID);
                _minSize=0.1*min(mesh.xMax(0)-mesh.xMin(0),min(mesh.xMax(1)-mesh.xMin(1),mesh.xMax(2)-mesh.xMin(2)));
                _maxSize=max(mesh.xMax(0)-mesh.xMin(0),max(mesh.xMax(1)-mesh.xMin(1),mesh.xMax(2)-mesh.xMin(2)));
                
                poly.init(gpo,"./polyCrystalInput.txt");
                
            }
            else
            {
                assert(0 && "MICROSTRUCTURE GENERATION IN INFINITE SPACE NOT SUPPORTED YET.");
            }
            
            
            unsigned int materialZ;
            EDR.readScalarInFile("./polyCrystalInput.txt","material",materialZ); // material by atomic number Z
            Material<Isotropic>::select(materialZ);
            
            
            EDR.readScalarInFile("./DDinput.txt","outputBinary",outputBinary);

            EDR.readScalarInFile("./microstructureInput.txt","targetIrradiationLoopDensity",targetIrradiationLoopDensity);
            EDR.readScalarInFile("./microstructureInput.txt","averageLoopSize",averageLoopSize);
            
            
        }
        
        
        double deltaHelicity(const std::deque<VectorDimD>& newPoints,
                             const VectorDimD& newBurgers) const
        {
        
            double h(0.0);
            assert(loopPoints.size()==loopBurgers.size());
            for(int k=0;k<loopPoints.size(); ++k)
            {
                h+=LinkingNumber<dim>::loopPairHelicity(loopPoints[k],loopBurgers[k],newPoints,newBurgers);
            }

            return h;
        }
        
        /**********************************************************************/
        std::pair<LatticeVector<dim>,int> randomPointInMesh() const
        {
            return poly.randomLatticePointInMesh();
            
//            VectorDimD P0=randomPoint();
//            auto searchResult=mesh.search(P0);
//            if(searchResult.first)
//            {// point inside
//                LatticeVector<dim> L0 = poly.grain(searchResult.second->region->regionID).snapToLattice(P0);
//                searchResult=mesh.searchRegionWithGuess(L0.cartesian(),searchResult.second);
//                if(searchResult.first)
//                {// point inside
//                    return std::make_pair(L0,searchResult.second->region->regionID);
//                }
//                else
//                {
//                    return randomPointInMesh();
//                }
//            }
//            else
//            {
//                return randomPointInMesh();
//            }
            
        }
        
        /**********************************************************************/
        double randomSize()
        {
            //            return _minSize+distribution(generator)*(_maxSize-_minSize);
            std::uniform_real_distribution<double> dist(_minSize,_maxSize);
            return dist(generator);
            
        }
        
        /**********************************************************************/
        int randomSign()
        {
//            std::random_device rd;  //Will be used to obtain a seed for the random number engine
//            std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
            std::uniform_int_distribution<> dis(0,1);
            return  dis(generator)*2-1;
            
        }
        
        
//        /**********************************************************************/
//        static std::deque<std::pair<int,VectorDimD>> boundaryProjection(const VectorDimD& P0,
//                                                                        const VectorDimD& P1,
//                                                                        const VectorDimD& D,
//                                                                        const PlaneMeshIntersectionContainerType& pp)
//        {
//            
//            std::deque<std::pair<int,VectorDimD>> temp;
//            
//            
//            const double dNorm(D.norm());
//            assert(dNorm>FLT_EPSILON);
//            const VectorDimD dir=D/dNorm;
//            // Let a point v on the boundary be written as v=P0+u1*(P1-P0)+u2*d
//            // then we have [P1-P0 d]*[u1 u2]^T=v-P0
//            
//            Eigen::Matrix<double,3,2> A;
//            A.col(0)=P1-P0;
//            A.col(1)=dir;
//            const Eigen::LLT<Eigen::Matrix<double,2,2>> llt(A.transpose()*A);
//            assert(llt.info()==Eigen::Success);
//            
//            
//            for(size_t m=0;m<pp.size();++m)
//            {
//                const Eigen::Matrix<double,2,1> x=llt.solve(A.transpose()*(pp[m].second-P0));
//                if(x(0)>FLT_EPSILON && x(0)<1.0-FLT_EPSILON && x(1)>FLT_EPSILON)
//                {
//                    temp.emplace_back(m,pp[m].second);
//                }
//            }
//            
//            return temp;
//            
//
//        }
        
        /**********************************************************************/
        static std::map<double,VectorDimD> boundaryProjection(const VectorDimD& P0,
                                                                        const VectorDimD& P1,
                                                                        const VectorDimD& D,
                                                                        const PlaneMeshIntersectionContainerType& pp)
        {
            
            
            
            const double dNorm(D.norm());
            assert(dNorm>FLT_EPSILON);
            const VectorDimD dir=D/dNorm;
            // Let a point v on the boundary be written as v=P0+u1*(P1-P0)+u2*d
            // then we have [P1-P0 d]*[u1 u2]^T=v-P0
            
            Eigen::Matrix<double,3,2> A;
            A.col(0)=P1-P0;
            A.col(1)=dir;
            const Eigen::LLT<Eigen::Matrix<double,2,2>> llt(A.transpose()*A);
            assert(llt.info()==Eigen::Success);
            
            std::map<double,VectorDimD> temp; // keep points sorted by parameter u1
            for(size_t m=0;m<pp.size();++m)
            {
                const Eigen::Matrix<double,2,1> x=llt.solve(A.transpose()*(pp[m].second-P0));
                if(x(0)>FLT_EPSILON && x(0)<1.0-FLT_EPSILON && x(1)>FLT_EPSILON)
                {
                    temp.emplace(x(0),pp[m].second);
                }
            }
            
            return temp;
            
            
        }
        
        /**********************************************************************/
        static std::pair<int,VectorDimD> boundaryProjection(const VectorDimD& P,
                                                            const VectorDimD& D,
                                                            const PlaneMeshIntersectionContainerType& pp)
        {
            const double dNorm(D.norm());
            assert(dNorm>FLT_EPSILON);
            const VectorDimD dir=D/dNorm;
            // line1 is P+u1*dir, u>0
            
            bool success=false;
            std::pair<int,VectorDimD> temp=std::make_pair(-1,VectorDimD::Zero());
            
            for(size_t k=0;k<pp.size();++k)
            {
                const size_t k1 = ((k==(pp.size()-1))? 0 : k+1);
                const VectorDimD& v0=pp[k].second;
                const VectorDimD& v1=pp[k1].second;
                // line2 is v0+u2*(v1-v0), 0<=u2<=1
                
                //P+u1*dir=v0+u2*(v1-v0)
                // [dir -(v1-v0)] [u1 u2] = [v0-P]
                // In least square sense
                // [dir -(v1-v0)]^T*[dir -(v1-v0)]* [u1 u2] = [v0-P]
                
                Eigen::Matrix<double,3,2> A;
                A.col(0)=dir;
                A.col(1)=-(v1-v0);
                
                const Eigen::Matrix<double,3,1> b=v0-P;
                
                const Eigen::LLT<Eigen::Matrix<double,2,2>> llt(A.transpose()*A);
//                std::cout<<"DO NOT USE LLT TO SEE IF SYSTEM HAS SOLUTION. See https://eigen.tuxfamily.org/dox/classEigen_1_1LDLT.html#a858dc77b65dd48248299bb6a6a758abf"<<std::endl;

                
                if(llt.info()==Eigen::Success)
                {
                    
                    const Eigen::Matrix<double,2,1> x=llt.solve(A.transpose()*b);
                    
                    if(x(0)>=0.0 && x(1)>=0.0 && x(1)<=1.0)
                    {
                        success=true;
                        temp=std::make_pair(k,v0+x(1)*(v1-v0));
                        break;
                    }
                }
                
            }
            
            assert(success);
            return temp;
        }
        
        /**********************************************************************/
        const double& minSize()
        {
            return _minSize;
        }
        
        /**********************************************************************/
        const double& maxSize()
        {
            return _maxSize;
        }
        
    };
    
}
#endif
