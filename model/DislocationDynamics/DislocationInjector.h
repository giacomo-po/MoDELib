/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationInjector_H_
#define model_DislocationInjector_H_

#include <Eigen/Dense>
#include <LatticeVector.h>
#include <SutherlandHodgman.h>


namespace model
{
    template <int dim>
    class DislocationInjectorBase
    {
        //        constexpr static int dim=DislocationNetworkType::dim;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        //      DislocationNetworkType& DN;
        
    public:
        

        /**********************************************************************/
        static std::vector<VectorDim> straightLineBoundaryClosure(const VectorDim& P0,
                                                                  const VectorDim& d,
                                                                  const VectorDim& n,
                                                                  const int& grainID,
                                                                  const SimplicialMesh<dim>& mesh)
        {
        
            return straightLineBoundaryClosure(P0,d,MeshPlane<dim>(mesh,grainID,P0,n),mesh);
        }
        
        /**********************************************************************/
        static std::vector<VectorDim> straightLineBoundaryClosure(const VectorDim& P0,
                                                                  const VectorDim& d,
                                                                  const MeshPlane<dim>& plane,
                                                                  const SimplicialMesh<dim>& mesh)
        {
            
//            std::cout<<"straightLineBoundaryClosure\n d="<< d.transpose()<<"\n n="<<plane.unitNormal.transpose()<<std::endl;
            
            assert(std::fabs(plane.unitNormal.dot(d))<FLT_EPSILON);
            assert(plane.contains(P0));
            
            const double maxSize=std::max(mesh.xMax(0)-mesh.xMin(0),std::max(mesh.xMax(1)-mesh.xMin(1),mesh.xMax(2)-mesh.xMin(2)));
            
            
            // Define line AB containing dislocaiton and piercing the mesh
            const VectorDim A(P0+3.0*maxSize*d);
            const VectorDim B(P0-3.0*maxSize*d);
            
//            MeshPlane<dim> plane(mesh,grainID,P0,n);
            //            const auto& segDeq(plane.meshIntersections);
            
            std::vector<VectorDim> nodePos;
            int nIntersections=0;
            for(const auto& pair : plane.meshIntersections)
            {
                
                SegmentSegmentDistance<dim> ssi(A,B,pair->P0,pair->P1);
                
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
                    if(ssi.dMin>FLT_EPSILON)
                    {// no intersection
                        if((pair->P0-nodePos.back()).norm()>FLT_EPSILON)
                        {
                            nodePos.push_back(pair->P0);
                        }
                    }
                    else
                    {// intersection
                        const VectorDim X((ssi.x0+ssi.x1)*0.5);
                        if((X-nodePos.back()).norm()>FLT_EPSILON)
                        {
                            nIntersections++;
                            nodePos.push_back(pair->P0);
                            if((X-pair->P0).norm()>FLT_EPSILON)
                            {
                                nodePos.push_back(X);
                            }
                        }
                    }
                }
                else
                {
                    
                }
            }
            
            if(nodePos.size()<3)
            {
                std::cout<<"nodePos.size="<<nodePos.size()<<std::endl;
                std::cout<<"plane.meshIntersections are:"<<std::endl;
                std::cout<<std::setprecision(15)<<std::scientific<<plane.meshIntersections<<std::endl;
                std::cout<<"P0="<<std::setprecision(15)<<std::scientific<<P0.transpose()<<std::endl;
                std::cout<<"d="<<std::setprecision(15)<<std::scientific<<d.transpose()<<std::endl;
                assert(false && "LOOP DOES NOT HAVE ENOUGH POINTS");
            }
            
            if(isRightHandedPolygon(nodePos,plane.unitNormal))
            {
//                std::cout<<"straightLineBoundaryClosure: right-handed"<<std::endl;
                return nodePos;
            }
            else
            {
//                std::cout<<"straightLineBoundaryClosure: NOT right-handed"<<std::endl;
                std::vector<VectorDim> revNodePos;
                for(typename std::vector<VectorDim>::reverse_iterator rIter=nodePos.rbegin();rIter!=nodePos.rend();++rIter)
                {
                    revNodePos.push_back(*rIter);
                }
                return revNodePos;
            }            
        }
        
        /**********************************************************************/
        static bool isRightHandedPolygon(const std::vector<VectorDim>& poly,const VectorDim& n)
        {
            if(poly.size())
            {
                VectorDim nA(VectorDim::Zero());
                const VectorDim P0(poly.front());
                for(size_t k=0;k<poly.size();++k)
                {
                    size_t k1(k==poly.size()-1? 0 : k+1);
                    nA+= 0.5*(poly[k]-P0).cross(poly[k1]-poly[k]);
                }
                return nA.dot(n)>0.0;
            }
            else
            {
                return false;
            }
        }
        
        
        /**********************************************************************/
//        template <typename DislocationNetworkType>
//        void insertRandomStraightDislocation(DislocationNetworkType& DN)
//        {
//            typedef typename DislocationNetworkType::NodeType NodeType;
//            const std::pair<LatticeVector<dim>,int> rp=DN.poly.randomLatticePointInMesh();
//            const int& grainID=rp.second;   // random grain ID
//            const LatticeVector<dim>& L0=rp.first; // random lattice position in the grain
//            const VectorDim P0(L0.cartesian());   // cartesian position of L0
//            std::mt19937 generator(std::chrono::system_clock::now().time_since_epoch().count());
//            std::uniform_int_distribution<> distribution(0,DN.poly.grain(grainID).slipSystems().size()-1);
//            const int rSS=distribution(generator); // a random SlipSystem ID
//            const SlipSystem& slipSystem=DN.poly.grain(grainID).slipSystems()[rSS];
//            const VectorDim b=slipSystem.s.cartesian();    // Burgers vector
//            const VectorDim n=slipSystem.n.cartesian().normalized(); // slip plane normal
//
//            std::uniform_real_distribution<> dis(0.0, 2.0*M_PI);
//            const double theta=dis(generator); // random angle of the dislocation line in the plane from screw orientation.
//            const VectorDim d=Eigen::AngleAxisd(theta,n)*b.normalized(); // a random direction in the glide plane
//
//            std::vector<VectorDim> nodePos(straightLineBoundaryClosure(P0,d,n,grainID,DN.mesh));
//            std::vector<std::shared_ptr<NodeType>> loopNodes;
//            for (const auto& pos : nodePos)
//            {
//                loopNodes.emplace_back(new NodeType(&DN,pos,VectorDim::Zero(),1.0));
//            }
//            DN.insertLoop(loopNodes,b,n,P0,grainID);
////            DN.clearDanglingNodes();
//        }
        
    };
    
    
    template <typename DislocationNetworkType>
    class DislocationInjector : public DislocationInjectorBase<DislocationNetworkType::dim>
    {
        static constexpr int dim=DislocationNetworkType::dim;
        typedef typename DislocationNetworkType::NodeType NodeType;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        
        typedef std::tuple<const Simplex<dim,dim-1>* const,const SlipSystem* const,double,size_t> SurfaceNucleationTupleType; // face, slipSystem, radius, sides

        
        DislocationNetworkType& DN;
        
    public:
        
        const int surfaceDislocationNucleationModel;
        const double criticalSurfaceDislocationNucleationShearStress;

        DislocationInjector(DislocationNetworkType& DN_in) :
        /* init */ DN(DN_in)
        /* init */,surfaceDislocationNucleationModel(TextFileParser("inputFiles/DD.txt").readScalar<int>("surfaceDislocationNucleationModel",true))
        /* init */,criticalSurfaceDislocationNucleationShearStress(std::fabs(TextFileParser("inputFiles/DD.txt").readScalar<double>("criticalSurfaceDislocationNucleationShearStress",true)))
        {
            
        }
        
        void injectDislocaitons()
        {
            //insertRandomStraightDislocation();
            surfaceDislocationNucleation();
        }

        /**********************************************************************/
        void surfaceDislocationNucleation()
        {

                switch (surfaceDislocationNucleationModel)
                {
                        
                    case 1:
                    {
                        std::cout<<"Nucleating surface dislocations: "<<std::flush;
                        const auto t0= std::chrono::system_clock::now();

                        const double R(50.0);
                        const size_t N(10);

                        std::set<double> mrss;
                        std::deque<SurfaceNucleationTupleType> simplexDeque; // face,slipSystemID,boundaryFace
                        
                        for (const auto& face : DN.mesh.template observer<dim-1>())
                        {// loop over mesh triangles
                            if(face.second->isBoundarySimplex())
                            {// triangle is on external boundary
                                const auto parents(face.second->parents());
                                assert(parents.size()==1);
                                const auto simplex(parents.begin()->second);
                                
                                const VectorDim pos(face.second->center()); // center of element in cartesian coordinates
                                MatrixDim sigma(DN.stress(pos)); // dislocation stress
                                if(DN.bvpSolver)
                                {
                                    const auto ele(DN.bvpSolver->finiteElement().elements().at(simplex->xID));
                                    sigma+=DN.bvpSolver->stress(ele,simplex->pos2bary(pos)); // FEM stress
                                }
                                if(DN.externalLoadController)
                                {
                                    sigma+=DN.externalLoadController->stress(pos); // external stress
                                }
                                
                                const auto& grain(DN.poly.grain(simplex->region->regionID));
                                std::map<double,const SlipSystem* const> rssMap;
                                for (const auto& slipSystem : grain.slipSystems())
                                {
                                    Eigen::Matrix<double,dim,1> n(slipSystem->unitNormal);
                                    Eigen::Matrix<double,dim,1> b(slipSystem->s.cartesian().normalized());
                                    if(n.cross(face.second->outNormal()).norm()>FLT_EPSILON)
                                    {// non-parallel slip system and face
                                        const double rss((sigma*n).dot(b)); // NOTE: a NEGATIVE rss tends to expand the loop on that slip system
                                        if(rss<-criticalSurfaceDislocationNucleationShearStress)
                                        {
                                            rssMap.emplace(rss,slipSystem.get());
//                                            simplexDeque.emplace_back(face.second,slipSystem.get(),R,N);
                                        }
                                    }
                                }
                                if(rssMap.size())
                                {// nucleate on slip system with most negative rss at that point
                                    simplexDeque.emplace_back(face.second,rssMap.begin()->second,R,N);
                                }
                            }
                        }
                        
                        size_t nucleated(0);
                        for(const auto& tup : simplexDeque)
                        {
                            nucleated+=nucleateSurfaceLoop(tup);
                        }
                        
                        model::cout<<nucleated<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
                        break;
                    }
                        
                    case 2:
                    {
                        //for (const auto& ele : DN.bvpSolver->finiteElement().elements())
                        std::cout<<"Nucleating sub-surface dislocations: "<<std::flush;
                        std::set<double> mrss;
                        for (const auto& simplexPair : DN.mesh.simplices())
                        {
                            const auto& simplexKey(simplexPair.first);
                            const auto& simplex(simplexPair.second);
                            if(simplex.isBoundarySimplex())
                            {
                                const auto& grain(DN.poly.grain(simplex.region->regionID));
                                const Eigen::Matrix<double,dim+1,1> bary(0.25*Eigen::Matrix<double,dim+1,1>::Ones()); // center of element in bary coordinates
                                const VectorDim pos(simplex.bary2pos(bary)); // center of element in cartesian coordinates
                                
                                // Stress at center of element
                                MatrixDim sigma(DN.stress(pos)); // dislocation stress
                                if(DN.bvpSolver)
                                {
                                    const auto ele(DN.bvpSolver->finiteElement().elements().at(simplexKey));
                                    sigma+=DN.bvpSolver->stress(ele,bary); // FEM stress
                                }
                                if(DN.externalLoadController)
                                {
                                    sigma+=DN.externalLoadController->stress(pos); // external stress
                                }
                                
                                // Compute RSS for each slip system and sort them
                                std::map<double,int> tauMap; // shear stress of each slip system, sorted from lowest to highest
                                
                                size_t k=0;
                                for (const auto& slipSystem : grain.slipSystems())
                                {
                                    Eigen::Matrix<double,dim,1> n(slipSystem->unitNormal);
                                    Eigen::Matrix<double,dim,1> b(slipSystem->s.cartesian().normalized());
//                                    const double rss(-(sigma*n).dot(b)); // tauRSS=-(sigma*n).dot(b), a positive tauRSS tends to expand the loop on that slip system
                                    const double rss((sigma*n).dot(b)); // NOTE: a NEGATIVE rss tends to expand the loop on that slip system
                                    tauMap.emplace(rss,k);
                                    k++;
//                                    mrss.insert(rss);
//                                    if(rss>criticalSurfaceDislocationNucleationShearStress)
//                                    {// nucleate on slip system of highest tauRSS
//                                        //                const double R(0.5*cbrt(8.4853*ele.second.simplex.vol0));
//                                        const double R(5.0);
//                                        //const auto& slipSystem(grain.slipSystems()[tauMap.rbegin()->second]);
//                                        LatticeVector<3> L(grain.snapToLattice(pos));
//                                        nucleateLoop(simplex,L,slipSystem,R);
//                                    }
                                }
                                
                                const auto& minRSS(*tauMap.begin());

                                mrss.insert(minRSS.first);
                                
                                if(minRSS.first<-criticalSurfaceDislocationNucleationShearStress)
                                {// nucleate on slip system of most negative  rss
                                    //                const double R(0.5*cbrt(8.4853*ele.second.simplex.vol0));
                                    const double R(5.0);
                                    const auto& slipSystem(grain.slipSystems()[minRSS.second]);
                                    LatticeVector<3> L(grain.snapToLattice(pos));
                                    nucleateElementLoop(simplex,L.cartesian(),slipSystem,R,slipSystem->s.cartesian().normalized(),4);
                                }
                                
                                
                            }
                        }
                        
                        std::cout<<" min rss="<<*mrss.begin()<<std::endl;

                        break;
                    }
                      
                        
                    default:
                    {
                        break;
                    }
                }
            

        }
        
        /**********************************************************************/
        bool nucleateSurfaceLoop(const SurfaceNucleationTupleType& tup)
        {
            const auto& meshFace(std::get<0>(tup));
            const auto& slipSystem(std::get<1>(tup));
            const auto& R(std::get<2>(tup));
            const auto& N(std::get<3>(tup));

            const long int glidePlaneIndex(slipSystem->n.closestPlaneIndexOfPoint(meshFace->center()));
            GlidePlaneKey<dim> glidePlaneKey(glidePlaneIndex,slipSystem->n);
            const auto glidePlane(DN.glidePlaneFactory.get(glidePlaneKey));
        
            PlanePlaneIntersection<dim> ppi(*glidePlane,Plane<3>(meshFace->center(),meshFace->outNormal()));
            const auto bndPair(ppi.snapToIntersection(meshFace->center()));
            
            if(bndPair.first)
            {
                const double dTheta=2.0*M_PI/N;
                std::vector<Eigen::Vector2d> circlePts2D;
                for(size_t i=0;i<N;++i)
                {// collect 2d bounding box (must be right-handed)
                    const MatrixDim rot(Eigen::AngleAxisd(i*dTheta,glidePlane->unitNormal));
                    circlePts2D.emplace_back(glidePlane->localPosition(bndPair.second+rot*(R*ppi.d)));
                }
                
                std::vector<Eigen::Vector2d> box2d;
                for(const auto& segment : glidePlane->meshIntersections)
                {// collect 2d bounding box (must be right-handed)
                    box2d.emplace_back(glidePlane->localPosition(segment->P0));
                }
                
                std::vector<Eigen::Vector2d> clippedCircle(SutherlandHodgman::clip(circlePts2D,box2d));
                
                std::vector<std::shared_ptr<NodeType>> loopNodes;
                for(const auto& pos : clippedCircle)
                {
                    loopNodes.emplace_back(new NodeType(&DN,glidePlane->globalPosition(pos),VectorDim::Zero(),1.0));
                }
                if(glidePlane->unitNormal.dot(slipSystem->unitNormal)>0.0)
                {
                    DN.insertLoop(loopNodes,slipSystem->s.cartesian(),glidePlane);
                }
                else
                {
                    DN.insertLoop(loopNodes,-slipSystem->s.cartesian(),glidePlane);
                }
                return true;
            }
            return false;
        }
        
        /**********************************************************************/
        void nucleateElementLoop(const Simplex<dim,dim>& simplex,
                          const VectorDim& P0,
                          //                  const LatticePlaneBase& n,
                          //                  const LatticeVector<3>& b,
                          const std::shared_ptr<SlipSystem>& slipSystem,
                          const double& R,
                          const VectorDim& d,
                                 const size_t N)
        {

            // Construct node posisions in loop
//            const VectorDim P0(L.cartesian());
//            const size_t N=8;
            const double dTheta=2.0*M_PI/N;
            const VectorDim v=R*d;
            std::vector<VectorDim> nodePos;
            for(size_t i=0;i<N;++i)
            {
                const double theta=i*dTheta;
                const MatrixDim rot(Eigen::AngleAxisd(theta,slipSystem->unitNormal));
                
                const VectorDim pos(P0+rot*v);
                if(DN.mesh.searchRegion(simplex.region->regionID,pos).first)
                {
                    nodePos.emplace_back(P0+rot*v);
                }
            }
            
            // If all nodes are inside grain create loop
//            if (allPointsInGrain(nodePos,DN.poly.grain(simplex.region->regionID).grainID)) // all points inside
                if (nodePos.size()>=3) // all points inside
            {
                std::cout<<"NUCLEATING LOOP"<<std::endl;
                
                std::vector<std::shared_ptr<NodeType>> loopNodes;
                for(const auto& pos : nodePos)
                {
                    loopNodes.emplace_back(new NodeType(&DN,pos,VectorDim::Zero(),1.0));
                }
                GlidePlaneKey<dim> loopPlaneKey(nodePos[0],slipSystem->n);
                DN.insertLoop(loopNodes,slipSystem->s.cartesian(),DN.glidePlaneFactory.get(loopPlaneKey));
            }
            else
            {
                std::cout<<"nodePos.size()="<<nodePos.size()<<std::endl;
            }
        }
        
//        bool allPointsInGrain(const std::vector<VectorDim>& points,const int& grainID) const
//        {
//            bool temp=true;
//            for(const auto& point : points)
//            {
//                temp*=DN.mesh.searchRegion(grainID,point).first;
//                if(!temp)
//                {
//                    break;
//                }
//            }
//            return temp;
//        }
        
        /**********************************************************************/
        void insertRandomStraightDislocation()
        {
            const std::pair<LatticeVector<dim>,int> rp=DN.poly.randomLatticePointInMesh();
            const int& grainID=rp.second;   // random grain ID
            const LatticeVector<dim>& L0=rp.first; // random lattice position in the grain
            const VectorDim P0(L0.cartesian());   // cartesian position of L0
            std::mt19937 generator(std::chrono::system_clock::now().time_since_epoch().count());
            std::uniform_int_distribution<> distribution(0,DN.poly.grain(grainID).slipSystems().size()-1);
            const int rSS=distribution(generator); // a random SlipSystem ID
            const auto& slipSystem=DN.poly.grain(grainID).slipSystems()[rSS];
            const VectorDim b=slipSystem->s.cartesian();    // Burgers vector
            const VectorDim n=slipSystem->n.cartesian().normalized(); // slip plane normal
            
            std::uniform_real_distribution<> dis(0.0, 2.0*M_PI);
            const double theta=dis(generator); // random angle of the dislocation line in the plane from screw orientation.
            const VectorDim d=Eigen::AngleAxisd(theta,n)*b.normalized(); // a random direction in the glide plane
            
            std::vector<VectorDim> nodePos(this->straightLineBoundaryClosure(P0,d,n,grainID,DN.mesh));
            std::vector<std::shared_ptr<NodeType>> loopNodes;
            for (const auto& pos : nodePos)
            {
                loopNodes.emplace_back(new NodeType(&DN,pos,VectorDim::Zero(),1.0));
            }
            GlidePlaneKey<dim> loopPlaneKey(nodePos[0],slipSystem->n);
            DN.insertLoop(loopNodes,b,DN.glidePlaneFactory.get(loopPlaneKey));
        }
                
    };
    
}
#endif
