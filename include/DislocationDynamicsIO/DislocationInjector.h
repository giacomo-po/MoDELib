/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationInjector_H_
#define model_DislocationInjector_H_




namespace model
{
    template <int dim>
    class DislocationInjector
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
            
            return nodePos;
        }
        
        
        /**********************************************************************/
        template <typename DislocationNetworkType>
        void insertRandomStraightDislocation(DislocationNetworkType& DN)
        {
            typedef typename DislocationNetworkType::NodeType NodeType;
            const std::pair<LatticeVector<dim>,int> rp=DN.poly.randomLatticePointInMesh();
            const int& grainID=rp.second;   // random grain ID
            const LatticeVector<dim>& L0=rp.first; // random lattice position in the grain
            const VectorDim P0(L0.cartesian());   // cartesian position of L0
            std::mt19937 generator(std::chrono::system_clock::now().time_since_epoch().count());
            std::uniform_int_distribution<> distribution(0,DN.poly.grain(grainID).slipSystems().size()-1);
            const int rSS=distribution(generator); // a random SlipSystem ID
            const SlipSystem& slipSystem=DN.poly.grain(grainID).slipSystems()[rSS];
            const VectorDim b=slipSystem.s.cartesian();    // Burgers vector
            const VectorDim n=slipSystem.n.cartesian().normalized(); // slip plane normal
            
            std::uniform_real_distribution<> dis(0.0, 2.0*std::numbers::pi);
            const double theta=dis(generator); // random angle of the dislocation line in the plane from screw orientation.
            const VectorDim d=Eigen::AngleAxisd(theta,n)*b.normalized(); // a random direction in the glide plane
            
            std::vector<VectorDim> nodePos(straightLineBoundaryClosure(P0,d,n,grainID,DN.mesh));
            std::vector<std::shared_ptr<NodeType>> loopNodes;
            for (const auto& pos : nodePos)
            {
                loopNodes.emplace_back(new NodeType(&DN,pos,VectorDim::Zero(),1.0));
            }
            DN.insertLoop(loopNodes,b,n,P0,grainID);
//            DN.clearDanglingNodes();
        }
        
    };
    
}
#endif
