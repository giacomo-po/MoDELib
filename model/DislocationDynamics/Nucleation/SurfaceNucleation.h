#ifndef model_SurfaceNucleation_H_
#define model_SurfaceNucleation_H_

#include <map>
#include <Eigen/Dense>
#include <model/LatticeMath/LatticeMath.h>


namespace model
{
    /**************************************************************************/
    /**************************************************************************/
    template<int SurfaceNucleationModel>
    struct SurfaceNucleation
    {
        template <typename DislocationNetworkType>
        static size_t nucleateDislocations(DislocationNetworkType& )
        {
            return 0;
        }
        
    };
    
    /**************************************************************************/
    /**************************************************************************/
    /*! Class template that nucleates shear loops underneath the surface of the mesh
     * based on a max rss criterion
     */
    template<>
    struct SurfaceNucleation<1>
    {
        
        /**********************************************************************/
        template <typename DislocationNetworkType>
        static size_t nucleateLoop(DislocationNetworkType& DN,
                                   const Simplex<DislocationNetworkType::dim,DislocationNetworkType::dim>& simplex,
                                   const LatticeVector<3>& L,
                                   const LatticePlaneBase& n,
                                   const LatticeVector<3>& b,
                                   const double R)
        {
            size_t numberNucleated=0;
            
            constexpr int dim=DislocationNetworkType::dim;
            
            LatticePlane plane(L,n);
            
            const Eigen::Matrix<double,dim,1> P0(L.cartesian());
            const size_t N=8;
            const double dTheta=2.0*M_PI/N;
            const Eigen::Matrix<double,dim,1> v=R*b.cartesian().normalized();
            
            std::deque<LatticeVector<3>> deq;
            
            for(size_t i=0;i<N;++i)
            {
                const double theta=i*dTheta;
                const Eigen::Matrix<double,dim,dim> rot(Eigen::AngleAxisd(theta,n.cartesian()));
                Eigen::Matrix<double,dim,1> P(P0+rot*v);
                P=plane.snapToLattice(P);
                if(DN.shared.mesh.searchWithGuess(P,&simplex).first)
                {
                    deq.emplace_back(P);
                }
            }
            
            if (deq.size()==N) // all points inside
            {
                std::cout<<"NUCLEATING LOOP"<<std::endl;
                
                std::deque<size_t> vIDs;
                for(size_t i=0;i<N;++i)
                {
                    vIDs.push_back(DN.insertVertex(deq[i]).first->first);
                }
                
                for(size_t i=0;i<N-1;++i)
                {
                    DN.connect(vIDs[i],vIDs[i+1],b);
                }
                DN.connect(vIDs[N-1],vIDs[0],b);
                
                ++numberNucleated;
            }
            return numberNucleated;
        }
        

        
        /**********************************************************************/
        template <typename DislocationNetworkType>
        static size_t nucleateDislocations(DislocationNetworkType& DN)
        {
            
            constexpr int dim=DislocationNetworkType::dim;
            
            size_t numberNucleated=0;
            const double tauCR=0.07; // critica nucleation stress
            const double R(5.0);     // nucleation radius
            
            
            for (const auto& ele : DN.shared.bvpSolver.finiteElement().elements())
            {
                if(ele.second.simplex.isBoundarySimplex())
                {
                    const Eigen::Matrix<double,dim+1,1> bary(0.25*Eigen::Matrix<double,dim+1,1>::Ones()); // center of the element in barycentric coordinates
                    const Eigen::Matrix<double,dim,1> pos(ele.second.simplex.bary2pos(bary));
                    const Eigen::Matrix<double,dim,dim> sigma(DN.shared.bvpSolver.stress(ele.second,bary)+DN.stress(pos));
                    
                    std::map<double,int> tauMap;

                    size_t k=0;
                    for (const auto& slipSystem : CrystalOrientation<dim>::slipSystems())
                    {
                        Eigen::Matrix<double,dim,1> n(slipSystem.n.cartesian().normalized());
                        Eigen::Matrix<double,dim,1> b(slipSystem.s.cartesian().normalized());
                        tauMap.emplace(-(sigma*n).dot(b),k); // tauRSS=-(sigma*n).dot(b) is the stress expanding loop of current slip system
                        k++;
                    }

                    if(tauMap.rbegin()->first>tauCR)
                    {
                        LatticeVector<3> L(LatticeBase<dim>::snapToLattice(pos));
                        const auto& slipSystem(CrystalOrientation<dim>::slipSystems()[tauMap.rbegin()->second]);
                        numberNucleated+=nucleateLoop(DN,ele.second.simplex,L,slipSystem.n,slipSystem.s,R);
                    }
                    
//                    LatticeVector<3> L(LatticeBase<dim>::snapToLattice(pos));
//                    for (const auto& slipSystem : CrystalOrientation<dim>::slipSystems())
//                    {
//                        Eigen::Matrix<double,dim,1> n(slipSystem.n.cartesian().normalized());
//                        Eigen::Matrix<double,dim,1> b(slipSystem.s.cartesian().normalized());
//                        if(-(sigma*n).dot(b)>tauCR)
//                        {
//                            numberNucleated+=nucleateLoop(DN,ele.second.simplex,L,slipSystem.n,slipSystem.s,R);
//                            
//                        }
//                    }
                    
                }
            }
            return numberNucleated;
        }
        
        
    };
    
    /**************************************************************************/
    /**************************************************************************/
    template<>
    struct SurfaceNucleation<2>
    {/*! Class template that nucleates half shear loops on the surface of the mesh
      * based on a max rss criterion
      */
        
        /**********************************************************************/
        template <typename DislocationNetworkType>
        static size_t nucleateHalfLoop(DislocationNetworkType& DN,
                                   const Simplex<DislocationNetworkType::dim,DislocationNetworkType::dim>& simplex,
                                   const int& boundaryFaceID,
//                                       const Simplex<3,2>& face,
                                   const LatticeVector<3>& L,
                                   const LatticePlaneBase& n,
                                   const LatticeVector<3>& b,
                                   const double R)
        {
            size_t numberNucleated=0;
            
            constexpr int dim=DislocationNetworkType::dim;
            
            LatticePlane plane(L,n);
            
            const Eigen::Matrix<double,dim,1> P0(L.cartesian());
            const size_t N=5;
            const double dTheta=M_PI/(N-1);
            Eigen::Matrix<double,dim,1> v=n.cartesian().cross(simplex.child(boundaryFaceID).outNormal());
            const double vNorm=v.norm();
            if(vNorm>FLT_EPSILON)
            {
                v*=R/vNorm;
                std::deque<LatticeVector<3>> deq;
                
                double sgn=1.0;
                const Eigen::Matrix<double,dim,dim> rotTest(Eigen::AngleAxisd(0.5*M_PI,n.cartesian())); // test if mid point along circle is inside or outside
                if(!DN.shared.mesh.searchWithGuess(P0+rotTest*v,&simplex).first)
                {
                    sgn=-1.0;
                }

                
                
                for(size_t i=0;i<N;++i)
                {
                    const double theta=i*dTheta;
                    const Eigen::Matrix<double,dim,dim> rot(Eigen::AngleAxisd(sgn*theta,n.cartesian()));
                    Eigen::Matrix<double,dim,1> P(P0+rot*v);
                    P=plane.snapToLattice(P);
                    if(DN.shared.mesh.searchWithGuess(P,&simplex).first)
                    {
                        deq.emplace_back(P);
                    }
                    else
                    {
                        std::cout<<"NOT CREATING NODE OUTSIDE"<<std::endl;
                    }
                }
                
                if (deq.size()==N) // all points inside
                {
                    std::cout<<"NUCLEATING HALF LOOP"<<std::endl;
                    
                    std::deque<size_t> vIDs;
                    for(size_t i=0;i<N;++i)
                    {
                        vIDs.push_back(DN.insertVertex(deq[i]).first->first);
                    }
                    
                    for(size_t i=0;i<N-1;++i)
                    {
                        DN.connect(vIDs[i],vIDs[i+1],b);
                    }
                    DN.connect(vIDs[N-1],vIDs[0],b);
                    
                    ++numberNucleated;
                }
            }
            

            return numberNucleated;
        }
        
        
        
        /**********************************************************************/
        template <typename DislocationNetworkType>
        static size_t nucleateDislocations(DislocationNetworkType& DN)
        {
            
            constexpr int dim=DislocationNetworkType::dim;
            
            size_t numberNucleated=0;
            const double tauCR=0.07; // critica nucleation stress
            const double R(10.0);     // nucleation radius
            
            
            for (const auto& ele : DN.shared.bvpSolver.finiteElement().elements())
            {
                const Simplex<dim,dim>& simplex(ele.second.simplex);
                
                if(simplex.isBoundarySimplex())
                {
                    
                    
                    const std::vector<int> boundaryFaces=simplex.boundaryFaces();
                    for (unsigned int f=0;f<boundaryFaces.size();++f) // loop over boundary faces of the current element
                    {
                        const int boundaryFaceID(boundaryFaces[f]);
                        Eigen::Matrix<double,dim+1,1> bary(1.0/3.0*Eigen::Matrix<double,dim+1,1>::Ones());
                        bary(boundaryFaceID)=0.0; // center of the face in barycentric coordinates
                        const Eigen::Matrix<double,dim,1> pos(simplex.bary2pos(bary));
                        const Eigen::Matrix<double,dim,dim> sigma(DN.shared.bvpSolver.stress(ele.second,bary)+DN.stress(pos));

                        std::map<double,int> tauMap;
                        
                        size_t k=0;
                        for (const auto& slipSystem : CrystalOrientation<dim>::slipSystems())
                        {
                            Eigen::Matrix<double,dim,1> n(slipSystem.n.cartesian().normalized());
                            Eigen::Matrix<double,dim,1> b(slipSystem.s.cartesian().normalized());
                            if(n.cross(simplex.child(boundaryFaceID).outNormal().normalized()).norm()>FLT_EPSILON)
                            {
                                tauMap.emplace(-(sigma*n).dot(b),k); // tauRSS=-(sigma*n).dot(b) is the stress expanding loop of current slip system
                            }
                            k++;
                        }
                        
                        if(tauMap.rbegin()->first>tauCR)
                        {
                            LatticeVector<3> L(LatticeBase<dim>::snapToLattice(pos-simplex.child(boundaryFaceID).outNormal().normalized()*2.0)); //shift center of face slightly inside mesh
                            const auto& slipSystem(CrystalOrientation<dim>::slipSystems()[tauMap.rbegin()->second]);
                            numberNucleated+=nucleateHalfLoop(DN,simplex,boundaryFaceID,L,slipSystem.n,slipSystem.s,R);
//                            numberNucleated+=nucleateHalfLoop(DN,simplex,simplex.child(boundaryFaceID),L,slipSystem.n,slipSystem.s,R);

                        }

                    
                    }

                    
                    
//                    const Eigen::Matrix<double,dim+1,1> bary(0.25*Eigen::Matrix<double,dim+1,1>::Ones()); // center of the element in barycentric coordinates
//                    const Eigen::Matrix<double,dim,1> pos(ele.second.simplex.bary2pos(bary));
//                    const Eigen::Matrix<double,dim,dim> sigma(DN.shared.bvpSolver.stress(ele.second,bary)+DN.stress(pos));
                    
                    
                    //                    LatticeVector<3> L(LatticeBase<dim>::snapToLattice(pos));
                    //                    for (const auto& slipSystem : CrystalOrientation<dim>::slipSystems())
                    //                    {
                    //                        Eigen::Matrix<double,dim,1> n(slipSystem.n.cartesian().normalized());
                    //                        Eigen::Matrix<double,dim,1> b(slipSystem.s.cartesian().normalized());
                    //                        if(-(sigma*n).dot(b)>tauCR)
                    //                        {
                    //                            numberNucleated+=nucleateLoop(DN,ele.second.simplex,L,slipSystem.n,slipSystem.s,R);
                    //
                    //                        }
                    //                    }
                    
                }
            }
            return numberNucleated;
        }
        
        
    };
    
    
    //    /**************************************************************************/
    //    /**************************************************************************/
    //    template<>
    //    struct SurfaceNucleation<2>
    //    {
    //    };
}


#endif


