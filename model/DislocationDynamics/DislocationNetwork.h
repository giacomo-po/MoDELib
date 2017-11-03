/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po       <gpo@ucla.edu>.
 * Copyright (C) 2011 by Benjamin Ramirez <ramirezbrf@gmail.com>
 * Copyright (C) 2011 by Mamdouh Mohamed  <msm07d@fsu.edu>
 * Copyright (C) 2011 by Tamer Crsoby     <tamercrosby@gmail.com>
 * Copyright (C) 2011 by Can Erel         <canerel55@gmail.com>
 * Copyright (C) 2011 by Yinan Cui        <cuiyinan@ucla.edu>
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DISLOCATIONNETWORK_H_
#define model_DISLOCATIONNETWORK_H_

#ifdef _MODEL_MPI_
#define _MODEL_DD_MPI_
#endif


#ifdef _OPENMP
#include <omp.h>
#define EIGEN_DONT_PARALLELIZE // disable Eigen Internal openmp Parallelization
#endif

#include <chrono>
#include <Eigen/Dense>

#include <model/LoopNetwork/LoopNetwork.h>
#include <model/Utilities/TerminalColors.h>
#include <model/IO/EigenDataReader.h>
#include <model/DislocationDynamics/DislocationNetworkTraits.h>
#include <model/DislocationDynamics/DislocationNetworkComponent.h>
#include <model/DislocationDynamics/DislocationNode.h>
#include <model/DislocationDynamics/DislocationSegment.h>
#include <model/DislocationDynamics/DislocationLoop.h>
#include <model/DislocationDynamics/GlidePlanes/GlidePlaneObserver.h>
#include <model/DislocationDynamics/DislocationNetworkRemesh.h>
#include <model/DislocationDynamics/Junctions/DislocationJunctionFormation.h>
//#include <model/DislocationDynamics/CrossSlip/DislocationCrossSlip.h>
#include <model/DislocationDynamics/Materials/Material.h>
#include <model/DislocationDynamics/IO/DislocationNetworkIO.h>
#include <model/DislocationDynamics/NearestNeighbor/DislocationParticle.h>
#include <model/DislocationDynamics/NearestNeighbor/DislocationStress.h>
#include <model/ParticleInteraction/ParticleSystem.h>
#include <model/MPI/MPIcout.h> // defines mode::cout
#include <model/ParticleInteraction/SingleFieldPoint.h>
#include <model/DislocationDynamics/DDtimeIntegrator.h>
#include <model/Threads/EqualIteratorRange.h>
#include <model/DislocationDynamics/BoundingLineSegments.h>
#include <model/DislocationDynamics/Polycrystals/GrainBoundaryTransmission.h>
//#include <model/DislocationDynamics/Polycrystals/GrainBoundaryDissociation.h>
#include <model/DislocationDynamics/BVP/BVPsolver.h>
#include <model/DislocationDynamics/Polycrystals/Polycrystal.h>
#include <model/DislocationDynamics/DislocationNodeContraction.h>


namespace model
{
    
    template <int _dim, short unsigned int corder, typename InterpolationType>
    class DislocationNetwork : public LoopNetwork<DislocationNetwork<_dim,corder,InterpolationType> >,
    /* base                 */ public GlidePlaneObserver<DislocationNetwork<_dim,corder,InterpolationType>>,
    /* base                 */ public ParticleSystem<DislocationParticle<_dim> >
    {
        
    public:
        
        static constexpr int dim=_dim; // make dim available outside class
        typedef DislocationNetwork<dim,corder,InterpolationType> DislocationNetworkType;
        typedef LoopNetwork<DislocationNetworkType> LoopNetworkType;
        typedef DislocationNode<dim,corder,InterpolationType> NodeType;
        typedef DislocationSegment<dim,corder,InterpolationType> LinkType;
        typedef DislocationNetworkComponent<NodeType,LinkType> DislocationNetworkComponentType;
        typedef NetworkComponent<NodeType,LinkType> NetworkComponentType;
        typedef Eigen::Matrix<double,dim,dim>	MatrixDimD;
        typedef Eigen::Matrix<double,dim,1>		VectorDim;
        typedef typename TypeTraits<DislocationNetworkType>::LoopType LoopType;
        typedef GlidePlaneObserver<DislocationNetworkType> GlidePlaneObserverType;
        typedef DislocationParticle<_dim> DislocationParticleType;
        typedef typename DislocationParticleType::StressField StressField;
        typedef typename DislocationParticleType::DisplacementField DisplacementField;
        typedef ParticleSystem<DislocationParticleType> ParticleSystemType;
        typedef typename ParticleSystemType::SpatialCellType SpatialCellType;
        typedef SpatialCellObserver<DislocationParticleType,_dim> SpatialCellObserverType;
        typedef BVPsolver<dim,2> BvpSolverType;
        typedef typename BvpSolverType::FiniteElementType FiniteElementType;
        typedef typename FiniteElementType::ElementType ElementType;
        typedef typename LoopNetworkType::IsNodeType IsNodeType;
        typedef DislocationNetworkIO<DislocationNetworkType> DislocationNetworkIOType;
        //        enum {NdofXnode=NodeType::NdofXnode};
        
        int timeIntegrationMethod;
        bool use_junctions;
        long int runID;
        double totalTime;
        double dt;
        double vMax;
        int Nsteps;
        MatrixDimD _plasticDistortion;
        MatrixDimD _plasticDistortionRate;
        MatrixDimD externalStress;
        bool use_boundary;
        unsigned int use_bvp;
        bool use_virtualSegments;
        SimplicialMesh<dim> mesh;
        Polycrystal<DislocationNetworkType> poly;
        BvpSolverType bvpSolver;
        
#ifdef DislocationNucleationFile
#include DislocationNucleationFile
#endif
        
    private:
        
        /**********************************************************************/
        void update_BVP_Solution()
        {/*! Updates bvpSolver using the stress and displacement fields of the
          *  current DD configuration.
          */
            const int quadraturePerTriangle=37;
            if(use_bvp)
            {
                if (!(runID%use_bvp))
                {// enter the if statement if use_bvp!=0 and runID is a multiple of use_bvp
                    model::cout<<"		Updating elastic bvp... "<<std::endl;
                    bvpSolver.template assembleAndSolve<DislocationNetworkType,quadraturePerTriangle>(*this);
                }
            }
        }
        
        /**********************************************************************/
        void computeNodaVelocities()
        {
            
            switch (timeIntegrationMethod)
            {
                case 0:
                    DDtimeIntegrator<0>::computeNodaVelocities(*this);
                    break;
                    
                    //                case 1:
                    //                    dt=DDtimeIntegrator<1>::integrate(*this);
                    //                    break;
                    
                default:
                    assert(0 && "time integration method not implemented");
                    break;
            }
            
            
            for(auto& node : this->nodes())
            {
                node.second->applyVelocityFilter(vMax);
            }
            
            model::cout<<std::setprecision(3)<<std::scientific<<"		dt="<<dt<<std::endl;
        }
        
        /**********************************************************************/
        void singleStep()
        {
            //! A simulation step consists of the following:
            model::cout<<blueBoldColor<< "runID="<<runID<<" (of "<<Nsteps<<")"
            /*                    */<< ", time="<<totalTime
            /*                    */<< ": nodes="<<this->nodes().size()
            /*                    */<< ", segments="<<this->links().size()
            /*                    */<< ", loopSegments="<<this->loopLinks().size()
            /*                    */<< ", loops="<<this->loops().size()
            /*                    */<< ", components="<<this->components().size()
            /*                    */<< defaultColor<<std::endl;
            
            //! 1- Check that all nodes are balanced
            //            checkBalance();
            
            //! 2 - Update quadrature points
            updateQuadraturePoints();
            
            //! 3- Calculate BVP correction
            update_BVP_Solution();
            
#ifdef DislocationNucleationFile
            if(use_bvp && !(runID%use_bvp))
            {
                nucleateDislocations(); // needs to be called before updateQuadraturePoints()
                updateQuadraturePoints();
            }
#endif
            
            computeNodaVelocities();
            
            
            //! 4- Solve the equation of motion
            
            
            //! 5- Compute time step dt (based on max nodal velocity) and increment totalTime
            // make_dt();
            
            if(DislocationNetworkIO<DislocationNetworkType>::outputElasticEnergy)
            {
                typedef typename DislocationParticleType::ElasticEnergy ElasticEnergy;
                this->template computeNeighborField<ElasticEnergy>();
            }
            
            //! 6- Output the current configuration before changing it
            //            output(runID);
            io().output(runID);
            
            
            //! 7- Moves DislocationNodes(s) to their new configuration using stored velocity and dt
            move(dt);
            
            //! 8- Update accumulated quantities (totalTime and plasticDistortion)
            totalTime+=dt;
            computePlasticDistortionRate();
            _plasticDistortion += _plasticDistortionRate*dt;
            
            
            //! 9- Contract segments of zero-length
            //            DislocationNetworkRemesh<DislocationNetworkType>(*this).contract0chordSegments();
            
            //! 10- Cross Slip (needs upated PK force)
            //            DislocationCrossSlip<DislocationNetworkType>(*this);
            
            
            //                        GrainBoundaryTransmission<DislocationNetworkType>(*this).transmit();
            GrainBoundaryTransmission<DislocationNetworkType>(*this);
            
            //
            //            GrainBoundaryDissociation<DislocationNetworkType>(*this).dissociate();
            
            //            poly.reConnectGrainBoundarySegments(*this); // this makes stressGauss invalid, so must follw other GB operations
            
            
            //! 11- detect loops that shrink to zero and expand as inverted loops
            //            DislocationNetworkRemesh<DislocationNetworkType>(*this).loopInversion(dt);
            
            //! 12- Form Junctions
            DislocationJunctionFormation<DislocationNetworkType>(*this).formJunctions(use_junctions,DDtimeIntegrator<0>::dxMax);
            
            //            // Remesh may contract juncitons to zero lenght. Remove those juncitons:
            //            DislocationJunctionFormation<DislocationNetworkType>(*this).breakZeroLengthJunctions();
            
            //! 13- Node redistribution
            DislocationNetworkRemesh<DislocationNetworkType>(*this).remesh(runID);
            
            //! 9- Contract segments of zero-length
            //            DislocationNetworkRemesh<DislocationNetworkType>(*this).contract0chordSegments();
            
            //! 14- If BVP solver is not used, remove DislocationSegment(s) that exited the boundary
            //            removeBoundarySegments();
            
            //            removeSmallComponents(3.0*dx,4);
            
            //            make_bndNormals();
            
            //! 16 - Increment runID counter
            ++runID;     // increment the runID counter
        }
        
    public:
        
        /**********************************************************************/
        DislocationNetwork(int& argc, char* argv[]) :
        /* init list  */ timeIntegrationMethod(0),
        /* init list  */ use_junctions(false),
        /* init list  */ runID(0),
        /* init list  */ totalTime(0.0),
        /* init list  */ dt(0.0),
        /* init list  */ vMax(0.0),
        /* init list  */ Nsteps(0),
        /* init list  */ _plasticDistortion(MatrixDimD::Zero()),
        /* init list  */ _plasticDistortionRate(MatrixDimD::Zero()),
        /* init list  */ externalStress(MatrixDimD::Zero()),
        /* init list  */ use_boundary(false),
        /* init list  */ use_bvp(0),
        /* init list  */ use_virtualSegments(true),
        /* init list  */ poly(mesh),
        /* init list  */ bvpSolver(mesh)
        {
            ParticleSystemType::initMPI(argc,argv);
            io().read("./","DDinput.txt");
            
            // Initializing configuration
            move(0.0);	// initial configuration
        }
        
        
//        /**********************************************************************/
//        bool contractYoungest(std::shared_ptr<NodeType> nA,
//                              std::shared_ptr<NodeType> nB)
//        {
//            return nA->sID<nB->sID? this->contractSecond(nA,nB) : this->contractSecond(nB,nA);
//        
//        }

        /**********************************************************************/
        bool contract(std::shared_ptr<NodeType> nA,
                      std::shared_ptr<NodeType> nB)
        {
            return DislocationNodeContraction<DislocationNetworkType>(*this).contract(nA,nB);
        }
        
//        /**********************************************************************/
//        bool contract(std::shared_ptr<NodeType> nA,
//                      std::shared_ptr<NodeType> nB)
//        {
//            
//            bool success=false;
//            
//            //std::pair<bool,VectorDim> c=std::make_pair(false,VectorDim::Zero());
//            
//            std::cout<<"DislocationNetwork::contract "<<nA->sID<<" "<<nB->sID<<std::endl;
//            
//
//            
//            if(nA->isGlissile() && nB->isGlissile())
//            {// both nodes are glissile
//                if(nA->isOnBoundingBox() && nB->isOnBoundingBox())
//                {// both nodes on bounding boxes. Intersect bounding boxes'
//                    
//                    std::cout<<"contractPoint case 1aa"<<std::endl;
////                    BoundingLineSegments<dim> temp(nA->boundingBoxSegments(),nB->boundingBoxSegments());
////                    std::cout<<"temp.size()="<<temp.size()<<std::endl;
////                    if(temp.size())
////                    {// a common portion of the boundary exists
////                        std::cout<<"contractPoint case 1"<<std::endl;
//
////                                    std::cout<<"ACCOUNT FOR GRAINBONUDARY SEGMENST WHEN CCHECKING LINKS"<<std::endl;
//                        
//                        std::cout<<"nA neighbors:"<<std::endl;
//                        bool nAisRemovable=true;
//                        for(const auto& pair : nA->neighbors())
//                        {
////                            std::cout<<"neighbor "<<nA->sID<<" "<<pair.first<<" ("<<std::get<0>(pair.second)->sID<<") "<<temp.contains(0.5*(std::get<0>(pair.second)->get_P()+nB->get_P())).first<<std::endl;
//                            if(pair.first!=nB->sID && std::get<1>(pair.second)->isBoundarySegment())
//                            {// boundary segments other than A->B must remain boundary
////                                std::cout<<"boundary neighbor "<<nA->sID<<" "<<pair.first<<" ("<<std::get<0>(pair.second)->sID<<") "<<temp.contains(0.5*(std::get<0>(pair.second)->get_P()+nB->get_P())).first<<std::endl;
////                                nAisRemovable*=temp.contains(0.5*(std::get<0>(pair.second)->get_P()+nB->get_P())).first;
//                                nAisRemovable*=std::get<1>(pair.second)->boundingBoxSegments().contains(0.5*(std::get<0>(pair.second)->get_P()+nB->get_P())).first;
//
//                            }
//                            
//                            if(pair.first!=nB->sID && std::get<1>(pair.second)->isGrainBoundarySegment())
//                            {// boundary segments other than A->B must remain boundary
//                                for(const auto& gb : std::get<1>(pair.second)->grainBoundaries())
//                                {
//                                    std::cout<<"grainBoundary neighbor "<<nA->sID<<" "<<pair.first<<" ("<<std::get<0>(pair.second)->sID<<") "<<gb->glidePlanes().begin()->second.contains(nB->get_P())<<std::endl;
//                                    nAisRemovable*=gb->glidePlanes().begin()->second.contains(nB->get_P());
//                                }
//                            }
//                        }
//                        
//                        std::cout<<"nB neighbors:"<<std::endl;
//                        bool nBisRemovable=true;
//                        for(const auto& pair : nB->neighbors())
//                        {
////                            std::cout<<"neighbor "<<nB->sID<<" "<<pair.first<<" ("<<std::get<0>(pair.second)->sID<<") "<<temp.contains(0.5*(std::get<0>(pair.second)->get_P()+nA->get_P())).first<<std::endl;
//                            if(pair.first!=nA->sID && std::get<1>(pair.second)->isBoundarySegment())
//                            {// boundary segments other than A->B must remain boundary
////                                std::cout<<"boundary neighbor "<<nB->sID<<" "<<pair.first<<" ("<<std::get<0>(pair.second)->sID<<") "<<temp.contains(0.5*(std::get<0>(pair.second)->get_P()+nA->get_P())).first<<std::endl;
////                                nBisRemovable*=temp.contains(0.5*(std::get<0>(pair.second)->get_P()+nA->get_P())).first;
//                                nBisRemovable*=std::get<1>(pair.second)->boundingBoxSegments().contains(0.5*(std::get<0>(pair.second)->get_P()+nA->get_P())).first;
//                            }
//                        
//                            if(pair.first!=nA->sID && std::get<1>(pair.second)->isGrainBoundarySegment())
//                            {// boundary segments other than A->B must remain boundary
//                                for(const auto& gb : std::get<1>(pair.second)->grainBoundaries())
//                                {
//                                    std::cout<<"grainBoundary neighbor "<<nB->sID<<" "<<pair.first<<" ("<<std::get<0>(pair.second)->sID<<") "<<gb->glidePlanes().begin()->second.contains(nA->get_P())<<std::endl;
//                                    nBisRemovable*=gb->glidePlanes().begin()->second.contains(nA->get_P());
//                                }
//                            }
//                        
//                        }
//                        
//                        if(nAisRemovable && nBisRemovable)
//                        {
//                            std::cout<<"contractPoint case 1a"<<std::endl;
////                            nA->set_P(std::get<0>(temp.snap(0.5*(nA->get_P()+nB->get_P()))));
////                            return this->contractSecond(nA,nB);
//                            return contractYoungest(nA,nB);
//                        }
//                        else if(nAisRemovable && !nBisRemovable)
//                        {
//                            std::cout<<"contractPoint case 1b"<<std::endl;
//                            return this->contractSecond(nB,nA);
//                        }
//                        else if(!nAisRemovable && nBisRemovable)
//                        {
//                            std::cout<<"contractPoint case 1c"<<std::endl;
//                            return this->contractSecond(nA,nB);
//                        }
//                        else
//                        {
//                            std::cout<<"contractPoint case 1d"<<std::endl;
//                            return false;
//                        }
//                        
//                        
////                        const auto tupA(temp.snap(nA->get_P()));
////                        const auto tupB(temp.snap(nB->get_P()));
////                        
////                        const size_t& bnIDA(std::get<1>(tupA));
////                        const size_t& bnIDB(std::get<1>(tupB));
////                        
////                        if(bnIDA==bnIDB)
////                        {// IDs of the bounding box lines are the same
////                            nA->set_P(std::get<0>(temp.snap(0.5*(nA->get_P()+nB->get_P()))));
////                            return this->contractSecond(nA,nB);
////                        }
////                        else
////                        {// IDs of the bounding box lines are different
////                            SegmentSegmentDistance<dim> ssd(temp[bnIDA].first,temp[bnIDA].second,
////                                                            temp[bnIDB].first,temp[bnIDB].second);
////                            if(ssd.dMin<FLT_EPSILON)
////                            {
////                                nA->set_P(0.5*(ssd.x0+ssd.x1));
////                                return this->contractSecond(nA,nB);
////                            }
////                            else
////                            {
////                                return false;
////                            }
////                        }
//                        
////                        const auto Av(temp.snapToVertex(nA->get_P()));
////                        const auto Bv(temp.snapToVertex(nB->get_P()));
////                        
////                        if(Av.first<FLT_EPSILON && Bv.first<FLT_EPSILON)
////                        {// both A and B are on vertices
////                            if((Av.second-Bv.second).squaredNorm()<FLT_EPSILON)
////                            {
////                                std::cout<<"contractPoint case 1aA"<<std::endl;
////                                return this->contractSecond(nA,nB);
////                            }
////                            else
////                            {
////                                std::cout<<"contractPoint case 1aB"<<std::endl;
////                                return false;
////                            }
////                        }
////                        else if(Av.first<FLT_EPSILON && !(Bv.first<FLT_EPSILON))
////                        {// A is on a vertex, leave A
////                            std::cout<<"contractPoint case 1aC"<<std::endl;
////                            return this->contractSecond(nA,nB);
////                        }
////                        else if(!(Av.first<FLT_EPSILON) && Bv.first<FLT_EPSILON)
////                        {// B is on a vertex, leave B
////                            std::cout<<"contractPoint case 1aD"<<std::endl;
////                            return this->contractSecond(nB,nA);
////                        }
////                        else
////                        {// neither A nor B are on vertices
////                            std::cout<<"contractPoint case 1aE"<<std::endl;
////                            nA->set_P(std::get<0>(temp.snap(0.5*(nA->get_P()+nB->get_P()))));
////                            return this->contractSecond(nA,nB);
////                        }
////                    }
////                    else
////                    {
////                        std::cout<<"contractPoint case 2"<<std::endl;
////                        return false;
////                    }
//                }
//                else if(nA->isOnBoundingBox() && !nB->isOnBoundingBox())
//                {// a on box, b is not
//                    if(nB->glidePlaneIntersections().size())
//                    {// other is confined internally by two or more planes
//                        std::cout<<"contractPoint case 3aa"<<std::endl;
////                        BoundingLineSegments<dim> temp(nA->boundingBoxSegments(),nB->glidePlaneIntersections()); // NO THIS IS WRONG. USE nB->boundingBoxSegments()
//                        BoundingLineSegments<dim> temp(nA->boundingBoxSegments(),nB->boundingBoxSegments());
//                        if(temp.size())
//                        {
////                            if(temp.contains())
//                            std::cout<<"contractPoint case 3"<<std::endl;
////                            const VectorDim X(std::get<0>(temp.snap(0.5*(nA->get_P()+nB->get_P())))); // point at which we will contract
//                            const VectorDim X(std::get<0>(temp.snap(nA->get_P()))); // point at which we will contract
//
//                            //const size_t X_ID(std::get<1>(nA->boundingBoxSegments().snap(X)));              // bonuding box line ID containing X
//                            const size_t A_ID(std::get<1>(nA->boundingBoxSegments().snap(nA->get_P())));    // bonuding box line ID containing A
//                            if(nA->boundingBoxSegments().contains(X,A_ID))
//                            {
//                                nA->set_P(X);
//                                return this->contractSecond(nA,nB);
//                            }
//                            else
//                            {
//                                return false;
//                            }
//                        }
//                        else
//                        {
//                            std::cout<<"contractPoint case 4"<<std::endl;
//                            return false;
//                        }
//                    }
//                    else
//                    {// other is confined by only one glide plane
//                        std::cout<<"contractPoint case 5aa"<<std::endl;
//                        BoundingLineSegments<dim> temp=nA->boundingBoxSegments();
//                        temp.updateWithGlidePlane(nB->glidePlane(0));
//                        std::cout<<"temp.size="<<temp.size()<<std::endl;
//                        if(temp.size())
//                        {
//                            std::cout<<"contractPoint case 5"<<std::endl;
//                            nA->set_P(std::get<0>(temp.snap(0.5*(nA->get_P()+nB->get_P()))));
//                            return this->contractSecond(nA,nB);
//                        }
//                        else
//                        {
//                            std::cout<<"contractPoint case 6"<<std::endl;
//                            return false;
//                        }
//                    }
//                }
//                else if(!nA->isOnBoundingBox() && nB->isOnBoundingBox())
//                {// b on box, a is not
//                    std::cout<<"contractPoint case 7"<<std::endl;
//                    return contract(nB,nA);
//                }
//                else
//                {// neither a nor b on bounding box
//                    if(nA->glidePlaneIntersections().size() && nB->glidePlaneIntersections().size())
//                    {// both nodes confined by more then one plane
//                        std::cout<<"contractPoint case 8aa"<<std::endl;
//                        BoundingLineSegments<dim> temp(nA->glidePlaneIntersections(),nB->glidePlaneIntersections());
//                        if(temp.size())
//                        {
//                            std::cout<<"contractPoint case 8"<<std::endl;
//                            nA->set_P(std::get<0>(temp.snap(0.5*(nA->get_P()+nB->get_P()))));
//                            return this->contractSecond(nA,nB);
//                        }
//                        else
//                        {
//                            std::cout<<"contractPoint case 9"<<std::endl;
//                            return false;
//                        }
//                    }
//                    else if(nA->glidePlaneIntersections().size() && !nB->glidePlaneIntersections().size())
//                    {// nA confined by more then one plane, nB confined by only one plane
//                        std::cout<<"contractPoint case 10aa"<<std::endl;
//                        
//                        PlaneLineIntersection<dim> pli(nB->glidePlane(0).P.cartesian(),
//                                                       nB->glidePlane(0).n.cartesian(),
//                                                       nA->glidePlaneIntersections()[0].first, // origin of line
//                                                       nA->glidePlaneIntersections()[0].second-nA->glidePlaneIntersections()[0].first // line direction
//                                                       );
//                        
//                        if(pli.type==PlaneLineIntersection<dim>::COINCIDENT)
//                        {// nothing to do, _glidePlaneIntersections remains unchanged
//                            std::cout<<"contractPoint case 10aa1"<<std::endl;
//                            nA->set_P(nA->snapToGlidePlaneIntersection(0.5*(nA->get_P()+nB->get_P())));
//                            return this->contractSecond(nA,nB);
//                        }
//                        else if(pli.type==PlaneLineIntersection<dim>::INCIDENT)
//                        {// _glidePlaneIntersections becomes a singular point
//                            std::cout<<"contractPoint case 10aa2"<<std::endl;
//                            nA->set_P(pli.P);
//                            return this->contractSecond(nA,nB);
//                        }
//                        else
//                        {// parallel planes, cannot contract
//                            std::cout<<"contractPoint case 10aa3"<<std::endl;
//                            return false;
//                        }
//                        
////                        BoundingLineSegments<dim> temp=nA->glidePlaneIntersections();
////                        temp.updateWithGlidePlane(nB->glidePlane(0));
////
////                        if(temp.size())
////                        {
////                            std::cout<<"contractPoint case 10"<<std::endl;
////                            nA->set_P(temp.snap(0.5*(nA->get_P()+nB->get_P()))); THE PROBLE HERE IS THAT WE NEED TO SNAP TO GLIDEPLANE INTERSECTION!!!
////                            return this->contractSecond(nA,nB);
////                        }
////                        else
////                        {
////                            std::cout<<"contractPoint case 11"<<std::endl;
////                            return false;
////                        }
//                    }
//                    else if(!nA->glidePlaneIntersections().size() && nB->glidePlaneIntersections().size())
//                    {
//                        std::cout<<"contractPoint case 12"<<std::endl;
//                        return contract(nB,nA);
//                    }
//                    else
//                    {// both nodes confined by only one plane
//                        std::cout<<"contractPoint case 13aa"<<std::endl;
//                        
//                        assert(nA->glidePlanes().size()==1);
//                        assert(nB->glidePlanes().size()==1);
//                        
//                        GlidePlaneObserver<DislocationNetworkType>* const gpo(nA->glidePlane(0).glidePlaneObserver);
//                        const PlanePlaneIntersection<dim>& ppi(gpo->glidePlaneIntersection(&nA->glidePlane(0),&nB->glidePlane(0)));
//                        
//                        //                        std::cout<<"gpA: "<<nA->glidePlane(0).P.cartesian().transpose()<<std::endl;
//                        //                        std::cout<<"gpA: "<<nA->glidePlane(0).n.cartesian().normalized().transpose()<<std::endl;
//                        //
//                        //                        std::cout<<"gpB: "<<nB->glidePlane(0).P.cartesian().transpose()<<std::endl;
//                        //                        std::cout<<"gpB: "<<nB->glidePlane(0).n.cartesian().normalized().transpose()<<std::endl;
//                        
//                        
//                        if(ppi.type==PlanePlaneIntersection<dim>::COINCIDENT)
//                        {
//                            std::cout<<"contractPoint case 13"<<std::endl;
//                            nA->set_P(nA->snapToGlidePlaneIntersection(0.5*(nA->get_P()+nB->get_P()))); // this snaps to the glide planes to kill numerical errors
//                            return this->contractSecond(nA,nB);
//                        }
//                        else if(ppi.type==PlanePlaneIntersection<dim>::INCIDENT)
//                        {
//                            //                            std::cout<<nA->sID<<" "<<nB->sID<<std::endl;
//                            std::cout<<"contractPoint case 13a"<<std::endl;
//                            //                            std::cout<<"norm d="<<ppi.d.norm()<<std::endl;
//                            const double u=(0.5*(nA->get_P()+nB->get_P())-ppi.P).dot(ppi.d);
//                            nA->set_P(ppi.P+u*ppi.d);
//                            return this->contractSecond(nA,nB);
//                        }
//                        else
//                        {
//                            std::cout<<"contractPoint case 14"<<std::endl;
//                            return false;
//                        }
//                    }
//                }
//            }
//            else if(nA->isGlissile() && !nB->isGlissile())
//            {// a is glissile, b is sessile
//                if(nA->isOnBoundingBox())
//                {// a is a boundary node. Check is the bonuding box contains b
//                    if(nA->boundingBoxSegments().contains(nB->get_P()).first)
//                    {
//                        std::cout<<"contractPoint case 15"<<std::endl;
//                        return this->contractSecond(nB,nA);
//                    }
//                    else
//                    {
//                        std::cout<<"contractPoint case 16"<<std::endl;
//                        return false;
//                    }
//                }
//                else
//                {// a is an inner node. Check if _glidePlaneIntersections contains b
//                    if(nA->glidePlaneIntersections().size())
//                    {
//                        if(nA->glidePlaneIntersections().contains(nB->get_P()).first)
//                        {
//                            std::cout<<"contractPoint case 17"<<std::endl;
//                            return this->contractSecond(nB,nA);
//                        }
//                        else
//                        {
//                            std::cout<<"contractPoint case 17a"<<std::endl;
//                            return false;
//                        }
//                    }
//                    else
//                    {
//                        assert(nA->glidePlanes().size()==1);
//                        if(nA->glidePlane(0).contains(nB->get_P()))
//                        {
//                            std::cout<<"contractPoint case 17b"<<std::endl;
//                            return this->contractSecond(nB,nA);
//                        }
//                        else
//                        {
//                            std::cout<<"contractPoint case 17c"<<std::endl;
//                            return false;
//                        }
//                    }
//                }
//            }
//            else if(!nA->isGlissile() && nB->isGlissile())
//            {
//                std::cout<<"contractPoint case 18"<<std::endl;
//                return contract(nB,nA);
//            }
//            else
//            {// both nodes are sessile
//                if((nA->get_P()-nB->get_P()).squaredNorm()<FLT_EPSILON)
//                {// contract only if coincident
//                    std::cout<<"contractPoint case 19"<<std::endl;
//                    nA->set_P(0.5*(nA->get_P()+nB->get_P()));
//                    return this->contractSecond(nA,nB);
//                }
//                else
//                {
//                    std::cout<<"contractPoint case 20"<<std::endl;
//                    return false;
//                }
//            }
//            
//        }
        
        /**********************************************************************/
        const double& get_dt() const
        {/*!\returns the current time step increment dt
          */
            return dt;
        }
        
        /**********************************************************************/
        void set_dt(const double& dt_in,const double& vMax_in)
        {/*!\param[in] dt_in
          * Sets the time increment dt to dt_in
          */
            dt=dt_in;
            vMax=vMax_in;
        }
        
        /**********************************************************************/
        const double& get_totalTime() const
        {/*! The elapsed simulation time step in dimensionless units
          */
            return totalTime;
        }
        
        /**********************************************************************/
        void assembleAndSolve()
        {/*! Performs the following operatons:
          */
            const auto t0= std::chrono::system_clock::now();
#ifdef _OPENMP
            const size_t nThreads = omp_get_max_threads();
#else
            const size_t nThreads = 1;
#endif
            
            //! -1 Compute the interaction StressField between dislocation particles
            model::cout<<"		Computing stress field at quadrature points ("<<nThreads<<" threads)..."<<std::flush;
            this->template computeNeighborField<StressField>();
            model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<defaultColor<<std::endl;
            
            //! -2 Loop over DislocationSegments and assemble stiffness matrix and force vector
            model::cout<<"		Computing segment stiffness matrices and force vectors ("<<nThreads<<" threads)..."<<std::flush;
            const auto t2= std::chrono::system_clock::now();
            typedef void (LinkType::*LinkMemberFunctionPointerType)(void); // define type of Link member function
            LinkMemberFunctionPointerType Lmfp(&LinkType::assemble); // Lmfp is a member function pointer to Link::assemble
            this->parallelExecute(Lmfp);
            model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t2)).count()<<" sec]."<<defaultColor<<std::endl;
            
            //! -3 Loop over DislocationSubNetworks, assemble subnetwork stiffness matrix and force vector, and solve
            model::cout<<"		Assembling NetworkComponents and solving "<<std::flush;
            const auto t3= std::chrono::system_clock::now();
            
            if(DislocationNetworkComponentType::use_directSolver)
            {
                
#ifdef _MODEL_PARDISO_SOLVER_
                model::cout<<"(PardisoLDLT "<<nThreads<<" threads)..."<<std::flush;
                for (const auto& networkComponent : this->components())
                {
                    DislocationNetworkComponentType(*networkComponent.second).directSolve();
                }
#else
                model::cout<<"(SimplicialLDLT "<<nThreads<<" threads)..."<<std::flush;
#ifdef _OPENMP // SimplicialLDLT is not multi-threaded. So parallelize loop over NetworkComponents.
#pragma omp parallel for
                for (unsigned int k=0;k<this->components().size();++k)
                {
                    auto snIter(this->components().begin());
                    std::advance(snIter,k);
                    DislocationNetworkComponentType(*snIter->second).directSolve();
                }
#else
                for (const auto& networkComponent : this->components())
                {
                    DislocationNetworkComponentType(*networkComponent.second).directSolve();
                }
#endif
#endif
            }
            else // iterative solver
            {
                model::cout<<"(MINRES "<<nThreads<<" threads)..."<<std::flush;
#ifdef _OPENMP // MINRES is not multi-threaded. So parallelize loop over NetworkComponents.
#pragma omp parallel for
                for (unsigned int k=0;k<this->components().size();++k)
                {
                    auto snIter(this->components().begin());
                    std::advance(snIter,k);
                    DislocationNetworkComponentType(*snIter->second).iterativeSolve();
                }
#else
                for (const auto& networkComponent : this->components())
                {
                    DislocationNetworkComponentType(*networkComponent.second).iterativeSolve();
                }
#endif
            }
            
            model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t3)).count()<<" sec]."<<defaultColor<<std::endl;
        }
        
        /**********************************************************************/
        void updateQuadraturePoints()
        {
            model::cout<<"		Updating quadrature points... "<<std::flush;
            const auto t0=std::chrono::system_clock::now();
            
            // Clear DislocationParticles
            this->clearParticles(); // this also destroys all cells
            
            // Populate DislocationParticles
            for(auto& linkIter : this->links())
            {
                linkIter.second->updateQuadraturePoints(*this);
            }
            model::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
        }
        
        
        /**********************************************************************/
        //		void move(const double & dt_in, const double & dt_old)
        void move(const double & dt_in)
        {/*! Moves all nodes in the DislocationNetwork using the stored velocity and current dt
          */
            model::cout<<"		Moving DislocationNodes (dt="<<dt_in<< ")... "<<std::flush;
            const auto t0= std::chrono::system_clock::now();
            for (auto& nodeIter : this->nodes())
            {
                nodeIter.second->move(dt_in,DDtimeIntegrator<0>::dxMax);
            }
            model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<defaultColor<<std::endl;
        }
        
        /**********************************************************************/
        void runSteps()
        {/*! Runs Nsteps simulation steps
          */
            const auto t0= std::chrono::system_clock::now();
            while (runID<Nsteps)
            {
                model::cout<<std::endl; // leave a blank line
                singleStep();
            }
            updateQuadraturePoints(); // necessary if quadrature data are necessary in main
            model::cout<<greenBoldColor<<std::setprecision(3)<<std::scientific<<Nsteps<< " simulation steps completed in "<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" [sec]"<<defaultColor<<std::endl;
        }
        
        /**********************************************************************/
        void computePlasticDistortionRate()
        {
            _plasticDistortionRate.setZero();
            for (const auto& linkIter : this->links())
            {
                _plasticDistortionRate+= linkIter.second->plasticDistortionRate();
            }
        }
        
        /**********************************************************************/
        const MatrixDimD& plasticDistortionRate() const
        {
            return _plasticDistortionRate;
        }
        
        /**********************************************************************/
        const MatrixDimD& plasticDistortion() const
        {
            return _plasticDistortion;
        }
        
        /**********************************************************************/
        MatrixDimD plasticStrainRate() const
        {/*!\returns the plastic distortion rate tensor generated by this
          * DislocaitonNetwork.
          */
            //const MatrixDimD temp(plasticDistortionRate());
            return (_plasticDistortionRate+_plasticDistortionRate.transpose())*0.5;
        }
        
        /**********************************************************************/
        std::tuple<double,double,double> networkLength() const
        {/*!\returns the total line length of the DislocationNetwork. The return
          * value is a tuple, where the first value is the length of bulk glissile
          * dislocations, the second value is the length of bulk sessile
          * dislocations, and the third value is the length accumulated on
          * the mesh boundary.
          */
            double bulkGlissileLength(0.0);
            double bulkSessileLength(0.0);
            double boundaryLength(0.0);
            
            for (const auto& linkIter : this->links())
            {
                const double temp(linkIter.second->arcLength());
                if(linkIter.second->isBoundarySegment())
                {
                    boundaryLength+=temp;
                }
                else
                {
                    if(linkIter.second->isSessile())
                    {
                        bulkSessileLength+=temp;
                    }
                    else
                    {
                        bulkGlissileLength+=temp;
                        
                    }
                }
                
            }
            return std::make_tuple(bulkGlissileLength,bulkSessileLength,boundaryLength);
            
        }
        
        /**********************************************************************/
        const long int& runningID() const
        {/*! The current simulation step ID.
          */
            return runID;
        }
        
        /**********************************************************************/
        MatrixDimD stress(const VectorDim& P) const
        {/*!\param[in] P position vector
          * \returns The stress field generated by the DislocationNetwork at P
          *
          * Note:
          */
            SingleFieldPoint<StressField> fieldPoint(P,true);
            this->template computeField<SingleFieldPoint<StressField>,StressField>(fieldPoint);
            return fieldPoint.field();
        }
        
        /**********************************************************************/
        std::pair<bool,const Simplex<dim,dim>*> pointIsInsideMesh(const VectorDim& P0, const Simplex<dim,dim>* const guess) const
        {/*!\param[in] P0 position vector
          * \param[in] guess pointer of the Simplex where the search starts
          * \returns true if P0 is inside the mesh
          */
            std::pair<bool,const Simplex<dim,dim>*> temp(true,NULL);
            if (use_boundary)
            {
                temp=mesh.searchWithGuess(P0,guess);
            }
            return temp;
        }
        
        /**********************************************************************/
        static unsigned int& userOutputColumn()
        {
            return DislocationNetworkIO<DislocationNetworkType>::_userOutputColumn;
        }
        
        /**********************************************************************/
        const ParticleSystemType& particleSystem() const
        {
            return *this;
        }
        
        /**********************************************************************/
        DislocationNetworkIOType io()
        {
            return DislocationNetworkIOType(*this);
        }
        
        /**********************************************************************/
        DislocationNetworkIOType io() const
        {
            return DislocationNetworkIOType(*this);
        }
        
    };
    
}
#endif

//        /**********************************************************************/
//        void make_bndNormals()
//        {
//            for(auto& node : this->nodes())
//            {
//                node.second->make_bndNormal();
//            }
//        }

//        /**********************************************************************/
//        void output(const size_t& fileID) const
//        {/*! Outputs DislocationNetwork data
//          */
//            if (!(runID%DislocationNetworkIO<DislocationNetworkType>::outputFrequency))
//            {
//                const auto t0=std::chrono::system_clock::now();
//#ifdef _MODEL_DD_MPI_
//                if(ModelMPIbase::mpiRank()==0)
//                {
//                    DislocationNetworkIO<DislocationNetworkType>::output(*this,fileID);
//                }
//#else
//                DislocationNetworkIO<DislocationNetworkType>::output(*this,fileID);
//#endif
//                model::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
//            }
//        }

//        /**********************************************************************/
//        void read(const std::string& inputDirectoryName_in, std::string inputFileName)
//        { // TO DO: move this to DislocationNetworkIO.h
//
//            std::ostringstream fullName;
//            fullName<<inputDirectoryName_in<<inputFileName;
//
//            model::cout<<greenBoldColor<<"Reading "<<fullName.str()<<"..."<<defaultColor<<std::endl;
//
//
//            // Create a file-reader object
//            EigenDataReader EDR;
//
//            // IO
//            EDR.readScalarInFile(fullName.str(),"outputFrequency",DislocationNetworkIO<DislocationNetworkType>::outputFrequency);
//            EDR.readScalarInFile(fullName.str(),"outputBinary",DislocationNetworkIO<DislocationNetworkType>::outputBinary);
//            EDR.readScalarInFile(fullName.str(),"outputGlidePlanes",DislocationNetworkIO<DislocationNetworkType>::outputGlidePlanes);
//            EDR.readScalarInFile(fullName.str(),"outputSpatialCells",DislocationNetworkIO<DislocationNetworkType>::outputSpatialCells);
//            EDR.readScalarInFile(fullName.str(),"outputPKforce",DislocationNetworkIO<DislocationNetworkType>::outputPKforce);
//            EDR.readScalarInFile(fullName.str(),"outputMeshDisplacement",DislocationNetworkIO<DislocationNetworkType>::outputMeshDisplacement);
//            EDR.readScalarInFile(fullName.str(),"outputElasticEnergy",DislocationNetworkIO<DislocationNetworkType>::outputElasticEnergy);
//
//            EDR.readScalarInFile(fullName.str(),"outputPlasticDistortion",DislocationNetworkIO<DislocationNetworkType>::outputPlasticDistortion);
//            if(DislocationNetworkIO<DislocationNetworkType>::outputPlasticDistortion)
//            {
//                DislocationNetworkIO<DislocationNetworkType>::_userOutputColumn+=9;
//            }
//            EDR.readScalarInFile(fullName.str(),"outputPlasticDistortionRate",DislocationNetworkIO<DislocationNetworkType>::outputPlasticDistortionRate);
//            if(DislocationNetworkIO<DislocationNetworkType>::outputPlasticDistortionRate)
//            {
//                DislocationNetworkIO<DislocationNetworkType>::_userOutputColumn+=9;
//            }
//
//            EDR.readScalarInFile(fullName.str(),"outputDislocationLength",DislocationNetworkIO<DislocationNetworkType>::outputDislocationLength);
//            if(DislocationNetworkIO<DislocationNetworkType>::outputDislocationLength)
//            {
//                DislocationNetworkIO<DislocationNetworkType>::_userOutputColumn+=3;
//            }
//
//            EDR.readScalarInFile(fullName.str(),"outputQuadratureParticles",DislocationNetworkIO<DislocationNetworkType>::outputQuadratureParticles);
//
//            // Parametrization exponent
//            EDR.readScalarInFile(fullName.str(),"parametrizationExponent",LinkType::alpha);
//            assert((LinkType::alpha)>=0.0 && "parametrizationExponent MUST BE >= 0.0");
//            assert((LinkType::alpha)<=1.0 && "parametrizationExponent MUST BE <= 1.0");
//
//            // Temperature. Make sure you initialize before calling Material<Isotropic>::select()
//            EDR.readScalarInFile(fullName.str(),"temperature",Material<Isotropic>::T); // temperature
//
//            // Material and crystal orientation
//            unsigned int materialZ;
//            EDR.readScalarInFile(fullName.str(),"material",materialZ); // material by atomic number Z
//            Material<Isotropic>::select(materialZ);
//
//            // quadPerLength
//            EDR.readScalarInFile(fullName.str(),"quadPerLength",LinkType::quadPerLength); // quadPerLength
//
//            // core size
//            EDR.readScalarInFile(fullName.str(),"coreSize",StressField::a); // core-width
//            assert((StressField::a)>0.0 && "coreSize MUST BE > 0.");
//            StressField::a2=StressField::a*StressField::a;
//
//            // multipole expansion
//            double cellSize(0.0);
//            EDR.readScalarInFile(fullName.str(),"dislocationCellSize",cellSize);
//            SpatialCellObserverType::setCellSize(cellSize);
//            EDR.readScalarInFile(fullName.str(),"use_DisplacementMultipole",DislocationDisplacement<dim>::use_multipole);
//            EDR.readScalarInFile(fullName.str(),"use_StressMultipole",DislocationStress<dim>::use_multipole);
//            EDR.readScalarInFile(fullName.str(),"use_EnergyMultipole",DislocationEnergy<dim>::use_multipole);
//
//            // Eternal Stress
//            EDR.readMatrixInFile(fullName.str(),"externalStress",externalStress);
//
//            //dt=0.0;
//            EDR.readScalarInFile(fullName.str(),"timeIntegrationMethod",timeIntegrationMethod);
//            EDR.readScalarInFile(fullName.str(),"dxMax",DDtimeIntegrator<0>::dxMax);
//            assert(DDtimeIntegrator<0>::dxMax>0.0);
//            //            EDR.readScalarInFile(fullName.str(),"shearWaveSpeedFraction",shearWaveSpeedFraction);
//            //            assert(shearWaveSpeedFraction>=0.0);
//            EDR.readScalarInFile(fullName.str(),"use_velocityFilter",NodeType::use_velocityFilter);
//            EDR.readScalarInFile(fullName.str(),"velocityReductionFactor",NodeType::velocityReductionFactor);
//            assert(NodeType::velocityReductionFactor>0.0 && NodeType::velocityReductionFactor<=1.0);
//
//            // Restart
//            EDR.readScalarInFile(fullName.str(),"startAtTimeStep",runID);
//            //            VertexReader<'F',201,double> vReader;
//            IDreader<'F',1,200,double> vReader;
//            Eigen::Matrix<double,1,200> temp(Eigen::Matrix<double,1,200>::Zero());
//
//
//            if (vReader.isGood(0,true))
//            {
//
//                vReader.read(0,true);
//
//                if(runID<0)
//                {
//                    if(vReader.size())
//                    {
//                        runID=vReader.rbegin()->first;
//                        temp=Eigen::Map<Eigen::Matrix<double,1,200>>(vReader.rbegin()->second.data());
//                    }
//                    else
//                    {
//                        runID=0;
//                    }
//                }
//                else
//                {
//                    const auto iter=vReader.find(runID);
//                    if(iter!=vReader.end())
//                    {// runID has been found
//                        temp=Eigen::Map<Eigen::Matrix<double,1,200>>(iter->second.data());
//                    }
//                    else
//                    {
//                        assert(0 && "runID NOT FOUND IN F/F_0.txt");
//                    }
//                }
//
//                dt=temp(1);
//
//            }
//            else
//            {
//                model::cout<<"could not read runID from F/F_0.txt"<<std::endl;
//                runID=0;
//                dt=10.0;
//            }
//
//            model::cout<<"dt="<<dt<<std::endl;
//
//            size_t curCol=0;
//            totalTime=temp(curCol);
//            curCol+=2;
//
//            if (DislocationNetworkIO<DislocationNetworkType>::outputPlasticDistortion)
//            {
//                std::cout<<"reading PD"<<std::endl;
//
//                for(int r=0;r<3;++r)
//                {
//                    for(int c=0;c<3;++c)
//                    {
//                        _plasticDistortion(r,c)=temp(curCol);
//                        curCol+=1;
//                    }
//                }
//            }
//
//            model::cout<<"starting at time step "<<runID<<std::endl;
//            model::cout<<"totalTime= "<<totalTime<<std::endl;
//            model::cout<<"plasticDistortion=\n "<<_plasticDistortion<<std::endl;
//
//            // time-stepping
//
//            //            EDR.readScalarInFile(fullName.str(),"useImplicitTimeIntegration",useImplicitTimeIntegration);
//            EDR.readScalarInFile(fullName.str(),"use_directSolver_DD",DislocationNetworkComponentType::use_directSolver);
//
//
//            EDR.readScalarInFile(fullName.str(),"Nsteps",Nsteps);
//            assert(Nsteps>=0 && "Nsteps MUST BE >= 0");
//
//            //            EDR.readScalarInFile(fullName.str(),"timeWindow",timeWindow);
//            //            assert(timeWindow>=0.0 && "timeWindow MUST BE >= 0");
//
//            // Check balance
////            EDR.readScalarInFile(fullName.str(),"check_balance",check_balance);
//
//            // JUNCTION FORMATION
//            EDR.readScalarInFile(fullName.str(),"use_junctions",use_junctions);
//            //            EDR.readScalarInFile(fullName.str(),"collisionTol",DislocationJunctionFormation<DislocationNetworkType>::collisionTol);
//
//            // VERTEX REDISTRIBUTION
//            EDR.readScalarInFile(fullName.str(),"use_redistribution",DislocationNetworkRemesh<DislocationNetworkType>::use_redistribution);
//            EDR.readScalarInFile(fullName.str(),"Lmin",DislocationNetworkRemesh<DislocationNetworkType>::Lmin);
//            assert(DislocationNetworkRemesh<DislocationNetworkType>::Lmin>=0.0);
//            assert(DislocationNetworkRemesh<DislocationNetworkType>::Lmin>=2.0*DDtimeIntegrator<0>::dxMax && "YOU MUST CHOOSE Lmin>2*dxMax.");
//            EDR.readScalarInFile(fullName.str(),"Lmax",DislocationNetworkRemesh<DislocationNetworkType>::Lmax);
//            assert(DislocationNetworkRemesh<DislocationNetworkType>::Lmax>DislocationNetworkRemesh<DislocationNetworkType>::Lmin);
//
//            // Cross-Slip
//            //            EDR.readScalarInFile(fullName.str(),"use_crossSlip",DislocationCrossSlip<DislocationNetworkType>::use_crossSlip);
//            //            if(DislocationCrossSlip<DislocationNetworkType>::use_crossSlip)
//            //            {
//            //                EDR.readScalarInFile(fullName.str(),"crossSlipDeg",DislocationCrossSlip<DislocationNetworkType>::crossSlipDeg);
//            //                assert(DislocationCrossSlip<DislocationNetworkType>::crossSlipDeg>=0.0 && DislocationCrossSlip<DislocationNetworkType>::crossSlipDeg <= 90.0 && "YOU MUST CHOOSE 0.0<= crossSlipDeg <= 90.0");
//            //                //                EDR.readScalarInFile(fullName.str(),"crossSlipLength",DislocationCrossSlip<DislocationNetworkType>::crossSlipLength);
//            //                //                assert(DislocationCrossSlip<DislocationNetworkType>::crossSlipLength>=DislocationNetworkRemesh<DislocationNetworkType>::Lmin && "YOU MUST CHOOSE crossSlipLength>=Lmin.");
//            //            }
//
//            // Mesh and BVP
//            EDR.readScalarInFile(fullName.str(),"use_boundary",use_boundary);
//            if (use_boundary)
//            {
//                //                EDR.readScalarInFile(fullName.str(),"use_meshRegions",use_meshRegions);
//
//                int meshID(0);
//                EDR.readScalarInFile(fullName.str(),"meshID",meshID);
//                mesh.readMesh(meshID);
//                assert(mesh.simplices().size() && "MESH IS EMPTY.");
//
//                // Initialize Polycrystal
//                poly.init(*this,fullName.str());
//
//                EDR.readScalarInFile(fullName.str(),"use_virtualSegments",use_virtualSegments);
//                if(use_virtualSegments)
//                {
//                    EDR.readScalarInFile(fullName.str(),"virtualSegmentDistance",LinkType::virtualSegmentDistance);
//                }
//
//                EDR.readScalarInFile(fullName.str(),"use_bvp",use_bvp);
//                if(use_bvp)
//                {
//                    EDR.readScalarInFile(fullName.str(),"use_directSolver_FEM",bvpSolver.use_directSolver);
//                    EDR.readScalarInFile(fullName.str(),"solverTolerance",bvpSolver.tolerance);
//                    bvpSolver.init(*this);
//                }
//            }
//            else{ // no boundary is used, DislocationNetwork is in inifinite medium
//                use_bvp=0;	// never comupute boundary correction
//            }
//
//            // Grain Boundary flags
//            //            EDR.readScalarInFile(fullName.str(),"use_GBdissociation",GrainBoundaryDissociation<DislocationNetworkType>::use_GBdissociation);
//            //            EDR.readScalarInFile(fullName.str(),"use_GBtransmission",GrainBoundaryTransmission<DislocationNetworkType>::use_GBtransmission);
//            //            EDR.readScalarInFile(fullName.str(),"use_GBdislocations",GrainBoundary<dim>::use_GBdislocations);
//
//            // Read Vertex and Edge information
//            io().readVertices(runID); // this requires mesh to be up-to-date
//            io().readEdges(runID);    // this requires mesh to be up-to-date
//
//            //#ifdef DislocationNucleationFile
//            //            EDR.readScalarInFile(fullName.str(),"nucleationFreq",nucleationFreq);
//            //#endif
//
//#ifdef _MODEL_MPI_
//            // Avoid that a processor starts writing before other are reading
//            MPI_Barrier(MPI_COMM_WORLD);
//#endif
//
//            // Initializing configuration
//            move(0.0);	// initial configuration
//        }

