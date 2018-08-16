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

// valgrind --leak-check=full --show-leak-kinds=all ./DDomp

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
#include <Eigen/StdVector>
#include <Eigen/StdDeque>


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
#include <model/DislocationDynamics/CrossSlip/DislocationCrossSlip.h>
#include <model/DislocationDynamics/Materials/Material.h>
#include <model/DislocationDynamics/IO/DislocationNetworkIO.h>
#include <model/DislocationDynamics/ElasticFields/DislocationParticle.h>
#include <model/DislocationDynamics/ElasticFields/DislocationStress.h>
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
#include <model/DislocationDynamics/ElasticFields/EshelbyInclusion.h>


#ifndef ExternalLoadControllerFile
#define ExternalLoadControllerFile <model/DislocationDynamics/ExternalLoadControllers/DummyExternalLoadController.h>
#endif
#include ExternalLoadControllerFile

//#include <model/DislocationDynamics/ExternalLoadController.h>
#include <model/DislocationDynamics/DislocationInjector.h>

namespace model
{
    
    template <int dim>
    struct DislocationNetworkBase
    {// Class storing objects that need to be destroyed last
        
        SimplicialMesh<dim> mesh;
        Polycrystal<dim> poly;
        BVPsolver<dim,2> bvpSolver;
        
        DislocationNetworkBase() :
        /* init list  */ poly(mesh),
        /* init list  */ bvpSolver(mesh)
        {
            
        }
        
    };
    
    template <int _dim, short unsigned int corder, typename InterpolationType>
    class DislocationNetwork : public DislocationNetworkBase<_dim>, // must be first in inheritance tree
    /* base                 */ public GlidePlaneObserver<_dim>,
    /* base                 */ public LoopNetwork<DislocationNetwork<_dim,corder,InterpolationType> >,
    /* base                 */ public ParticleSystem<DislocationParticle<_dim> >,
    /* base                 */ public std::map<size_t,EshelbyInclusion<_dim>>
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
        typedef GlidePlaneObserver<dim> GlidePlaneObserverType;
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
        typedef Polycrystal<dim> PolycrystalType;
        typedef ExternalLoadController<dim> ExternalLoadControllerType;
        //        enum {NdofXnode=NodeType::NdofXnode};
        typedef std::map<size_t,EshelbyInclusion<_dim>> EshelbyInclusionContainerType;
        
        typedef NetworkLinkObserver<LinkType> NetworkLinkObserverType;
        typedef typename NetworkLinkObserverType::LinkContainerType NetworkLinkContainerType;
        
        
        int timeIntegrationMethod;
        //        bool use_analyticalStraightStress;
        //        bool use_junctions;
        int maxJunctionIterations;
        long int runID;
        double totalTime;
        double dt;
        double vMax;
        int Nsteps;
        MatrixDimD _plasticDistortion;
        MatrixDimD _plasticDistortionRate;
        int ddSolverType;
        bool computeDDinteractions;
        int crossSlipModel;
        bool use_boundary;
        unsigned int use_bvp;
        bool use_virtualSegments;
        //        SimplicialMesh<dim> mesh;
        //        PolycrystalType poly;
        //        BvpSolverType bvpSolver;
        // MatrixDimD externalStress;
        bool use_externalStress;
        //        bool use_userStress;
        bool use_extraStraightSegments;
        ExternalLoadControllerType extStressController;
        std::deque<StressStraight<dim>,Eigen::aligned_allocator<StressStraight<dim>>> ssdeq;
        
        int  outputFrequency;
        bool outputBinary;
        bool outputGlidePlanes;
        bool outputSpatialCells;
        bool outputPKforce;
        bool outputElasticEnergy;
        bool outputMeshDisplacement;
        bool outputFEMsolution;
        bool outputDislocationLength;
        bool outputPlasticDistortion;
        bool outputPlasticDistortionRate;
        bool outputQuadratureParticles;
        bool outputLinkingNumbers;
        bool outputLoopLength;
        bool outputSegmentPairDistances;
        unsigned int _userOutputColumn;
        bool use_stochasticForce;
        int dislocationImages_x;
        int dislocationImages_y;
        int dislocationImages_z;
        double surfaceAttractionDistance;
        std::string folderSuffix;
        
#ifdef DislocationNucleationFile
#include DislocationNucleationFile
#endif
        
    private:
        
        /**********************************************************************/
        void updateLoadControllers()
        {/*! Updates bvpSolver using the stress and displacement fields of the
          *  current DD configuration.
          */
            const int quadraturePerTriangle=37;
            if(use_bvp)
            {
                if (!(runID%use_bvp))
                {// enter the if statement if use_bvp!=0 and runID is a multiple of use_bvp
                    model::cout<<"		Updating elastic bvp... "<<std::endl;
                    this->bvpSolver.template assembleAndSolve<DislocationNetworkType,quadraturePerTriangle>(*this);
                }
            }
            if (use_externalStress)
            {
                extStressController.update(*this);
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
            
            //            if(NodeType::use_velocityFilter)
            //            {
            //                assert(0 && "velocityFilter not implemented yet.");
            //                //            for(auto& node : this->nodes())
            //                //            {
            //                //                node.second->applyVelocityFilter(vMax);
            //                //            }
            //            }
            
            
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
            updateLoadControllers();
            
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
            
            if(outputElasticEnergy)
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
            DislocationCrossSlip<DislocationNetworkType>(*this);
            
            
            //                        GrainBoundaryTransmission<DislocationNetworkType>(*this).transmit();
            GrainBoundaryTransmission<DislocationNetworkType>(*this).directTransmit();
            
            //
            //            GrainBoundaryDissociation<DislocationNetworkType>(*this).dissociate();
            
            //            poly.reConnectGrainBoundarySegments(*this); // this makes stressGauss invalid, so must follw other GB operations
            
            
            //! 11- detect loops that shrink to zero and expand as inverted loops
            //            DislocationNetworkRemesh<DislocationNetworkType>(*this).loopInversion(dt);
            
            //! 12- Form Junctions
            DislocationJunctionFormation<DislocationNetworkType>(*this).formJunctions(maxJunctionIterations,DDtimeIntegrator<0>::dxMax);
            
            //            // Remesh may contract juncitons to zero lenght. Remove those juncitons:
            //            DislocationJunctionFormation<DislocationNetworkType>(*this).breakZeroLengthJunctions();
            
            //! 13- Node redistribution
            DislocationNetworkRemesh<DislocationNetworkType>(*this).remesh(runID);
            
//            mergeLoopsAtNodes();
            
            //            DislocationInjector<DislocationNetworkType>(*this).insertRandomStraightDislocation();
            
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
        //        /* init list  */ use_analyticalStraightStress(true),
        //        /* init list  */ use_junctions(false),
        maxJunctionIterations(1),
        /* init list  */ runID(0),
        /* init list  */ totalTime(0.0),
        /* init list  */ dt(0.0),
        /* init list  */ vMax(0.0),
        /* init list  */ Nsteps(0),
        /* init list  */ _plasticDistortion(MatrixDimD::Zero()),
        /* init list  */ _plasticDistortionRate(MatrixDimD::Zero()),
        ///* init list  */ externalStress(MatrixDimD::Zero()),
        /* init list  */ ddSolverType(0),
        /* init list  */ computeDDinteractions(true),
        /* init list  */ crossSlipModel(0),
        /* init list  */ use_boundary(false),
        /* init list  */ use_bvp(0),
        /* init list  */ use_virtualSegments(true),
        //        /* init list  */ poly(mesh),
        //        /* init list  */ bvpSolver(mesh),
        /* init list  */ use_externalStress(false),
        /* init list  */ use_extraStraightSegments(false),
        /* init list  */ extStressController(),
        /* init list  */ ssdeq(),
        /* init list  */ outputFrequency(1),
        /* init list  */ outputBinary(0),
        /* init list  */ outputGlidePlanes(false),
        /* init list  */ outputSpatialCells(false),
        /* init list  */ outputPKforce(false),
        /* init list  */ outputElasticEnergy(false),
        /* init list  */ outputMeshDisplacement(false),
        /* init list  */ outputFEMsolution(false),
        /* init list  */ outputDislocationLength(false),
        /* init list  */ outputPlasticDistortion(false),
        /* init list  */ outputPlasticDistortionRate(false),
        /* init list  */ outputQuadratureParticles(false),
        /* init list  */ outputLinkingNumbers(false),
        /* init list  */ outputLoopLength(false),
        /* init list  */ outputSegmentPairDistances(false),
        /* init list  */ _userOutputColumn(3),
        /* init list  */ use_stochasticForce(false),
        /* init list  */ surfaceAttractionDistance(0.0),
        //        /* init list  */ use_userStress(false),
        /* init list  */ folderSuffix("")
        {
            
            if(argc>1)
            {
                folderSuffix=argv[1];
                std::cout<<"folderSuffix="<<folderSuffix<<std::endl;
            }
            
            ParticleSystemType::initMPI(argc,argv);
            io().read("./","DDinput.txt");
            
            // Initializing configuration
            move(0.0);	// initial configuration
        }
        
        /**********************************************************************/
        const EshelbyInclusionContainerType& eshelbyInclusions() const
        {
            return *this;
        }

        EshelbyInclusionContainerType& eshelbyInclusions()
        {
            return *this;
        }

        
        /**********************************************************************/
        bool contract(std::shared_ptr<NodeType> nA,
                      std::shared_ptr<NodeType> nB)
        {
            return DislocationNodeContraction<DislocationNetworkType>(*this).contract(nA,nB);
        }
        
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
#ifdef _OPENMP
            const size_t nThreads = omp_get_max_threads();
#else
            const size_t nThreads = 1;
#endif
            
            //! -1 Compute the interaction StressField between dislocation particles
            
            
            if(corder==0)
            {// For straight segments use analytical expression of stress field
                const auto t0= std::chrono::system_clock::now();
                model::cout<<"		Collecting StressStraight objects: "<<std::flush;
                
                std::deque<StressStraight<dim>,Eigen::aligned_allocator<StressStraight<dim>>> straightSegmentsDeq;
                size_t currentSize=0;
                if(computeDDinteractions)
                {
                    for(const auto& link : this->networkLinks())
                    {
                        link.second->addToStressStraight(straightSegmentsDeq);
                    }
                    
                    
                    
                    currentSize=straightSegmentsDeq.size();
                    
                    const VectorDim meshSize(this->mesh.xMax()-this->mesh.xMin());
                    
                    for(int i=-dislocationImages_x;i<=dislocationImages_x;++i)
                    {
                        for(int j=-dislocationImages_y;j<=dislocationImages_y;++j)
                        {
                            for(int k=-dislocationImages_z;k<=dislocationImages_z;++k)
                            {
                                
                                const Eigen::Matrix<int,3,1> cellID((Eigen::Matrix<int,3,1>()<<i,j,k).finished());
                                
                                if( cellID.squaredNorm()!=0) //skip current cell
                                {
                                    for (size_t c=0;c<currentSize;++c)
                                    {
                                        const VectorDim P0=straightSegmentsDeq[c].P0+(meshSize.array()*cellID.cast<double>().array()).matrix();
                                        const VectorDim P1=straightSegmentsDeq[c].P1+(meshSize.array()*cellID.cast<double>().array()).matrix();
                                        
                                        straightSegmentsDeq.emplace_back(P0,P1,straightSegmentsDeq[c].b);
                                    }
                                }
                                
                                
                                
                            }
                        }
                    }
                    
                    
                }
                model::cout<< straightSegmentsDeq.size()<<" straight segments ("<<currentSize<<"+"<<straightSegmentsDeq.size()-currentSize<<" images)"<<std::flush;
                model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<defaultColor<<std::endl;
                
                const auto t1= std::chrono::system_clock::now();
                model::cout<<"		Computing analytical stress field at quadrature points ("<<nThreads<<" threads) "<<std::flush;
#ifdef _OPENMP
                //                const size_t nThreads = omp_get_max_threads();
                EqualIteratorRange<typename NetworkLinkContainerType::iterator> eir(this->links().begin(),this->links().end(),nThreads);
#pragma omp parallel for
                for(size_t thread=0;thread<eir.size();thread++)
                {
                    for (auto& linkIter=eir[thread].first;linkIter!=eir[thread].second;linkIter++)
                    {
                        linkIter->second->assemble(straightSegmentsDeq);
                    }
                }
#else
                for (auto& linkIter : this->links())
                {
                    linkIter.second->assemble(straightSegmentsDeq);
                }
#endif
                model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t1)).count()<<" sec]."<<defaultColor<<std::endl;
                
                
                
                
                
                
                //#ifdef _OPENMP
                //#pragma omp parallel for
                //#endif
                //                for (unsigned int k=0; k<this->particleSystem().size();++k)
                //                {
                //                    if(this->particleSystem()[k].template fieldPointBase<StressField>().enabled)
                //                    {
                //                        this->particleSystem()[k].template fieldPointBase<StressField>() += StressField::addSourceContribution(this->particleSystem()[k],straightSegmentsDeq);
                //                    }
                //                }
                //                model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<defaultColor<<std::endl;
            }
            else
            {// For curved segments use quandrature integration of stress field
                
                const auto t0= std::chrono::system_clock::now();
                
                if(computeDDinteractions)
                {
                    
                    if(dislocationImages_x!=0 || dislocationImages_y!=0 || dislocationImages_z!=0)
                    {
                        assert(0 && "FINISH HERE");
                    }
                    
                    model::cout<<"		Computing numerical stress field at quadrature points ("<<nThreads<<" threads)..."<<std::flush;
                    if (use_extraStraightSegments)
                    {
                        this->template computeNeighborField<StressField>(ssdeq);
                    }
                    else
                    {
                        this->template computeNeighborField<StressField>();
                    }
                    model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<defaultColor<<std::endl;
                }
                
                //! -2 Loop over DislocationSegments and assemble stiffness matrix and force vector
                const auto t1= std::chrono::system_clock::now();
                model::cout<<"		Computing segment stiffness matrices and force vectors ("<<nThreads<<" threads)..."<<std::flush;
                typedef void (LinkType::*LinkMemberFunctionPointerType)(void); // define type of Link member function
                LinkMemberFunctionPointerType Lmfp(&LinkType::assemble); // Lmfp is a member function pointer to Link::assemble
                this->parallelExecute(Lmfp);
                model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t1)).count()<<" sec]."<<defaultColor<<std::endl;
                
                
            }
            
            //! -3 Loop over DislocationSubNetworks, assemble subnetwork stiffness matrix and force vector, and solve
            const auto t2= std::chrono::system_clock::now();
            model::cout<<"		Assembling NetworkComponents and solving "<<std::flush;
            
            
            
            switch (ddSolverType)
            {
                case 1: // iterative solver
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
                    break;
                }
                    
                case 2: // direct solver
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
                    break;
                }
                    
                default: // lumped solver
                {
                    model::cout<<"(lumpedSolver "<<nThreads<<" threads)..."<<std::flush;
#ifdef _OPENMP
#pragma omp parallel for
                    for (unsigned int k=0;k<this->components().size();++k)
                    {
                        auto snIter(this->components().begin());
                        std::advance(snIter,k);
                        DislocationNetworkComponentType(*snIter->second).lumpedSolve(runID);
                    }
#else
                    for (const auto& networkComponent : this->components())
                    {
                        DislocationNetworkComponentType(*networkComponent.second).lumpedSolve(runID);
                    }
#endif
                    break;
                }
            }
            
            //            if(DislocationNetworkComponentType::use_directSolver)
            //            {
            //
            //
            //            }
            //            else // iterative solver
            //            {
            //
            //            }
            
            model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t2)).count()<<" sec]."<<defaultColor<<std::endl;
        }
        
        /**********************************************************************/
        void updateQuadraturePoints()
        {
            if(corder>0)
            {// quadrature points only used for curved segments
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
        }
        
        
        /**********************************************************************/
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
        
//        /**********************************************************************/
//        void mergeLoopsAtNodes()
//        {
//            
//            std::map<std::pair<size_t,size_t>,std::set<size_t>> loopPairMap;
//            
//            for(const auto& node : this->nodes())
//            {
//                const auto nodeLoops=node.second->loops();
//                
//                for(const auto& loop1 : nodeLoops)
//                {
//                    for(const auto& loop2 : nodeLoops)
//                    {
//                        if(loop1!=loop2)
//                        {
//                            if(   ((loop1->Burgers()+loop2->Burgers()).norm()<FLT_EPSILON || (loop1->Burgers()-loop2->Burgers()).norm()<FLT_EPSILON)
//                               && ((loop1->glidePlane.unitNormal+loop2->glidePlane.unitNormal).norm()<FLT_EPSILON || (loop1->glidePlane.unitNormal-loop2->glidePlane.unitNormal).norm()<FLT_EPSILON)
//                               && (loop1->grain.grainID==loop2->grain.grainID)
//                               )
//                            {
//                                
//                                std::pair<size_t,size_t> key=std::make_pair(std::min(loop1->sID,loop2->sID),std::max(loop1->sID,loop2->sID));
//                                
//                                loopPairMap[key].insert(node.second->sID);
//                            }
//                            
//                        }
//                    }
//                }
//                
//            }
//            
//            
//            for(const auto& pair : loopPairMap)
//            {
//                std::cout<<"loops "<<pair.first.first<<" "<<pair.first.second<<" meet at:"<<std::endl;
//                for(const int& nodeID : pair.second)
//                {
//                    std::cout<<"meetNode "<<nodeID<<std::endl;
//                }
//            }
//            
//        }
        
        /**********************************************************************/
        const long int& runningID() const
        {/*! The current simulation step ID.
          */
            return runID;
        }
        
        
        /**********************************************************************/
        const unsigned int& userOutputColumn()  const
        {
            return _userOutputColumn;
        }
        
        /**********************************************************************/
        MatrixDimD stress(const VectorDim& P) const
        {/*!\param[in] P position vector
          * \returns The stress field generated by the DislocationNetwork at P
          *
          * Note:
          */
            
            assert(false && "REWORK THIS FOR STRAIGHT SEGMENTS");
            
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
                temp=this->mesh.searchWithGuess(P0,guess);
            }
            return temp;
        }
        
        /**********************************************************************/
        const ParticleSystemType& particleSystem() const
        {
            return *this;
        }
        
        /**********************************************************************/
        ParticleSystemType& particleSystem()
        {
            return *this;
        }
        
        /**********************************************************************/
        DislocationNetworkIOType io()
        {
            return DislocationNetworkIOType(*this,folderSuffix);
        }
        
        /**********************************************************************/
        DislocationNetworkIOType io() const
        {
            return DislocationNetworkIOType(*this,folderSuffix);
        }
        
    };
    
}
#endif


