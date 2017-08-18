/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2011 by Benjamin Ramirez <ramirezbrf@gmail.com>,
 * Copyright (C) 2011 by Mamdouh Mohamed  <msm07d@fsu.edu>,
 * Copyright (C) 2011 by Tamer Crsoby     <tamercrosby@gmail.com>,
 * Copyright (C) 2011 by Can Erel         <canerel55@gmail.com>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */



// profiling
// http://stackoverflow.com/questions/375913/what-can-i-use-to-profile-c-code-in-linux/378024#378024
// valgrind --tool=callgrind ./(Your binary)
// It will generate a file called callgrind.out.x. You can then use kcachegrind tool to read this file. It will give you a graphical analysis of things with results like which lines cost how much.

// http://code.google.com/p/jrfonseca/wiki/Gprof2Dot

// valgrind --tool=callgrind ./DDserial
// KCachegrind callgrind.out.1378

// BEING MODIFIED
// 3 crossSlip

// TO DO
// -5 use modelMacros instead of assert
// -4 CHANGE definition of modelMacros in debug mode BECAUSE assert(x) CAN BE REMOVED BY NDEBUG!!! MUST USE mode::assert_fail
// -3 MAKE isIsolated and isBalanced data members of NetworkNode that are modified by addToNeighbors and removeFromNeigbors
// There is a bug with implicit integration. Nodes may exit the domain boundary.

// -1 - remove DislocationEnergyRules and EdgeConfig !!!
// 0 - Finish DepthFirst class. Don't allow to search/execute if N=0, so that N=1 means that node only, n=2 means first neighbor. Or change SpineNodeBase_CatmullRom::TopologyChangeActions
// 1- Implement operator << SpatialCell
// 2- remove AddressBook, wherever possible, Done in Node chain
// 40 - clean MultiExpand in Network Layer, remove get_r from there
// 35- Simplify Neighborhood structure
// 18- Should define linear=1, quadratic=2, cubic=3 and use polyDegree instead of corder. Put corder in SplineEnums
// 37- IS PLANAR SHOULD RETURN 0 IF IS A LINE!!!!! CHANGE ALSO IN SPLINESEGMENTBASE
// 38- IS PLANAR SHOULD RETURN the normal as pair<bool,normal>
// 37- NetworkNode, initializations from expansion and contraction
// IMPLEMENT NEIGHBOR ITERATORS IN NETWORKNODE
// CHANGE CONST VERSION OF EDGEFINDER/VERTEXFINDER PASS this IN MEMBER FUNCTION, REMOVE TEMPLATE SPECIALIZATION AND MAKE STATIC FUNCTIONS
// 9- cellSize should depend on applied load. Or better the number of cell neighbors used in each cell should depend on the applied stress to that cell
// 1- BIN WRITER/READER. DEFINE DISLOCATIONSEGMENT::OUTDATA as Eigen::Matrix<double,1,9>. Then in Network add the function friend void operator<< (SequentialBinFile<template Char,???>, ...)


// ISSUES
// 3- SplineIntersection::planePlaneType uses the wrong tolerance (10x)
// 10- SplineNetworkBase_common :		assert(iter->second.size()==2);
//     SplineIntersection: //			assert(T0.norm()>tol && "SplineIntersection<3,3,1,2>: T0 too small.");  // NOT RIGHT FOR DIPOLAR LOOPS
//									assert(T1.norm()>tol && "SplineIntersection<3,3,1,2>: T1 too small.");	// NOT RIGHT FOR DIPOLAR LOOPS
// 14- If T0=0 or T1=0, then rl(0)=0/0=NaN !!!!! This is not true since rl still tends to a finite vector. Remove class Parametric curve and implement special case of rl at 0 and 1 for vanishing nodal tangents


#ifndef model_DISLOCATIONNETWORK_H_
#define model_DISLOCATIONNETWORK_H_

//#ifdef _MODEL_DD_MPI_
//#define _MODEL_MPI_  // required in ParticleSystem.h
//#endif

#ifdef _MODEL_MPI_
#define _MODEL_DD_MPI_
#endif


#ifdef _OPENMP
#include <omp.h>
#define EIGEN_DONT_PARALLELIZE // disable Eigen Internal openmp Parallelization
#endif

#include <chrono>
#include <Eigen/Dense>

//#include <model/Network/Network.h>
#include <model/LoopNetwork/LoopNetwork.h>
#include <model/Utilities/TerminalColors.h>
#include <model/IO/EigenDataReader.h>
#include <model/DislocationDynamics/DislocationConsts.h>
#include <model/DislocationDynamics/DislocationNetworkTraits.h>
#include <model/DislocationDynamics/DislocationNetworkComponent.h>
#include <model/DislocationDynamics/DislocationNode.h>
#include <model/DislocationDynamics/DislocationSegment.h>
#include <model/DislocationDynamics/DislocationLoop.h>
#include <model/DislocationDynamics/DislocationSharedObjects.h>
#include <model/DislocationDynamics/GlidePlanes/GlidePlaneObserver.h>
#include <model/DislocationDynamics/DislocationNetworkRemesh.h>
//#include <model/DislocationDynamics/Junctions/DislocationJunctionFormation.h>
//#include <model/DislocationDynamics/CrossSlip/DislocationCrossSlip.h>
#include <model/DislocationDynamics/Materials/Material.h>
#include <model/DislocationDynamics/DislocationNetworkIO.h>
#include <model/DislocationDynamics/NearestNeighbor/DislocationParticle.h>
#include <model/DislocationDynamics/NearestNeighbor/DislocationStress.h>
#include <model/ParticleInteraction/ParticleSystem.h>
#include <model/MPI/MPIcout.h> // defines mode::cout
#include <model/ParticleInteraction/SingleFieldPoint.h>
#include <model/DislocationDynamics/DislocationNodeContraction.h>
#include <model/Threads/EqualIteratorRange.h>
//#include <model/DislocationDynamics/Polycrystals/GrainBoundaryTransmission.h>
//#include <model/DislocationDynamics/Polycrystals/GrainBoundaryDissociation.h>


namespace model
{
    
    template <short unsigned int _dim, short unsigned int corder, typename InterpolationType,
    /*	   */ template <short unsigned int, size_t> class QuadratureRule>
    class DislocationNetwork :
    /* inheritance          */ public LoopNetwork<DislocationNetwork<_dim,corder,InterpolationType,QuadratureRule> >,
    /* inheritance          */ public GlidePlaneObserver<typename TypeTraits<DislocationNetwork<_dim,corder,InterpolationType,QuadratureRule> >::LinkType>,
    /* inheritance          */ public ParticleSystem<DislocationParticle<_dim> >
    {
        
    public:
        
        enum {dim=_dim}; // make dim available outside class
        
        typedef DislocationNetwork<dim,corder,InterpolationType,QuadratureRule> DislocationNetworkType;
        typedef LoopNetwork<DislocationNetworkType> LoopNetworkType;
//        typedef DislocationNetworkType Derived; // define Derived to use NetworkTypedefs.h
//#include <model/Network/NetworkTypedefs.h>
        typedef DislocationNode<dim,corder,InterpolationType,QuadratureRule> NodeType; 		// Define "LinkType" so that NetworkTypedefs.h can be used
        typedef DislocationSegment<dim,corder,InterpolationType,QuadratureRule> LinkType; 		// Define "LinkType" so that NetworkTypedefs.h can be used
        typedef DislocationNetworkComponent<NodeType,LinkType> DislocationNetworkComponentType;
        typedef NetworkComponent<NodeType,LinkType> NetworkComponentType;
//        typedef NetworkComponentObserver<NetworkComponentType> NetworkComponentObserverType;
//        typedef typename NetworkComponentObserverType::NetworkComponentContainerType NetworkComponentContainerType;
        typedef Eigen::Matrix<double,dim,dim>	MatrixDimD;
        typedef Eigen::Matrix<double,dim,1>		VectorDimD;
        typedef GlidePlaneObserver<LinkType> GlidePlaneObserverType;
        typedef DislocationParticle<_dim> DislocationParticleType;
        typedef typename DislocationParticleType::StressField StressField;
        typedef typename DislocationParticleType::DisplacementField DisplacementField;
        typedef ParticleSystem<DislocationParticleType> ParticleSystemType;
        typedef typename ParticleSystemType::SpatialCellType SpatialCellType;
        typedef SpatialCellObserver<DislocationParticleType,_dim> SpatialCellObserverType;
        typedef DislocationSharedObjects<dim> DislocationSharedObjectsType;
        typedef typename DislocationSharedObjectsType::BvpSolverType BvpSolverType;
        typedef typename DislocationSharedObjectsType::BvpSolverType::FiniteElementType FiniteElementType;
        typedef typename FiniteElementType::ElementType ElementType;
        typedef typename LoopNetworkType::IsNodeType IsNodeType;
//        typedef DislocationNodeContraction<DislocationNetworkType> ContractionType;
        
        enum {NdofXnode=NodeType::NdofXnode};
        
#ifdef DislocationNucleationFile
#include DislocationNucleationFile
#endif
        
    private:
        
        bool check_balance;
        short unsigned int use_redistribution;
        bool use_junctions;
        double shearWaveSpeedFraction;
        
        long int runID;
        
//        bool use_crossSlip;
        
        double totalTime;
        
        double dx, dt;
        double vmax;
        
        
        /**********************************************************************/
        void formJunctions()
        {/*! Performs dislocation junction formation if use_junctions==true
          */
            if (use_junctions)
            {
//                DislocationJunctionFormation<DislocationNetworkType>(*this).formJunctions(dx);
            }
        }
        
        /**********************************************************************/
        void make_dt()
        {/*! Computes the time step size \f$dt\f$ for the current simulation step,
          *  based on maximum nodal velocity \f$v_{max}\f$.
          *
          *  The time step is calculated according to:
          *	\f[
          *  dt=
          *  \begin{cases}
          *		\frac{dx}{v_{max}} & v_{max} > fc_s\\
          *      \frac{dx}{fc_s} & v_{max} \le fc_s\\
          *  \end{cases}
          *	\f]
          *  where \f$c_s\f$ is the shear velocity and \f$f=0.1\f$ is a constant.
          */
            model::cout<<"		Computing dt..."<<std::flush;
            const auto t0=std::chrono::system_clock::now();
            
            //			double vmax(0.0);
            vmax=0.0;
            int nVmean=0;
            double vmean=0.0;
            double dt_mean=0.0;
            
            for (const auto& nodeIter : this->nodes())
            {
                if(!nodeIter.second->isBoundaryNode() && !nodeIter.second->isConnectedToBoundaryNodes() && nodeIter.second->confiningPlanes().size()<3)
                {
                    const double vNorm(nodeIter.second->get_V().norm());
                    vmean +=vNorm;
                    nVmean++;
                    if (vNorm>vmax)
                    {
                        vmax=vNorm;
                    }
                }
            }
            vmean/=nVmean;
            
            //double shearWaveSpeedFraction(0.01);
            //short unsigned int shearWaveExp=1;
            if (vmax > Material<Isotropic>::cs*shearWaveSpeedFraction)
            {
                dt=dx/vmax;
                
            }
            else
            {
                //dt=dx/std::pow(shared.material.cs*shearWaveSpeedFraction,shearWaveExp+1)*std::pow(vmax,shearWaveExp);
                //dt=dx/(shared.material.cs*shearWaveSpeedFraction)*std::pow(vmax/(shared.material.cs*shearWaveSpeedFraction),1);
                dt=dx/(Material<Isotropic>::cs*shearWaveSpeedFraction);
            }
            
            if (vmean > Material<Isotropic>::cs*shearWaveSpeedFraction)
            {
                dt_mean=dx/vmean;
            }
            else
            {
                dt_mean=dx/(Material<Isotropic>::cs*shearWaveSpeedFraction);
            }
            
            model::cout<<std::setprecision(3)<<std::scientific<<" vmax="<<vmax;
            model::cout<<std::setprecision(3)<<std::scientific<<" vmax/cs="<<vmax/Material<Isotropic>::cs;
            model::cout<<std::setprecision(3)<<std::scientific<<" dt="<<dt;
            model::cout<<std::setprecision(3)<<std::scientific<<" eta_dt="<<dt/dt_mean;
            model::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
        }
        
//        /**********************************************************************/
//        void crossSlip()
//        {/*! Performs dislocation cross-slip if use_crossSlip==true
//          */
//            if(use_crossSlip)
//            {
//                const auto t0= std::chrono::system_clock::now();
//                model::cout<<"		performing cross-slip ... ("<<std::flush;
//                size_t crossSlipEvents(DislocationCrossSlip<DislocationNetworkType>(*this).crossSlip());
//                model::cout<<crossSlipEvents<<" cross-slip events) ";
//                model::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
//            }
//        }
        
        /**********************************************************************/
        void update_BVP_Solution()
        {
            
            //! Update the BonudaryDislocationNetwork
            if(shared.use_bvp)
            {
                //                if(shared.use_virtualSegments) OLD APPROACH WITH BOUNDARY DISLOCATION NETWORK
                //                {
                //                    const auto t1= std::chrono::system_clock::now();
                //                    model::cout<<"		Updating virtual segments..."<<std::flush;
                //                    shared.bdn.update(*this);
                //                    model::cout<<" ("<<shared.bdn.size()<<" boundary segments)"
                //                    /*       */<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t1)).count()<<" sec]."<<defaultColor<<std::endl;
                //                }
                
                // enter the if statement if use_bvp!=0 and runID is a multiple of use_bvp
                if (!(runID%shared.use_bvp))
                {
                    model::cout<<"		Updating elastic bvp... "<<std::endl;
                    //                shared.bvpSolver.template assembleAndSolve<DislocationNetworkType,4>(*this);
                    shared.bvpSolver.template assembleAndSolve<DislocationNetworkType,37>(*this);
                }
                
            }
            
        }
        
        /**********************************************************************/
        void singleStep()
        {
            //! A simulation step consists of the following:
            model::cout<<blueBoldColor<< "runID="<<runID<<" (of "<<Nsteps<<")"
            /*                    */<< ", time="<<totalTime
            /*                    */<< ": nodeOrder="<<this->nodes().size()
            /*                    */<< ", linkOrder="<<this->links().size()
            /*                    */<< ", components="<<this->components().size()
            /*                    */<< defaultColor<<std::endl;
            
            //! 1- Check that all nodes are balanced
//            checkBalance();
            
            //! 2 - Update quadrature points
            updateQuadraturePoints();
            
            //! 3- Calculate BVP correction
            update_BVP_Solution();
            
#ifdef DislocationNucleationFile
            if(shared.use_bvp && !(runID%shared.use_bvp))
            {
                nucleateDislocations(); // needs to be called before updateQuadraturePoints()
                updateQuadraturePoints();
            }
#endif
            
            //! 4- Solve the equation of motion
            assembleAndSolve();
            
            //! 5- Compute time step dt (based on max nodal velocity) and increment totalTime
            make_dt();
            
            if(DislocationNetworkIO<DislocationNetworkType>::outputElasticEnergy)
            {
                typedef typename DislocationParticleType::ElasticEnergy ElasticEnergy;
                this->template computeNeighborField<ElasticEnergy>();
            }
            
            //! 6- Output the current configuration before changing it
            output(runID);
            
            
            //! 7- Moves DislocationNodes(s) to their new configuration using stored velocity and dt
            move(dt);

            //! 8- Update accumulated quantities (totalTime and plasticDistortion)
            totalTime+=dt;
            const MatrixDimD pdr(plasticDistortionRate()); // move may limit velocity, so compute pdr after move
            plasticDistortion += pdr*dt;
            
            
            //! 9- Contract segments of zero-length
//            DislocationNetworkRemesh<DislocationNetworkType>(*this).contract0chordSegments();
            
            //! 10- Cross Slip (needs upated PK force)
//            DislocationCrossSlip<DislocationNetworkType>(*this);

            
//            GrainBoundaryTransmission<DislocationNetworkType>(*this).transmit();
//            
//            GrainBoundaryDissociation<DislocationNetworkType>(*this).dissociate();

//            shared.poly.reConnectGrainBoundarySegments(*this); // this makes stressGauss invalid, so must follw other GB operations

            
            //! 11- detect loops that shrink to zero and expand as inverted loops
//            DislocationNetworkRemesh<DislocationNetworkType>(*this).loopInversion(dt);
            
            //! 12- Form Junctions
            formJunctions();
            
//            // Remesh may contract juncitons to zero lenght. Remove those juncitons:
//            DislocationJunctionFormation<DislocationNetworkType>(*this).breakZeroLengthJunctions();
            
            //! 13- Node redistribution
//            remesh();
            
            //! 9- Contract segments of zero-length
//            DislocationNetworkRemesh<DislocationNetworkType>(*this).contract0chordSegments();
            
            //! 14- If BVP solver is not used, remove DislocationSegment(s) that exited the boundary
//            removeBoundarySegments();

//            removeSmallComponents(3.0*dx,4);
            
            make_bndNormals();
            
            //! 16 - Increment runID counter
            ++runID;     // increment the runID counter
        }
        
        /**********************************************************************/
        void make_bndNormals()
        {
            for(auto& node : this->nodes())
            {
                node.second->make_bndNormal();
            }
        }
        
//        /**********************************************************************/
//        void removeBoundarySegments()
//        {/*! Removes DislocationSegment(s) on the mesh boundary
//          */
//            if (shared.use_boundary && !shared.use_bvp)
//            {
//                const auto t0= std::chrono::system_clock::now();
//                model::cout<<"		Removing DislocationSegments on mesh boundary... ";
//                typedef bool (LinkType::*link_member_function_pointer_type)(void) const;
//                link_member_function_pointer_type boundarySegment_Lmfp;
//                boundarySegment_Lmfp=&LinkType::is_boundarySegment;
//                this->template disconnect_if<1>(boundarySegment_Lmfp);
//                model::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
//            }
//        }
        
//        /**********************************************************************/
//        /**********************************************************************/
//        void removeSmallComponents(const double& smallcritvalue,
//                                   const size_t& maxNodeSize)
//        {
//            const auto t0= std::chrono::system_clock::now();
//            model::cout<<"        removeSmallComponents "<<std::flush;
//            size_t removed=0;
//            std::deque<int> nodesToBeRemoved;
//            for (typename NetworkComponentContainerType::iterator snIter=this->ABbegin(); snIter!=this->ABend();++snIter)
//            {
//                if (DislocationNetworkComponentType(*snIter->second).isSmall(smallcritvalue,maxNodeSize))
//                {
//                    const auto nodes=DislocationNetworkComponentType(*snIter->second).networkComponent().nodes();
//                    for(const auto& node : nodes)
//                    {
//                        if(this->node(node.first).first)
//                        {
//                            nodesToBeRemoved.push_back(node.first);
//                        }
//                    }
//                }
//            }
//            for(int k:nodesToBeRemoved)
//            {
//                this->template removeVertex<false>(k);
//                removed++;
//                //
//                //                model::cout<<"  node "<<K<<" is deleted!!!!"<<std::endl;
//            }
//            model::cout<<std::setprecision(3)<<std::scientific<<removed<<" removed) "<<magentaColor<<"["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<defaultColor<<std::endl;
//        }
        
//            for(const auto& segment : this->links())
//            {
//            if(segment.second.is_boundarySegment())
//            {
//                const LatticeVector<dim>& b(segment.second.flow);
//                const VectorDimD P0(some position here);
//                const VectorDimD P1(some position here);
//                const size_t nodeID0(DN.insertVertex(P0).first->first);
//                const size_t nodeID1(DN.insertVertex(P1).first->first);
//                this->connect(nodeID0,nodeID1,b);
//            }
//            }
        
    public:
        
        //! The object conataing shared data
        DislocationSharedObjectsType shared;
        
        //! The number of simulation steps taken by the next call to runSteps()
        int Nsteps;
        
//        //! The simulation time run run foby the next call to runTime()
//        double timeWindow;
        
        //! The accumulated plastic distortion
        MatrixDimD plasticDistortion;
        
        /**********************************************************************/
        DislocationNetwork(int& argc, char* argv[]) :
        /* init list  */ check_balance(true),
        /* init list  */ use_redistribution(0),
        /* init list  */ use_junctions(false),
        //		/* init list  */ useImplicitTimeIntegration(false),
        /* init list  */ shearWaveSpeedFraction(0.0001),
        /* init list  */ runID(0),
        /* init list  */ totalTime(0.0),
        /* init list  */ dx(0.0),
        /* init list  */ dt(0.0),
        /* init list  */ vmax(0.0),
        /* init list  */ Nsteps(0),
//        /* init list  */ timeWindow(0.0),
        /* init list  */ plasticDistortion(MatrixDimD::Zero())
        {
            ParticleSystemType::initMPI(argc,argv);
            read("./","DDinput.txt");
        }
        
        /**********************************************************************/
        ~DislocationNetwork()
        {
//            for(auto& linkIter : this->links())
//            {
//                linkIter.second->quadratureParticleContainer.clear();
//            }
            
            this->clearParticles();
        }
        
        /**********************************************************************/
        size_t contractWithConstraintCheck(const IsNodeType& Ni,
                                           const IsNodeType& Nj)
        {
            return ContractionType(*this).contractWithConstraintCheck(Ni,Nj);
        }
        
        
        /**********************************************************************/
        void remesh()
        {
            if (use_redistribution)
            {
                if(!(runID%use_redistribution))
                {
                    DislocationNetworkRemesh<DislocationNetworkType>(*this).remesh();
                }
            }
        }
        
        /**********************************************************************/
        const double& get_dt() const
        {/*! The current simulation time step in dimensionless units
          */
            return dt;
        }
        
        /**********************************************************************/
        const double& get_totalTime() const
        {/*! The elapsed simulation time step in dimensionless units
          */
            return totalTime;
        }
        
        /**********************************************************************/
        void read(const std::string& inputDirectoryName_in, std::string inputFileName)
        { // TO DO: move this to DislocationNetworkIO.h
            
            std::ostringstream fullName;
            fullName<<inputDirectoryName_in<<inputFileName;
            
            model::cout<<greenColor<<"Reading "<<fullName.str()<<"..."<<defaultColor<<std::endl;
            
            
            // Create a file-reader object
            EigenDataReader EDR;
            
            // IO
            EDR.readScalarInFile(fullName.str(),"outputFrequency",DislocationNetworkIO<DislocationNetworkType>::outputFrequency);
            EDR.readScalarInFile(fullName.str(),"outputBinary",DislocationNetworkIO<DislocationNetworkType>::outputBinary);
            EDR.readScalarInFile(fullName.str(),"outputGlidePlanes",DislocationNetworkIO<DislocationNetworkType>::outputGlidePlanes);
            EDR.readScalarInFile(fullName.str(),"outputSpatialCells",DislocationNetworkIO<DislocationNetworkType>::outputSpatialCells);
            EDR.readScalarInFile(fullName.str(),"outputPKforce",DislocationNetworkIO<DislocationNetworkType>::outputPKforce);
            EDR.readScalarInFile(fullName.str(),"outputMeshDisplacement",DislocationNetworkIO<DislocationNetworkType>::outputMeshDisplacement);
            EDR.readScalarInFile(fullName.str(),"outputElasticEnergy",DislocationNetworkIO<DislocationNetworkType>::outputElasticEnergy);
            
            
            EDR.readScalarInFile(fullName.str(),"outputPlasticDistortion",DislocationNetworkIO<DislocationNetworkType>::outputPlasticDistortion);
            if(DislocationNetworkIO<DislocationNetworkType>::outputPlasticDistortion)
            {
                DislocationNetworkIO<DislocationNetworkType>::_userOutputColumn+=9;
            }
            EDR.readScalarInFile(fullName.str(),"outputPlasticDistortionRate",DislocationNetworkIO<DislocationNetworkType>::outputPlasticDistortionRate);
            if(DislocationNetworkIO<DislocationNetworkType>::outputPlasticDistortionRate)
            {
                DislocationNetworkIO<DislocationNetworkType>::_userOutputColumn+=9;
            }
            
            EDR.readScalarInFile(fullName.str(),"outputDislocationLength",DislocationNetworkIO<DislocationNetworkType>::outputDislocationLength);
            if(DislocationNetworkIO<DislocationNetworkType>::outputDislocationLength)
            {
                DislocationNetworkIO<DislocationNetworkType>::_userOutputColumn+=3;
            }
            
            EDR.readScalarInFile(fullName.str(),"outputQuadratureParticles",DislocationNetworkIO<DislocationNetworkType>::outputQuadratureParticles);
            
            // Parametrization exponent
            EDR.readScalarInFile(fullName.str(),"parametrizationExponent",LinkType::alpha);
            assert((LinkType::alpha)>=0.0 && "parametrizationExponent MUST BE >= 0.0");
            assert((LinkType::alpha)<=1.0 && "parametrizationExponent MUST BE <= 1.0");
            
            // Temperature. Make sure you initialize before calling Material<Isotropic>::select()
            EDR.readScalarInFile(fullName.str(),"temperature",Material<Isotropic>::T); // temperature
            
            // Material and crystal orientation
            unsigned int materialZ;
            EDR.readScalarInFile(fullName.str(),"material",materialZ); // material by atomic number Z
            Material<Isotropic>::select(materialZ);
//            MatrixDimD C2Gtemp;
//            EDR.readMatrixInFile(fullName.str(),"C2G",C2Gtemp); // crystal-to-global orientation
//            Material<Isotropic>::rotateCrystal(C2Gtemp);
            
            // quadPerLength
            EDR.readScalarInFile(fullName.str(),"quadPerLength",LinkType::quadPerLength); // quadPerLength
            
            // core size
            EDR.readScalarInFile(fullName.str(),"coreSize",StressField::a); // core-width
            assert((StressField::a)>0.0 && "coreSize MUST BE > 0.");
            StressField::a2=StressField::a*StressField::a;
            
            // multipole expansion
            double cellSize(0.0);
            EDR.readScalarInFile(fullName.str(),"dislocationCellSize",cellSize);
            SpatialCellObserverType::setCellSize(cellSize);
            EDR.readScalarInFile(fullName.str(),"use_DisplacementMultipole",DislocationDisplacement<dim>::use_multipole);
            EDR.readScalarInFile(fullName.str(),"use_StressMultipole",DislocationStress<dim>::use_multipole);
            EDR.readScalarInFile(fullName.str(),"use_EnergyMultipole",DislocationEnergy<dim>::use_multipole);
            
            // Eternal Stress
            EDR.readMatrixInFile(fullName.str(),"externalStress",shared.externalStress);
            
            // Restart
            EDR.readScalarInFile(fullName.str(),"startAtTimeStep",runID);
            VertexReader<'F',201,double> vReader;
            Eigen::Matrix<double,1,200> temp(Eigen::Matrix<double,1,200>::Zero());
            
            
            if (vReader.isGood(0,true))
            {
                
                vReader.read(0,true);
                
                if(runID<0)
                {
                    if(vReader.size())
                    {
                        runID=vReader.rbegin()->first;
                        temp=vReader.rbegin()->second;
                    }
                    else
                    {
                        runID=0;
                    }
                }
                else
                {
                    const auto iter=vReader.find(runID);
                    if(iter!=vReader.end())
                    {// runID has been found
                        temp=iter->second;
                    }
                    else
                    {
                        assert(0 && "runID NOT FOUND IN F/F_0.txt");
                    }
                }
            }
            else
            {
                model::cout<<"could not read runID from F/F_0.txt"<<std::endl;
                runID=0;
            }
            
            size_t curCol=0;
            totalTime=temp(curCol);
            curCol+=2;
            
            if (DislocationNetworkIO<DislocationNetworkType>::outputPlasticDistortion)
            {
                std::cout<<"reading PD"<<std::endl;
                
                for(int r=0;r<3;++r)
                {
                    for(int c=0;c<3;++c)
                    {
                        plasticDistortion(r,c)=temp(curCol);
                        curCol+=1;
                    }
                }
            }
            
            model::cout<<"starting at time step "<<runID<<std::endl;
            model::cout<<"totalTime= "<<totalTime<<std::endl;
            model::cout<<"plasticDistortion=\n "<<plasticDistortion<<std::endl;
            
            // time-stepping
            dt=0.0;
            EDR.readScalarInFile(fullName.str(),"dx",dx);
            assert(dx>0.0);
            //            EDR.readScalarInFile(fullName.str(),"shearWaveSpeedFraction",shearWaveSpeedFraction);
            //            assert(shearWaveSpeedFraction>=0.0);
            EDR.readScalarInFile(fullName.str(),"use_velocityFilter",NodeType::use_velocityFilter);
            EDR.readScalarInFile(fullName.str(),"velocityReductionFactor",NodeType::velocityReductionFactor);
            assert(NodeType::velocityReductionFactor>0.0 && NodeType::velocityReductionFactor<=1.0);
            EDR.readScalarInFile(fullName.str(),"bndDistance",NodeType::bndDistance);
            //            EDR.readScalarInFile(fullName.str(),"useImplicitTimeIntegration",useImplicitTimeIntegration);
            EDR.readScalarInFile(fullName.str(),"use_directSolver_DD",DislocationNetworkComponentType::use_directSolver);
            
            
            EDR.readScalarInFile(fullName.str(),"Nsteps",Nsteps);
            assert(Nsteps>=0 && "Nsteps MUST BE >= 0");
            
//            EDR.readScalarInFile(fullName.str(),"timeWindow",timeWindow);
//            assert(timeWindow>=0.0 && "timeWindow MUST BE >= 0");
            
            // Check balance
            EDR.readScalarInFile(fullName.str(),"check_balance",check_balance);
            
            // JUNCTION FORMATION
            EDR.readScalarInFile(fullName.str(),"use_junctions",use_junctions);
//            EDR.readScalarInFile(fullName.str(),"collisionTol",DislocationJunctionFormation<DislocationNetworkType>::collisionTol);
            
            // VERTEX REDISTRIBUTION
            EDR.readScalarInFile(fullName.str(),"use_redistribution",use_redistribution);
            EDR.readScalarInFile(fullName.str(),"Lmin",DislocationNetworkRemesh<DislocationNetworkType>::Lmin);
            assert(DislocationNetworkRemesh<DislocationNetworkType>::Lmin>=0.0);
            assert(DislocationNetworkRemesh<DislocationNetworkType>::Lmin>=2.0*dx && "YOU MUST CHOOSE Lmin>2*dx.");
            EDR.readScalarInFile(fullName.str(),"Lmax",DislocationNetworkRemesh<DislocationNetworkType>::Lmax);
            assert(DislocationNetworkRemesh<DislocationNetworkType>::Lmax>DislocationNetworkRemesh<DislocationNetworkType>::Lmin);
            
            // Cross-Slip
//            EDR.readScalarInFile(fullName.str(),"use_crossSlip",DislocationCrossSlip<DislocationNetworkType>::use_crossSlip);
//            if(DislocationCrossSlip<DislocationNetworkType>::use_crossSlip)
//            {
//                EDR.readScalarInFile(fullName.str(),"crossSlipDeg",DislocationCrossSlip<DislocationNetworkType>::crossSlipDeg);
//                assert(DislocationCrossSlip<DislocationNetworkType>::crossSlipDeg>=0.0 && DislocationCrossSlip<DislocationNetworkType>::crossSlipDeg <= 90.0 && "YOU MUST CHOOSE 0.0<= crossSlipDeg <= 90.0");
//                //                EDR.readScalarInFile(fullName.str(),"crossSlipLength",DislocationCrossSlip<DislocationNetworkType>::crossSlipLength);
//                //                assert(DislocationCrossSlip<DislocationNetworkType>::crossSlipLength>=DislocationNetworkRemesh<DislocationNetworkType>::Lmin && "YOU MUST CHOOSE crossSlipLength>=Lmin.");
//            }
            
            // Mesh and BVP
            EDR.readScalarInFile(fullName.str(),"use_boundary",shared.use_boundary);
            if (shared.use_boundary)
            {
//                EDR.readScalarInFile(fullName.str(),"use_meshRegions",shared.use_meshRegions);
                
                int meshID(0);
                EDR.readScalarInFile(fullName.str(),"meshID",meshID);
                shared.mesh.readMesh(meshID);
                assert(shared.mesh.simplices().size() && "MESH IS EMPTY.");
                
                // Initialize Polycrystal
                shared.poly.init(fullName.str());
                
                EDR.readScalarInFile(fullName.str(),"use_virtualSegments",shared.use_virtualSegments);
                if(shared.use_virtualSegments)
                {
                    EDR.readScalarInFile(fullName.str(),"virtualSegmentDistance",LinkType::virtualSegmentDistance);
                }
                
                EDR.readScalarInFile(fullName.str(),"use_bvp",shared.use_bvp);
                if(shared.use_bvp)
                {
                    EDR.readScalarInFile(fullName.str(),"use_directSolver_FEM",shared.bvpSolver.use_directSolver);
                    EDR.readScalarInFile(fullName.str(),"solverTolerance",shared.bvpSolver.tolerance);
                    shared.bvpSolver.init(*this);
                }
            }
            else{ // no boundary is used, DislocationNetwork is in inifinite medium
                shared.use_bvp=0;	// never comupute boundary correction
            }
            
            // Grain Boundary flags
//            EDR.readScalarInFile(fullName.str(),"use_GBdissociation",GrainBoundaryDissociation<DislocationNetworkType>::use_GBdissociation);
//            EDR.readScalarInFile(fullName.str(),"use_GBtransmission",GrainBoundaryTransmission<DislocationNetworkType>::use_GBtransmission);
//            EDR.readScalarInFile(fullName.str(),"use_GBdislocations",GrainBoundary<dim>::use_GBdislocations);

            // Read Vertex and Edge information
            DislocationNetworkIO<DislocationNetworkType>::readVertices(*this,runID); // this requires mesh to be up-to-date
            DislocationNetworkIO<DislocationNetworkType>::readEdges(*this,runID);    // this requires mesh to be up-to-date
            
            //#ifdef DislocationNucleationFile
            //            EDR.readScalarInFile(fullName.str(),"nucleationFreq",nucleationFreq);
            //#endif
            
#ifdef _MODEL_MPI_
            // Avoid that a processor starts writing before other are reading
            MPI_Barrier(MPI_COMM_WORLD);
#endif
            
            // Initializing configuration
            move(0.0);	// initial configuration
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
            model::cout<<"		Computing dislocation-dislocation interactions ("<<nThreads<<" threads)..."<<std::flush;
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
                model::cout<<"(PardisoLDLT)..."<<std::flush;
                for (const auto& networkComponent : this->components())
                {
                    DislocationNetworkComponentType(*networkComponent.second).directSolve();
                }
#else
                model::cout<<"(SimplicialLDLT)..."<<std::flush;
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
                model::cout<<"(MINRES)..."<<std::flush;
#ifdef _OPENMP // SimplicialLDLT is not multi-threaded. So parallelize loop over NetworkComponents.
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
        void output(const size_t& fileID) const
        {/*! Outputs DislocationNetwork data
          */
            if (!(runID%DislocationNetworkIO<DislocationNetworkType>::outputFrequency))
            {
                const auto t0=std::chrono::system_clock::now();
#ifdef _MODEL_DD_MPI_
                if(ModelMPIbase::mpiRank()==0)
                {
                    DislocationNetworkIO<DislocationNetworkType>::output(*this,fileID);
                }
#else
                DislocationNetworkIO<DislocationNetworkType>::output(*this,fileID);
#endif
                model::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
            }
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
                nodeIter.second->move(dt_in);
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
                
//        /**********************************************************************/
//        void checkBalance() const
//        {/*! Checks that each DislocationNode is balanced, and asserts otherwise
//          * Exceptions are nodes with only one neighbors (FR source)
//          * and nodes on the boundary.
//          */
//            if(check_balance)
//            {
//                model::cout<<"		Checking node balance..."<<std::flush;
//                const auto t0= std::chrono::system_clock::now();
//                for (const auto& nodeIter : this->nodes())
//                {
//                    if (nodeIter.second->neighbors().size()>2)
//                    {
//                        const bool nodeIsBalanced(nodeIter.second->is_balanced());
//                        if (!nodeIsBalanced && nodeIter.second->meshLocation()==insideMesh)
//                        {
//                            model::cout<<"Node "<<nodeIter.second->sID<<" is not balanced:"<<std::endl;
////                            model::cout<<"    outflow="<<nodeIter.second->outFlow().transpose()<<std::endl;
////                            model::cout<<"     inflow="<<nodeIter.second->inFlow().transpose()<<std::endl;
//                            assert(0 && "NODE IS NOT BALANCED");
//                        }
//                    }
//                }
//                model::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
//            }
//        }
        
        
        /**********************************************************************/
        MatrixDimD plasticDistortionRate() const
        {
            MatrixDimD temp(MatrixDimD::Zero());
            for (const auto& linkIter : this->links())
            {
                temp+= linkIter.second->plasticDistortionRate();
            }
            return temp;
        }
        
        /**********************************************************************/
        MatrixDimD plasticStrainRate() const
        {/*!\returns the plastic distortion rate tensor generated by this
          * DislocaitonNetwork.
          */
            const MatrixDimD temp(plasticDistortionRate());
            return (temp+temp.transpose())*0.5;
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
                if(linkIter.second->is_boundarySegment())
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
        const double& vMax() const
        {/*! The maximum vertex velocity.
          */
            return vmax;
        }
        
        /**********************************************************************/
        MatrixDimD stress(const VectorDimD& P) const
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
        std::pair<bool,const Simplex<dim,dim>*> pointIsInsideMesh(const VectorDimD& P0, const Simplex<dim,dim>* const guess) const
        {/*!\param[in] P0 position vector
          * \param[in] guess pointer of the Simplex where the search starts
          * \returns true if P0 is inside the mesh
          */
            std::pair<bool,const Simplex<dim,dim>*> temp(true,NULL);
            if (shared.use_boundary)
            {
                temp=shared.mesh.searchWithGuess(P0,guess);
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
        
//            /**********************************************************************/
//            template <typename FEMfunctionType,typename ...Args>
//            void femIntegralLoop(FEMfunctionType& femObj,const Args&... args) const
//            {
//                for (const auto& segment : this->links())
//                {
//                    const Simplex<dim,dim>* guess=segment.second.source->includingSimplex();
//                    for (size_t k=0;k<segment.second.qOrder;++k)
//                    {
//                        const auto& x=segment.second.rgauss.col(k);
//                        
//                        std::pair<bool,const ElementType*> elem=DN.bvpSolver.fe2().searcWithGuess(x,guess);
//                        
//                        if(elem.first)
//                        {
//                            guess=&elem.second->simplex;
//                            
//                            const Eigen::Matrix<int,1,dim+1> bary=elem.second->simplex.pos2bary(x);
//                            femObj(segment.second,k,elem.second,bary,args...);
//                        }
//                    }
//                }
//            }
    };
    
} // namespace model
#endif

