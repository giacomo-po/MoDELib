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
// Motion of BoundarySubnetworks

// TO DO
// -5 use modelMacros instead of assert
// -4 CHANGE definition of modelMacros in debug mode BECAUSE assert(x) CAN BE REMOVED BY NDEBUG!!! MUST USE mode::assert_fail
// -3 MAKE isIsolated and isBalanced data members of NetworkNode that are modified by addToNeighbors and removeFromNeigbors
// implement IndependentPath
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

#ifdef _MODEL_DD_MPI_
#define _MODEL_MPI_  // required in ParticleSystem.h
#endif

#ifdef _OPENMP
#include <omp.h>
//#define EIGEN_DONT_PARALLELIZE // disable Eigen Internal openmp Parallelization
#endif


#include <Eigen/Dense>

#include <model/Network/Network.h>
#include <model/Utilities/TerminalColors.h>
#include <model/Utilities/EigenDataReader.h>
#include <model/DislocationDynamics/DislocationConsts.h>
#include <model/DislocationDynamics/DislocationNetworkTraits.h>
#include <model/DislocationDynamics/DislocationNetworkComponent.h>
#include <model/DislocationDynamics/DislocationNode.h>
#include <model/DislocationDynamics/DislocationSegment.h>
#include <model/DislocationDynamics/DislocationSharedObjects.h>
#include <model/DislocationDynamics/GlidePlanes/GlidePlaneObserver.h>
#include <model/DislocationDynamics/Remeshing/DislocationNetworkRemesh.h>
#include <model/DislocationDynamics/Junctions/DislocationJunctionFormation.h>
#include <model/DislocationDynamics/CrossSlip/DislocationCrossSlip.h>
#include <model/DislocationDynamics/Materials/Material.h>
#include <model/DislocationDynamics/IO/DislocationNetworkIO.h>
#include <model/DislocationDynamics/NearestNeighbor/DislocationParticle.h>
#include <model/DislocationDynamics/NearestNeighbor/DislocationStress.h>
#include <model/ParticleInteraction/ParticleSystem.h>
#include <model/MPI/MPIcout.h> // defines mode::cout



namespace model
{
 
	template <short unsigned int _dim, short unsigned int corder, typename InterpolationType,
	/*	   */ short unsigned int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule>
	class DislocationNetwork :
    /* inheritance          */ public Network<DislocationNetwork<_dim,corder,InterpolationType,qOrder,QuadratureRule> >,
	/* inheritance          */ public GlidePlaneObserver<typename TypeTraits<DislocationNetwork<_dim,corder,InterpolationType,qOrder,QuadratureRule> >::LinkType>,
    /* inheritance          */ public ParticleSystem<DislocationParticle<_dim> >
    {
		
    public:
        
        enum {dim=_dim}; // make dim available outside class
        
		typedef DislocationNetwork<dim,corder,InterpolationType,qOrder,QuadratureRule> DislocationNetworkType;
		typedef DislocationNetworkType Derived; // define Derived to use NetworkTypedefs.h
#include <model/Network/NetworkTypedefs.h>
        typedef DislocationNetworkComponent<NodeType,LinkType> DislocationNetworkComponentType;
		typedef Eigen::Matrix<double,dim,dim>	MatrixDimD;
		typedef Eigen::Matrix<double,dim,1>		VectorDimD;
		typedef GlidePlaneObserver<LinkType> GlidePlaneObserverType;
        typedef DislocationParticle<_dim> DislocationParticleType;
        typedef typename DislocationParticleType::StressField StressField;
        typedef ParticleSystem<DislocationParticleType> ParticleSystemType;
        typedef SpatialCellObserver<DislocationParticleType,_dim> SpatialCellObserverType;
		enum {NdofXnode=NodeType::NdofXnode};
		      
//#ifdef UpdateBoundaryConditionsFile
//#include UpdateBoundaryConditionsFile
//#endif
        
#ifdef DislocationNucleationFile
#include DislocationNucleationFile
#endif
        
	private:
        
		short unsigned int use_redistribution;
		bool use_junctions;
        bool useImplicitTimeIntegration;
        double equilibriumVelocity;
        
		unsigned int runID;
        
		bool use_crossSlip;
		double crossSlipDeg;
		double crossSlipLength;
		
		double totalTime;
		
		double dx, dt;
        double vmax;

        /**********************************************************************/
		void formJunctions()
        {/*! Performs dislocation junction formation if use_junctions==true
          */
			if (use_junctions)
            {
				double t0=clock();
				model::cout<<"		Forming Junctions: found ("<<std::flush;
                const double avoidNodeIntersection(0.05);
				DislocationJunctionFormation<DislocationNetworkType>(*this).formJunctions(dx,avoidNodeIntersection);
				model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]."<<defaultColor<<std::endl;
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
			double t0=clock();
			
            //			double vmax(0.0);
            vmax=0.0;
			
			for (typename NetworkNodeContainerType::iterator nodeIter=this->nodeBegin();nodeIter!=this->nodeEnd();++nodeIter)
            {
				const double vNorm(nodeIter->second->get_V().norm());
				if (vNorm>vmax)
                {
					vmax=vNorm;
				}
			}
			
            //double equilibriumVelocity(0.01);
			//short unsigned int shearWaveExp=1;
			if (vmax > Material<Isotropic>::cs*equilibriumVelocity)
            {
				dt=dx/vmax;
			}
			else
            {
                //dt=dx/std::pow(shared.material.cs*equilibriumVelocity,shearWaveExp+1)*std::pow(vmax,shearWaveExp);
				//dt=dx/(shared.material.cs*equilibriumVelocity)*std::pow(vmax/(shared.material.cs*equilibriumVelocity),1);
				dt=dx/(Material<Isotropic>::cs*equilibriumVelocity);
			}
			model::cout<<std::setprecision(3)<<std::scientific<<" vmax="<<vmax;
			model::cout<<std::setprecision(3)<<std::scientific<<" dt="<<dt;
			model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]."<<defaultColor<<std::endl;
		}
        
        /**********************************************************************/
		void crossSlip()
        {/*! Performs dislocation cross-slip if use_crossSlip==true
          */
			if(use_crossSlip)
            {
				double t0=clock();
				model::cout<<"		Performing Cross Slip ... "<<std::flush;
				size_t crossSlipEvents(DislocationCrossSlip<DislocationNetworkType>(*this).crossSlip(crossSlipDeg,crossSlipLength));
				model::cout<<crossSlipEvents<<" cross slip events found ";
				model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]."<<defaultColor<<std::endl;
			}
		}
        
        /**********************************************************************/
        void update_BVP_Solution()
        {
            // enter the if statement if use_bvp!=0 and runID is a multiple of use_bvp
            if (shared.use_bvp && !(runID%shared.use_bvp))
            {
                    double t0=clock();
                    model::cout<<"		Updating bvp stress ... "<<std::flush;
                    shared.bvpSolver.template assembleAndSolve<DislocationNetworkType,4>(*this);
                    std::cout<<"FINISH HERE";
                    model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]."<<defaultColor<<std::endl;
			}
        }
        
        /**********************************************************************/
		void singleStep()
        {
			//! A simulation step consists of the following:
			model::cout<<blueBoldColor<< "runID="<<runID
            /*                    */<< ", time="<<totalTime
            /*                    */<< ": nodeOrder="<<this->nodeOrder()
            /*                    */<< ", linkOrder="<<this->linkOrder()
            /*                    */<< ", components="<<this->Naddresses()
			/*                    */<< defaultColor<<std::endl;
            
#ifdef DislocationNucleationFile
            if(shared.use_bvp && !(runID%shared.use_bvp))
            {
                nucleateDislocations(); // needs to be called before updateQuadraturePoints()
                removeBoundarySegments();
            }
#endif
            
			//! 1- Check that all nodes are balanced
			checkBalance();
			
			//! 1 - Update quadrature points
			updateQuadraturePoints();
            
			//! 2- Calculate BVP correction
            update_BVP_Solution();
            
			//! 3- Solve the equation of motion
			assembleAndSolve();
            
			//! 4- Compute time step dt (based on max nodal velocity) and increment totalTime
			make_dt();
			totalTime+=dt;
            
            const MatrixDimD pdr(plasticDistortionRate());
            plasticDistortion += pdr*dt;
            
            if(DislocationNetworkIO<DislocationNetworkType>::outputElasticEnergy)
            {
                typedef typename DislocationParticleType::ElasticEnergy ElasticEnergy;
                this->template computeNeighborField<ElasticEnergy>();
            }
            
            //! 11- Output the current configuration before changing it
            output();
            
			//! 5- Moves DislocationNodes(s) to their new configuration using stored velocity and dt
			move(dt,0.0);
            
            //bool useImplicitTimeIntegration(false); // THIS COULD BE DECIDED ON THE FLY BASED ON DISTRIBUTION OF VELOCITIES
            if(useImplicitTimeIntegration)
            {
                updateQuadraturePoints();
                assembleAndSolve(); // this also stores the new velocity in each node
                
                for (typename NetworkNodeContainerType::iterator nodeIter=this->nodeBegin();nodeIter!=this->nodeEnd();++nodeIter)
                {
					nodeIter->second->implicitStep(); // average velocities
				}
                
                const double dt_old(dt); // store current dt
                make_dt();      // compute dt again with average velocity
                totalTime += dt-dt_old; // correct accumulated totalTime
                plasticDistortion += pdr*(dt-dt_old); // // correct accumulated plasticDistortion
                move(dt,dt_old); // move again (internally this subtracts DislocationNode::vOld*dt_old)
            }
            
			DislocationNetworkRemesh<DislocationNetworkType>(*this).contract0chordSegments();
			
			//! 6- Moves DislocationNodes(s) to their new configuration using stored velocity and dt
            DislocationNetworkRemesh<DislocationNetworkType>(*this).loopInversion(dt);
			
			//! 7- If soft boundaries are used, remove DislocationSegment(s) that exited the boundary
			removeBoundarySegments();
			
			//! 8- Form Junctions
			formJunctions();
			
			//! 9- Node redistribution
			remesh();
			removeBoundarySegments();
            
			//! 10- Cross Slip
			crossSlip(); // do crossSlip after remesh so that cross-slip points are not removed
            
            updateQuadraturePoints(); // necessary if quadrature data are necessary in main
            
			//! 12 - Increment runID counter
			++runID;     // increment the runID counter
		}
		
		/**********************************************************************/
		void removeBoundarySegments()
        {/*! Removes DislocationSegment(s) on the mesh boundary
          */
			if (shared.boundary_type==softBoundary)
            {
				double t0=clock();
				model::cout<<"		Removing DislocationSegments outside mesh boundary... ";
				typedef bool (LinkType::*link_member_function_pointer_type)(void) const;
				link_member_function_pointer_type boundarySegment_Lmfp;
				boundarySegment_Lmfp=&LinkType::is_boundarySegment;
				this->template disconnect_if<1>(boundarySegment_Lmfp);
				model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]."<<defaultColor<<std::endl;
			}
		}
        
        /**********************************************************************/
		void segmentMeshCollision()
        {
			typedef std::map<double,VectorDimD> MultiIntersectionType; // keeps intersections sorted for increasing parameter
			typedef std::map<LinkIDType,MultiIntersectionType> MultiIntersectionContainerType; // Container of MultiIntersectionType for different links
			MultiIntersectionContainerType multiIntersectionContainer;
			for (typename NetworkLinkContainerType::iterator linkIter=this->linkBegin();linkIter!=this->linkEnd();++linkIter)
            {
				multiIntersectionContainer[linkIter->second->nodeIDPair()]=linkIter->second->boundaryCollision();
			}
			for (typename MultiIntersectionContainerType::const_iterator iter=multiIntersectionContainer.begin();iter!=multiIntersectionContainer.end();++iter)
            {
				this->multiExpand(iter->first.first,iter->first.second,iter->second);
			}
		}
        
	public:
        //		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
		
        typedef DislocationSharedObjects<LinkType> DislocationSharedObjectsType;
        
        typedef typename DislocationSharedObjectsType::BvpSolverType BvpSolverType;
        
		DislocationSharedObjectsType shared;
		
		//! The number of simulation steps taken by the next call to runByStep()
		int Nsteps;
		double timeWindow;
        
        MatrixDimD plasticDistortion;
        
		/**********************************************************************/
        DislocationNetwork(int& argc, char* argv[]) :
        /* init list  */ use_redistribution(0),
		/* init list  */ use_junctions(false),
		/* init list  */ useImplicitTimeIntegration(false),
        /* init list  */ equilibriumVelocity(0.01),
		/* init list  */ runID(0),
		/* init list  */ use_crossSlip(false),
		/* init list  */ crossSlipDeg(45.0),
		/* init list  */ crossSlipLength(100.0),
		/* init list  */ totalTime(0.0),
		/* init list  */ dx(0.0),
        /* init list  */ dt(0.0),
        /* init list  */ vmax(0.0),
        /* init list  */ Nsteps(0),
		/* init list  */ timeWindow(0.0),
        /* init list  */ plasticDistortion(MatrixDimD::Zero())
        {
            ParticleSystemType::initMPI(argc,argv);
            read("./","DDinput.txt");
        }
        
		/**********************************************************************/
        ~DislocationNetwork()
        {
            for (typename NetworkLinkContainerType::iterator linkIter =this->linkBegin();
                 /*                                       */ linkIter!=this->linkEnd();
                 /*                                       */ linkIter++)
            {
				linkIter->second->quadratureParticleContainer.clear();
			}
            
            this->clearParticles();
        }
        
        
		/**********************************************************************/
		void remesh()
        {
			if (use_redistribution)
            {
				if(!(runID%use_redistribution))
                {
					double t0=clock();
					model::cout<<"		Remeshing network... "<<std::flush;
					DislocationNetworkRemesh<DislocationNetworkType>(*this).remesh();
					model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]."<<defaultColor<<std::endl;
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
            
            // Create a file-reader object
            EigenDataReader EDR;
            
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
            MatrixDimD C2Gtemp;
            EDR.readMatrixInFile(fullName.str(),"C2G",C2Gtemp); // crystal-to-global orientation
            Material<Isotropic>::rotateCrystal(C2Gtemp);
            
            
            // Min SubNetwork::nodeOrder for Assemble and solve
            int minSNorderForSolve_temp(0);
            EDR.readScalarInFile(fullName.str(),"minSNorderForSolve",minSNorderForSolve_temp); // material by atomic number Z
            assert(minSNorderForSolve_temp>=0 && "minSNorderForSolve must be >=0");
            shared.minSNorderForSolve=(size_t)minSNorderForSolve_temp;
            
            // QuadratureParticle
            EDR.readScalarInFile(fullName.str(),"coreSize",StressField::a); // core-width
            assert((StressField::a)>0.0 && "coreSize MUST BE > 0.");
            StressField::a2=StressField::a*StressField::a;
            LinkType::coreLsquared=StressField::a2;
            //            EDR.readScalarInFile(fullName.str(),"useMultipoleStress",DislocationQuadratureParticle<dim,cellSize>::useMultipoleStress); // useMultipoleStress
            
            // Multipole Expansion
            EDR.readScalarInFile(fullName.str(),"dislocationCellSize",SpatialCellObserverType::cellSize); // cellSize
            
            //            EDR.readScalarInFile(fullName.str(),"nearCellStressApproximation",DislocationParticleType::nearCellStressApproximation); // useMultipoleStress
            //            EDR.readScalarInFile(fullName.str(),"farCellStressApproximation",DislocationParticleType::farCellStressApproximation); // useMultipoleStress
            //            assert((DislocationParticleType::farCellStressApproximation >= DislocationParticleType::nearCellStressApproximation) && "NEAR-FIELD APPROXIMATION IS COARSER THAN FAR-FIELD APPROXIMATION");
            
            
            EDR.readMatrixInFile(fullName.str(),"externalStress",shared.externalStress);
			
            // Implicit time integration
            EDR.readScalarInFile(fullName.str(),"useImplicitTimeIntegration",useImplicitTimeIntegration);
            
			// Restart
            EDR.readScalarInFile(fullName.str(),"startAtTimeStep",runID);
            
            // IO
            EDR.readScalarInFile(fullName.str(),"outputFrequency",DislocationNetworkIO<DislocationNetworkType>::outputFrequency);
            EDR.readScalarInFile(fullName.str(),"outputBinary",DislocationNetworkIO<DislocationNetworkType>::outputBinary);
            EDR.readScalarInFile(fullName.str(),"outputGlidePlanes",DislocationNetworkIO<DislocationNetworkType>::outputGlidePlanes);
            EDR.readScalarInFile(fullName.str(),"outputSpatialCells",DislocationNetworkIO<DislocationNetworkType>::outputSpatialCells);
            EDR.readScalarInFile(fullName.str(),"outputPKforce",DislocationNetworkIO<DislocationNetworkType>::outputPKforce);
            EDR.readScalarInFile(fullName.str(),"outputMeshDisplacement",DislocationNetworkIO<DislocationNetworkType>::outputMeshDisplacement);
            EDR.readScalarInFile(fullName.str(),"outputElasticEnergy",DislocationNetworkIO<DislocationNetworkType>::outputElasticEnergy);
            
			// Mesh and BVP
			EDR.readScalarInFile(fullName.str(),"boundary_type",shared.boundary_type);
			if (shared.boundary_type)
            {
                int meshID(0);
                EDR.readScalarInFile(fullName.str(),"meshID",meshID);
                shared.mesh.readMesh(meshID);
                assert(shared.mesh.size() && "MESH IS EMPTY.");
                
				EDR.readScalarInFile(fullName.str(),"use_bvp",shared.use_bvp);
				if(shared.use_bvp)
                {
                    shared.bvpSolver.init();
				}
			}
			else{ // no boundary is used, DislocationNetwork is in inifinite medium
				shared.use_bvp=0;	// never comupute boundary correction
			}
			
			dt=0.0;
			EDR.readScalarInFile(fullName.str(),"dx",dx);
			assert(dx>0.0);
            EDR.readScalarInFile(fullName.str(),"equilibriumVelocity",equilibriumVelocity);
			assert(equilibriumVelocity>=0.0);
            
			
			EDR.readScalarInFile(fullName.str(),"Nsteps",Nsteps);
			assert(Nsteps>=0 && "Nsteps MUST BE >= 0");
			
			EDR.readScalarInFile(fullName.str(),"timeWindow",timeWindow);
			assert(timeWindow>=0.0 && "timeWindow MUST BE >= 0");
            
            // JUNCTION FORMATION
            EDR.readScalarInFile(fullName.str(),"use_junctions",use_junctions);
            EDR.readScalarInFile(fullName.str(),"collisionTol",DislocationJunctionFormation<DislocationNetworkType>::collisionTol);
            
            // VERTEX REDISTRIBUTION
            EDR.readScalarInFile(fullName.str(),"use_redistribution",use_redistribution);
            EDR.readScalarInFile(fullName.str(),"Lmin",DislocationNetworkRemesh<DislocationNetworkType>::Lmin);
            assert(DislocationNetworkRemesh<DislocationNetworkType>::Lmin>=0.0);
            //            DislocationNetworkRemesh<DislocationNetworkType>::Lmin=2.0*dx;
            assert(DislocationNetworkRemesh<DislocationNetworkType>::Lmin>=2.0*dx && "YOU MUST CHOOSE Lmin>2*dx.");
            EDR.readScalarInFile(fullName.str(),"Lmax",DislocationNetworkRemesh<DislocationNetworkType>::Lmax);
            assert(DislocationNetworkRemesh<DislocationNetworkType>::Lmax>DislocationNetworkRemesh<DislocationNetworkType>::Lmin);
            EDR.readScalarInFile(fullName.str(),"thetaDeg",DislocationNetworkRemesh<DislocationNetworkType>::thetaDeg);
            assert(DislocationNetworkRemesh<DislocationNetworkType>::thetaDeg>=0.0);
            assert(DislocationNetworkRemesh<DislocationNetworkType>::thetaDeg<=90.0);
            
            // Cross-Slip
            EDR.readScalarInFile(fullName.str(),"use_crossSlip",use_crossSlip);
            if(use_crossSlip)
            {
                EDR.readScalarInFile(fullName.str(),"crossSlipDeg",crossSlipDeg);
                assert(crossSlipDeg>=0.0 && crossSlipDeg <= 90.0 && "YOU MUST CHOOSE 0.0<= crossSlipDeg <= 90.0");
                EDR.readScalarInFile(fullName.str(),"crossSlipLength",crossSlipLength);
                assert(crossSlipLength>=DislocationNetworkRemesh<DislocationNetworkType>::Lmin && "YOU MUST CHOOSE crossSlipLength>=Lmin.");
				//				assert(crossSlipLength<Lmin && "YOU MUST CHOOSE crossSlipLength<Lmin."); // Because otherwise cross-slip points would go outside segment
            }
			
            // Read Vertex and Edge information
            DislocationNetworkIO<DislocationNetworkType>::readVertices(*this,runID); // this requires mesh to be up-to-date
            DislocationNetworkIO<DislocationNetworkType>::readEdges(*this,runID);    // this requires mesh to be up-to-date
            
            
            if (shared.use_bvp && (shared.boundary_type==softBoundary))
            { // MOVE THIS WITH REST OB BVP STUFF
                //                shared.vbsc.read(runID,&shared);
                shared.vbsc.initializeVirtualSegments(*this);
            }
            
            
            //#ifdef DislocationNucleationFile
            //            EDR.readScalarInFile(fullName.str(),"nucleationFreq",nucleationFreq);
            //#endif
            
#ifdef _MODEL_MPI_
            // Avoid that a processor starts writing before other are reading
            MPI_Barrier(MPI_COMM_WORLD);
#endif
			
			// Initializing initial configuration
//			model::cout<<redBoldColor<<"runID "<<runID<<" (initial configuration). nodeOrder="<<this->nodeOrder()<<", linkOrder="<<this->linkOrder()<<defaultColor<<std::endl;
			move(0.0,0.0);	// initial configuration
//			output();	// initial configuration, this overwrites the input file
//			if (runID==0) // not a restart
//            {
//				remesh();	// expand initial FR sources
//			}
//			updateQuadraturePoints();
//			++runID;     // increment the runID counter
        }
		
		/**********************************************************************/
		void assembleAndSolve()
        {/*! Performs the following operatons:
          */
            
            model::cout<<"		Assembling edge stiffness and force vectors..."<<std::flush;
			double t0=clock();
            
            //! -1 Compute the interaction StressField between dislocation particles
            this->template computeNeighborField<StressField>();
            
			//! -2 Loop over DislocationSegments and assemble stiffness matrix and force vector
			typedef void (LinkType::*LinkMemberFunctionPointerType)(void); // define type of Link member function
			LinkMemberFunctionPointerType Lmfp(&LinkType::assemble); // Lmfp is a member function pointer to Link::assemble
			this->parallelExecute(Lmfp);
			model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]."<<defaultColor<<std::endl;
            
			//! -3 Loop over DislocationSubNetworks, assemble subnetwork stiffness matrix and force vector, and solve
			model::cout<<"		Solving..."<<std::flush;
			t0=clock();
            
#ifdef _OPENMP
#pragma omp parallel for
            for (unsigned int k=0;k<this->Naddresses();++k)
            {
                typename NetworkComponentContainerType::iterator snIter(this->ABbegin()); //  data within a parallel region is private to each thread
                std::advance(snIter,k);
                if (snIter->second->nodeOrder()>=shared.minSNorderForSolve)
                {
                    DislocationNetworkComponentType(*snIter->second).sparseSolve();
                }
            }
#else
			for (typename NetworkComponentContainerType::iterator snIter=this->ABbegin(); snIter!=this->ABend();++snIter)
            {
                if (snIter->second->nodeOrder()>=shared.minSNorderForSolve)
                {
                    DislocationNetworkComponentType(*snIter->second).sparseSolve();
                }
			}
#endif
			model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]."<<defaultColor<<std::endl;
		}
		
		/**********************************************************************/
		void output() const
        {/*! Outputs DislocationNetwork data
          */
            double t0=clock();
            if (!(runID%DislocationNetworkIO<DislocationNetworkType>::outputFrequency))
            {
#ifdef _MODEL_DD_MPI_
				if(ModelMPIbase::mpiRank()==0)
                {
                    DislocationNetworkIO<DislocationNetworkType>::output(*this,runID);
				}
#else
                DislocationNetworkIO<DislocationNetworkType>::output(*this,runID);
#endif
			}
			model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]."<<defaultColor<<std::endl;
		}
		
		
		/**********************************************************************/
		void updateQuadraturePoints()
        {
            model::cout<<"		Updating quadrature points... "<<std::flush;
			double t0=clock();
            
            // Clear DislocationParticles
            this->clearParticles();
            
            // Populate DislocationParticles
            for (typename NetworkLinkContainerType::iterator linkIter =this->linkBegin();
                 /*                                       */ linkIter!=this->linkEnd();
                 /*                                     */ ++linkIter)
            {
                linkIter->second->updateQuadraturePoints(*this);
			}
			model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]."<<defaultColor<<std::endl;
		}
		
		
		/**********************************************************************/
		void move(const double & dt_in, const double & dt_old)
        {/*! Moves all nodes in the DislocationNetwork using the stored velocity and current dt
          */
			model::cout<<"		Moving DislocationNodes (dt="<<dt_in<< ")... "<<std::flush;
			double t0=clock();
            //			typedef void (NodeType::*NodeMemberFunctionPointerType)(const double&); // define type of Link member function
            //			NodeMemberFunctionPointerType Nmfp(&NodeType::move); // Lmfp is a member function pointer to Link::assemble
            //			double t0=clock();
            //			this->parallelExecute(Nmfp,dt_in); // NOT POSSIBLE TO PARALLELIZE SINCE move exectuts tansmit::makeT
            
			for (typename NetworkNodeContainerType::iterator nodeIter=this->nodeBegin();nodeIter!=this->nodeEnd();++nodeIter)
            {
				nodeIter->second->move(dt_in,dt_old);
			}
			model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]."<<defaultColor<<std::endl;
		}
		
		/**********************************************************************/
		void runSteps()
        {/*! Runs Nsteps simulation steps
          */
			double ts(clock());
			for (int k=0;k<Nsteps;++k)
            {
				model::cout<<std::endl; // leave a blank line
				model::cout<<blueBoldColor<<"Step "<<k+1<<" of "<<Nsteps<<defaultColor<<std::endl;
				singleStep();
			}
			model::cout<<greenBoldColor<<std::setprecision(3)<<std::scientific<<Nsteps<< " simulation steps completed in "<<(clock()-ts)/CLOCKS_PER_SEC<<" [sec]"<<defaultColor<<std::endl;
		}
		
		/**********************************************************************/
		void runTime()
        {/*! Runs a number simulation steps corresponding to a total
          * dimensionless time timeWindow
          */
			double ts(clock());
			double elapsedTime(0.0);
			while (elapsedTime<timeWindow)
            {
				model::cout<<std::endl; // leave a blank line
				model::cout<<blueBoldColor<<"Time "<<elapsedTime<<" of "<<timeWindow<<defaultColor<<std::endl;
				singleStep();
				elapsedTime+=dt;
			}
			model::cout<<greenBoldColor<<std::setprecision(3)<<std::scientific<<timeWindow<< " simulation time completed in "<<(clock()-ts)/CLOCKS_PER_SEC<<" [sec]"<<defaultColor<<std::endl;
		}
		
		/**********************************************************************/
		void checkBalance() const // TO DO: MOVE THIS TO NETWORK LAYER
        {/*! Checks that each DislocationNode is balanced, and asserts otherwise
          * Exceptions are nodes with only one neighbors (FR source)
          * and nodes on the boundary.
          */
			for (typename NetworkNodeContainerType::const_iterator nodeIter =this->nodeBegin();
                 /*                                             */ nodeIter!=this->nodeEnd();
                 /*                                             */ nodeIter++)
            {
				if (nodeIter->second->neighborhood().size()>2)
                {
                    const bool nodeIsBalanced(nodeIter->second->is_balanced());
                    if (!nodeIsBalanced && nodeIter->second->meshLocation()==insideMesh)
                    {
                        model::cout<<"Node "<<nodeIter->second->sID<<" is not balanced:"<<std::endl;
                        model::cout<<"    outflow="<<nodeIter->second->outFlow().transpose()<<std::endl;
                        model::cout<<"     inflow="<<nodeIter->second->inFlow().transpose()<<std::endl;
                        assert(0 && "NODE IS NOT BALANCED");
                    }
				}
			}
		}
				
		/**********************************************************************/
		double energy()
        {/*!\returns The total elastic energy of the dislocation network
          */
			double temp(0.0);
			for (typename NetworkLinkContainerType::iterator linkIter=this->linkBegin();linkIter!=this->linkEnd();++linkIter)
            {
				temp+= linkIter->second->energy();
			}
			return 0.5*temp;
		}
		
		/**********************************************************************/
		//template <bool useFullField=true>
		MatrixDimD stress(const VectorDimD& Rfield, const bool& useFullField=true) const
        {
			MatrixDimD temp=MatrixDimD::Zero();
			if(useFullField)
            {
				for (typename NetworkLinkContainerType::const_iterator linkIter=this->linkBegin();linkIter!=this->linkEnd();++linkIter)
                {
					temp+= linkIter->second->stress_source(Rfield);
				}
			}
			else
            {
                assert(0 && "RE-ENABLE THIS WITH NEW CELL CLASS");
                //				SpatialCellObserverType sCO;
                //				typename SpatialCellObserverType::SharedPtrType sharedCellPointer(sCO.getCellByPosition(Rfield));
                //				for(typename SpatialCellType::ParticleContainerType::const_iterator particleIter =sharedCellPointer->neighborParticleContainer.begin();
                //					/*                                                         */ particleIter!=sharedCellPointer->neighborParticleContainer.end();
                //					/*                                                         */ ++particleIter){
                //					temp+=(*particleIter)->stress_at(Rfield);
                //				}
			}
            if (shared.use_bvp && (shared.boundary_type==1))
            {
                temp+= shared.vbsc.stress(Rfield);
            }
			return temp;
		}
        
        
//        /**********************************************************************/
//        void getStress(std::deque<SimpleFieldPoint>& fieldPoints)
//        {
//            DN.updateQuadraturePoints();
//            DN.computeField<SimpleFieldPoint,StressField>(fieldPoints);
//        }
        
		
		/**********************************************************************/
		VectorDimD displacement(const VectorDimD & Rfield,const VectorDimD & S) const
        {
			VectorDimD temp(VectorDimD::Zero());
			for (typename NetworkLinkContainerType::const_iterator linkIter =this->linkBegin();
                 /*                                             */ linkIter!=this->linkEnd();
                 /*                                             */ linkIter++)
            {
				temp+= linkIter->second->displacement(Rfield,S);
			}
			return temp;
		}
		
        /**********************************************************************/
		MatrixDimD plasticDistortionRate() const
        {
			MatrixDimD temp(MatrixDimD::Zero());
			for (typename NetworkLinkContainerType::const_iterator linkIter =this->linkBegin();
                 /*                                             */ linkIter!=this->linkEnd();
                 /*                                             */ linkIter++)
            {
				temp+= linkIter->second->plasticDistortionRate();
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
		MatrixDimD latticeRotation(const VectorDimD & Rfield) const
        {
			MatrixDimD temp(MatrixDimD::Zero());
			for (typename NetworkLinkContainerType::const_iterator linkIter =this->linkBegin();
                 /*                                             */ linkIter!=this->linkEnd();
                 /*                                             */ linkIter++)
            {
				temp+= linkIter->second->lattice_rotation_source(Rfield);
			}
			return temp;
		}
		
        /**********************************************************************/
		MatrixDimD elasticDistortion(const VectorDimD & Rfield) const
        {
			MatrixDimD temp(MatrixDimD::Zero());
			for (typename NetworkLinkContainerType::const_iterator linkIter =this->linkBegin();
                 /*                                             */ linkIter!=this->linkEnd();
                 /*                                             */ linkIter++)
            {
				temp+= linkIter->second->displacement_gradient_source(Rfield);
			}
			return temp;
		}
		
        /**********************************************************************/
		double networkLength() const
        {/*!\returns the line length of *this DislocationNetwork.
          */
			double temp(0.0);
			for (typename NetworkLinkContainerType::const_iterator linkIter =this->linkBegin();
                 /*                                             */ linkIter!=this->linkEnd();
                 /*                                             */ linkIter++)
            {
				temp+= linkIter->second->template arcLength<qOrder,QuadratureRule>();
			}
			return temp;
		}
        
        /**********************************************************************/
        const unsigned int& runningID() const
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
		
	};
    
//    // decalre static data
//    template <short unsigned int _dim, short unsigned int corder, typename InterpolationType,
//	/*	   */ short unsigned int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule>
//    bool DislocationNetwork<_dim,corder,InterpolationType,qOrder,QuadratureRule>::useImplicitTimeIntegration=false;
//    
//    template <short unsigned int _dim, short unsigned int corder, typename InterpolationType,
//	/*	   */ short unsigned int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule>
//    double DislocationNetwork<_dim,corder,InterpolationType,qOrder,QuadratureRule>::equilibriumVelocity=0.01;
    
} // namespace model
#endif

//		/**********************************************************************/
//		MatrixDimD stressFromGlidePlane(const Eigen::Matrix<double,dim+1,1>& key, const VectorDimD& Rfield) const
//        {
//			//GlidePlaneObserver<dim,LinkType> gpObsever;
//			MatrixDimD temp(MatrixDimD::Zero());
//			typedef typename GlidePlaneObserver<LinkType>::GlidePlaneType GlidePlaneType;
//			std::pair<bool, const GlidePlaneType* const> isGp(this->isGlidePlane(key));
//			if (isGp.first)
//            {
//				temp=isGp.second->stress(Rfield);
//			}
//            if (shared.use_bvp && (shared.boundary_type==1))
//            {
//                temp+= shared.vbsc.stressFromGlidePlane(key,Rfield);
//            }
//			return temp;
//		}

