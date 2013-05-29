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



// BEING MODIFIED
// 1- BIN WRITER/READER. DEFINE DISLOCATIONSEGMENT::OUTDATA as Eigen::Matrix<double,1,9>. Then in Network add the function friend void operator<< (SequentialBinFile<template Char,???>, ...)
// 2- Zero node Tangent Problem
// BVP: integration method
// 3 crossSlip
// Motion of BoundarySubnetworks

// TO DO
// -3 MAKE isIsolated and isBalanced data members of NetworkNode that are modified by addToNeighbors and removeFromNeigbors
// implement IndependentPath

// -1 - remove DislocationEnergyRules and EdgeConfig !!!
// 0 - Finish DepthFirst class. Don't allow to search/execute if N=0, so that N=1 means that node only, n=2 means first neighbor. Or change SpineNodeBase_CatmullRom::TopologyChangeActions
// 1- Implement operator << SpatialCell
// 2- remove AddressBook, wherever possible, Done in Node chain
// 40 - clean MultiExpand in Network Layer, remove get_r from there
// 35- Simplify Neighborhood structure
// 25- Remove template parameter alpha and make template member function
// 14- RENAME ORIGINAL MESH FOLDER /M
// 16- READ/WRITE IN BINARY FORMAT
// 18- Should define linear=1, quadratic=2, cubic=3 and use polyDegree instead of corder. Put corder in SplineEnums
// 37- IS PLANAR SHOULD RETURN 0 IF IS A LINE!!!!! CHANGE ALSO IN SPLINESEGMENTBASE
// 38- IS PLANAR SHOULD RETURN the normal as pair<bool,normal>
// 37- NetworkNode, initializations from expansion and contraction
// IMPLEMENT NEIGHBOR ITERATORS IN NETWORKNODE
// CHANGE CONST VERSION OF EDGEFINDER/VERTEXFINDER PASS this IN MEMBER FUNCTION, REMOVE TEMPLATE SPECIALIZATION AND MAKE STATIC FUNCTIONS
// 9- cellSize should depend on applied load. Or better the number of cell neighbors used in each cell should depend on the applied stress to that cell


// ISSUES
// 3- SplineIntersection::planePlaneType uses the wrong tolerance (10x)
// 10- SplineNetworkBase_common :		assert(iter->second.size()==2);
//     SplineIntersection: //			assert(T0.norm()>tol && "SplineIntersection<3,3,1,2>: T0 too small.");  // NOT RIGHT FOR DIPOLAR LOOPS
//									assert(T1.norm()>tol && "SplineIntersection<3,3,1,2>: T1 too small.");	// NOT RIGHT FOR DIPOLAR LOOPS
// 14- If T0=0 or T1=0, then rl(0)=0/0=NaN !!!!! This is not true since rl still tends to a finite vector. Remove class Parametric curve and implement special case of rl at 0 and 1 for vanishing nodal tangents


#ifndef model_DISLOCATIONNETWORK_H_
#define model_DISLOCATIONNETWORK_H_



#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#endif

#define EIGEN_DONT_PARALLELIZE // disable Eigen Internal openmp Parallelization
#include <Eigen/Dense>

#include <math.h>
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

#ifdef _MODEL_DD_MPI_
#define _MODEL_MPI_  // required in ParticleSystem.h
#include <model/ParticleInteraction/ParticleSystem.h> // the main object from pil library
#endif




namespace model {
	
	
	/**************************************************************************/
	/**************************************************************************/
	template <short unsigned int _dim, short unsigned int corder, typename InterpolationType,
	/*	   */ double & alpha, short unsigned int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule>
	class DislocationNetwork :
    /* inheritance          */ public Network<DislocationNetwork<_dim,corder,InterpolationType,alpha,qOrder,QuadratureRule> >,
	/* inheritance          */ public GlidePlaneObserver<typename TypeTraits<DislocationNetwork<_dim,corder,InterpolationType,alpha,qOrder,QuadratureRule> >::LinkType>
#ifdef _MODEL_DD_MPI_
    /* inheritance          */, private ParticleSystem<DislocationQuadratureParticle<_dim> >
#endif
    {
		
    public:
        
        enum {dim=_dim}; // make dim available outside class
        
		typedef DislocationNetwork<dim,corder,InterpolationType,alpha,qOrder,QuadratureRule> DislocationNetworkType;
		typedef DislocationNetworkType Derived; // define Derived to use NetworkTypedefs.h
#include <model/Network/NetworkTypedefs.h>
        typedef DislocationNetworkComponent<NodeType,LinkType> DislocationNetworkComponentType;
		typedef Eigen::Matrix<double,dim,dim>	MatrixDimD;
		typedef Eigen::Matrix<double,dim,1>		VectorDimD;
		typedef typename LinkType::DislocationQuadratureParticleType DislocationQuadratureParticleType;
		typedef DislocationCell<dim> SpatialCellType;
		typedef SpatialCellObserver<SpatialCellType,dim> SpatialCellObserverType;
		typedef typename SpatialCellObserverType::CellMapType CellMapType;
		typedef GlidePlaneObserver<LinkType> GlidePlaneObserverType;
#ifdef _MODEL_DD_MPI_
        typedef ParticleSystem<DislocationQuadratureParticle<_dim> > ParticleSystemType;
#endif
		enum {NdofXnode=NodeType::NdofXnode};
		
#ifdef UpdateBoundaryConditionsFile
#include UpdateBoundaryConditionsFile
#endif
        
#ifdef DislocationNucleationFile
#include DislocationNucleationFile
//         int nucleationFreq;
#endif
		
		
	private:
        
        
        
		
		short unsigned int use_redistribution;
		bool use_junctions;
		static bool useImplicitTimeIntegration;
        
		unsigned int runID;
				
		bool use_crossSlip;
		double crossSlipDeg;
		double crossSlipLength;
		
		double totalTime;
		
		double dx, dt;
        double vmax;
        
		

        
        
		/* formJunctions ******************************************************/
		void formJunctions()
        {/*! Performs dislocation junction formation if use_junctions==true
          */
			if (use_junctions)
            {
				double t0=clock();
				std::cout<<"		Forming Junctions: found (";
                const double avoidNodeIntersection(0.05);
				DislocationJunctionFormation<DislocationNetworkType>(*this).formJunctions(dx,avoidNodeIntersection);
				std::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]."<<defaultColor<<std::endl;
			}
		}
		
		/* make_dt ************************************************************/
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
			std::cout<<"		Computing dt..."<<std::flush;
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
			
			double shearWaveFraction(0.1);
			//short unsigned int shearWaveExp=1;
			if (vmax > Material<Isotropic>::cs*shearWaveFraction)
            {
				dt=dx/vmax;
			}
			else
            {
                //dt=dx/std::pow(shared.material.cs*shearWaveFraction,shearWaveExp+1)*std::pow(vmax,shearWaveExp);
				//dt=dx/(shared.material.cs*shearWaveFraction)*std::pow(vmax/(shared.material.cs*shearWaveFraction),1);
				dt=dx/(Material<Isotropic>::cs*shearWaveFraction);
			}
			std::cout<<std::setprecision(3)<<std::scientific<<" vmax="<<vmax;
			std::cout<<std::setprecision(3)<<std::scientific<<" dt="<<dt;
			std::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]."<<defaultColor<<std::endl;
		}
				
		/* crossSlip **********************************************************/
		void crossSlip()
        {/*! Performs dislocation cross-slip if use_crossSlip==true
          */
			if(use_crossSlip)
            {
				double t0=clock();
				std::cout<<"		Performing Cross Slip ... "<<std::flush;
				size_t crossSlipEvents(DislocationCrossSlip<DislocationNetworkType>(*this).crossSlip(crossSlipDeg,crossSlipLength));
				std::cout<<crossSlipEvents<<" cross slip events found ";
				std::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]."<<defaultColor<<std::endl;
			}
		}
        
        /* update_BVP_Solution ************************************************/
        void update_BVP_Solution(const bool& updateUserBC)
        {
            if (shared.use_bvp)
            {
				double t0=clock();
				std::cout<<"		Updating bvp stress ... ";
				if(!(runID%shared.use_bvp))
                {
					shared.domain.update_BVP_Solution(updateUserBC,this);
				}
				for (typename NetworkNodeContainerType::iterator nodeIter=this->nodeBegin();nodeIter!=this->nodeEnd();++nodeIter)
                { // THIS SHOULD BE PUT IN THE MOVE FUNCTION OF DISLOCATIONNODE
					nodeIter->second->updateBvpStress();
				}
				std::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]."<<defaultColor<<std::endl;
			}
        }
		        
        /*  singleStep*********************************************************/
		void singleStep(const bool& updateUserBC=false)
        {
			//! A simulation step consists of the following:
			std::cout<<blueBoldColor<< "runID="<<runID
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
            update_BVP_Solution(updateUserBC);
			
#ifdef _MODEL_DD_MPI_
            MPI_Barrier(MPI_COMM_WORLD);
#endif

			//! 3- Solve the equation of motion
			assembleAndSolve();
			
    
			//! 4- Compute time step dt (based on max nodal velocity) and increment totalTime
			make_dt();
			totalTime+=dt;
            
            const MatrixDimD pdr(plasticDistortionRate());
            plasticDistortion += pdr*dt;
            
            
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
//			loopInversion();
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
		
		/* remeshByContraction ************************************************/
		void removeBoundarySegments()
        {/*! Removes DislocationSegment(s) on the mesh boundary
          */
			if (shared.boundary_type==softBoundary)
            {
				double t0=clock();
				std::cout<<"		Removing Segments outside Mesh Boundaries... ";
				typedef bool (LinkType::*link_member_function_pointer_type)(void) const;
				link_member_function_pointer_type boundarySegment_Lmfp;
				boundarySegment_Lmfp=&LinkType::is_boundarySegment;
				this->template disconnect_if<1>(boundarySegment_Lmfp);
				std::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]."<<defaultColor<<std::endl;
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
        
		/**********************************************************************/
		/* PUBLIC SECTION *****************************************************/
	public:
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
		
		DislocationSharedObjects<LinkType> shared;
		
		//! The number of simulation steps taken by the next call to runByStep()
		int Nsteps;
		double timeWindow;
        
        MatrixDimD plasticDistortion;
        
		/* Constructor ********************************************************/
        DislocationNetwork(int& argc, char* argv[])
#ifdef _MODEL_DD_MPI_
        /* base initialization */ : ParticleSystemType(argc,argv)
#endif
        
        {
//        	if (argc==1){
                // argv[0] is by default the name of the executable so use default name for inputfile
                read("./","DDinput.txt");
//            }
//            else{
//                // argv[1] is assumed to be the filename with working directory the current directory
//                read("./",argv[1]);
//            }
            
            plasticDistortion.setZero();
        }
        
		/* remesh *************************************************************/
		void remesh()
        {
			if (use_redistribution)
            {
				if(!(runID%use_redistribution))
                {
					double t0=clock();
					std::cout<<"		remeshing network... "<<std::flush;
					DislocationNetworkRemesh<DislocationNetworkType>(*this).remesh();
					std::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]."<<defaultColor<<std::endl;
				}
			}
		}
		
		/* get_dt *************************************************************/
		const double& get_dt() const
        {/*! The current simulation time step in dimensionless units
          */
			return dt;
		}
        
        /* get_totalTime ******************************************************/
		const double& get_totalTime() const
        {/*! The elapsed simulation time step in dimensionless units
          */
			return totalTime;
		}
		
		/* read ***************************************************************/
		void read(const std::string& inputDirectoryName_in, std::string inputFileName)
        { // TO DO: move this to DislocationNetworkIO.h
			
            std::ostringstream fullName;
            fullName<<inputDirectoryName_in<<inputFileName;
            
            // Create a file-reader object
            EigenDataReader EDR;
            
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
            EDR.readScalarInFile(fullName.str(),"coreWidthSquared",DislocationQuadratureParticle<dim>::a2); // core-width
            assert((DislocationQuadratureParticle<dim>::a2)>0.0 && "coreWidthSquared MUST BE > 0.");
            LinkType::coreLsquared=DislocationQuadratureParticle<dim>::a2;
            //            EDR.readScalarInFile(fullName.str(),"useMultipoleStress",DislocationQuadratureParticle<dim,cellSize>::useMultipoleStress); // useMultipoleStress
            
            // Multipole Expansion
            EDR.readScalarInFile(fullName.str(),"dislocationCellSize",SpatialCellObserverType::cellSize); // cellSize
            EDR.readScalarInFile(fullName.str(),"nearCellStressApproximation",DislocationQuadratureParticle<dim>::nearCellStressApproximation); // useMultipoleStress
            EDR.readScalarInFile(fullName.str(),"farCellStressApproximation",DislocationQuadratureParticle<dim>::farCellStressApproximation); // useMultipoleStress
            assert((DislocationQuadratureParticle<dim>::farCellStressApproximation >= DislocationQuadratureParticle<dim>::nearCellStressApproximation) && "NEAR-FIELD APPROXIMATION IS COARSER THAN FAR-FIELD APPROXIMATION");
            
            
            EDR.readMatrixInFile(fullName.str(),"externalStress",shared.externalStress);
			
            // Implicit time integration
            EDR.readScalarInFile(fullName.str(),"useImplicitTimeIntegration",useImplicitTimeIntegration);
        
			// Restart
            EDR.readScalarInFile(fullName.str(),"startAtTimeStep",runID);
			
            // IO
            EDR.readScalarInFile(fullName.str(),"outputFrequency",DislocationNetworkIO<DislocationNetworkType>::outputFrequency);
            EDR.readScalarInFile(fullName.str(),"outputGlidePlanes",DislocationNetworkIO<DislocationNetworkType>::outputGlidePlanes);
            EDR.readScalarInFile(fullName.str(),"outputSpatialCells",DislocationNetworkIO<DislocationNetworkType>::outputSpatialCells);
            EDR.readScalarInFile(fullName.str(),"outputPKforce",DislocationNetworkIO<DislocationNetworkType>::outputPKforce);
            EDR.readScalarInFile(fullName.str(),"outputMeshDisplacement",DislocationNetworkIO<DislocationNetworkType>::outputMeshDisplacement);
            
            
			// BVP
			EDR.readScalarInFile(fullName.str(),"boundary_type",shared.boundary_type);
			if (shared.boundary_type)
            {
				EDR.readScalarInFile(fullName.str(),"use_bvp",shared.use_bvp);
                shared.domain.readMesh();
				if(shared.use_bvp)
                {
					shared.domain.readInputBCs();
				}
			}
			else{ // no boundary is used, this means dislocation network in inifinite medium
				shared.use_bvp=0;	// never comupute boundary correction
			}
			
			dt=0.0;
			EDR.readScalarInFile(fullName.str(),"dx",dx);
			assert(dx>0.0);
			
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
            assert(DislocationNetworkRemesh<DislocationNetworkType>::Lmin>2.0*dx && "YOU MUST CHOOSE Lmin>2*dx.");
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
//            readNodes(runID);
            DislocationNetworkIO<DislocationNetworkType>::readVertices(*this,runID);
            //readLinks(runID);
            DislocationNetworkIO<DislocationNetworkType>::readEdges(*this,runID);

            
            if (shared.use_bvp && (shared.boundary_type==softBoundary))
            { // MOVE THIS WITH REST OB BVP STUFF
                //                shared.vbsc.read(runID,&shared);
                shared.vbsc.initializeVirtualSegments(*this);
            }
            
            
//#ifdef DislocationNucleationFile
//            EDR.readScalarInFile(fullName.str(),"nucleationFreq",nucleationFreq);
//#endif

            
			
			// Initializing initial configuration
			std::cout<<redBoldColor<<"runID "<<runID<<" (initial configuration). nodeOrder="<<this->nodeOrder()<<", linkOrder="<<this->linkOrder()<<defaultColor<<std::endl;
			move(0.0,0.0);	// initial configuration
			output();	// initial configuration, this overwrites the input file
			if (runID==0) // not a restart
            { 
				remesh();	// expand initial FR sources
			}
			updateQuadraturePoints();
			++runID;     // increment the runID counter
		}
		
		/* solve **************************************************************/
		void assembleAndSolve()
        {/*! Assemble and solve equation system
          */
            
			//! 1- Loop over DislocationSegments and assemble stiffness matrix and force vector
			std::cout<<"		Assembling edge stiffness and force vectors..."<<std::flush;
			typedef void (LinkType::*LinkMemberFunctionPointerType)(void); // define type of Link member function
			LinkMemberFunctionPointerType Lmfp(&LinkType::assemble); // Lmfp is a member function pointer to Link::assemble
			double t0=clock();
			this->parallelExecute(Lmfp);
			std::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]."<<defaultColor<<std::endl;
            
			//! 2- Loop over DislocationSubNetworks, assemble subnetwork stiffness matrix and force vector, and solve
			std::cout<<"		Solving..."<<std::flush;
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
			std::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]."<<defaultColor<<std::endl;
		}
		
		/* output *************************************************************/
		void output() const
        {/*! Outputs DislocationNetwork data
          */
            double t0=clock();
            if (!(runID%DislocationNetworkIO<DislocationNetworkType>::outputFrequency))
            {
#ifdef _MODEL_DD_MPI_
				if(this->mpiRank==0)
                {
                    DislocationNetworkIO<DislocationNetworkType>::outputTXT(*this,runID);
				}
#else
                DislocationNetworkIO<DislocationNetworkType>::outputTXT(*this,runID);
#endif
			}
			std::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]."<<defaultColor<<std::endl;
		}
		
		
		/* updateQuadraturePoints *********************************************/
		void updateQuadraturePoints()
        {
			double t0=clock();
			std::cout<<"		Updating Quadrature Points... "<<std::flush;
			for (typename NetworkLinkContainerType::iterator linkIter=this->linkBegin();linkIter!=this->linkEnd();++linkIter)
            {
				linkIter->second->quadratureParticleVector.clear(); // first clear all quadrature points for all segments
			}
            for (typename NetworkLinkContainerType::iterator linkIter=this->linkBegin();linkIter!=this->linkEnd();++linkIter)
            {
				Quadrature<1,qOrder>::execute(linkIter->second,&LinkType::updateQuadGeometryKernel); // then update again
			}
            for (typename CellMapType::const_iterator cellIter=SpatialCellObserverType::cellBegin();cellIter!=SpatialCellObserverType::cellEnd();++cellIter)
            {
                cellIter->second->computeCenterStress();
            }
			std::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]."<<defaultColor<<std::endl;
		}
		
		
		/* move ***************************************************************/
		void move(const double & dt_in, const double & dt_old)
        {/*! Moves all nodes in the DislocationNetwork using the stored velocity and current dt
          */
			std::cout<<"		Moving Dislocation Nodes (dt="<<dt_in<< ")... "<<std::flush;
			double t0=clock();
            //			typedef void (NodeType::*NodeMemberFunctionPointerType)(const double&); // define type of Link member function
            //			NodeMemberFunctionPointerType Nmfp(&NodeType::move); // Lmfp is a member function pointer to Link::assemble
            //			double t0=clock();
            //			this->parallelExecute(Nmfp,dt_in);
            
			for (typename NetworkNodeContainerType::iterator nodeIter=this->nodeBegin();nodeIter!=this->nodeEnd();++nodeIter)
            {
				nodeIter->second->move(dt_in,dt_old);
			}
			std::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]."<<defaultColor<<std::endl;
		}
		
		/**********************************************************************/
		void runByStep(const bool& updateUserBC=false)
        {/*! Runs Nsteps simulation steps
          */
			double ts(clock());
			for (int k=0;k<Nsteps;++k)
            {
				std::cout<<std::endl; // leave a blank line
				std::cout<<blueBoldColor<<"Step "<<k+1<<" of "<<Nsteps<<defaultColor<<std::endl;
				singleStep(updateUserBC);
			}
			std::cout<<greenBoldColor<<std::setprecision(3)<<std::scientific<<Nsteps<< " simulation steps completed in "<<(clock()-ts)/CLOCKS_PER_SEC<<" [sec]"<<defaultColor<<std::endl;
		}
		
		
		/**********************************************************************/
		void runByTime(const bool& updateUserBC=false)
        {/*! Runs a number simulation steps corresponding to a total
          * dimensionless time timeWindow
          */
			double ts(clock());
			double elapsedTime(0.0);
			while (elapsedTime<timeWindow)
            {
				std::cout<<std::endl; // leave a blank line
				std::cout<<blueBoldColor<<"Time "<<elapsedTime<<" of "<<timeWindow<<defaultColor<<std::endl;
				singleStep(updateUserBC);
				elapsedTime+=dt;
			}
			std::cout<<greenBoldColor<<std::setprecision(3)<<std::scientific<<timeWindow<< " simulation time completed in "<<(clock()-ts)/CLOCKS_PER_SEC<<" [sec]"<<defaultColor<<std::endl;
		}
		
		
		/**********************************************************************/
		void checkBalance() const // TO DO: MOVE THIS TO NETWORK LAYER
        {/*! Checks that each DislocationNode is balanced, and asserts otherwise
          * Exceptions are nodes with only one neighbors (FR source)
          * and nodes on the boundary.
          */
			for (typename NetworkNodeContainerType::const_iterator nodeIter=this->nodeBegin();nodeIter!=this->nodeEnd();++nodeIter)
            {
				if (nodeIter->second->neighborhood().size()>2)
                {
                    const bool nodeIsBalanced(nodeIter->second->is_balanced());
                    if (!nodeIsBalanced && nodeIter->second->nodeMeshLocation==insideMesh)
                    {
                        std::cout<<"Node "<<nodeIter->second->sID<<" is not balanced:"<<std::endl;
                        std::cout<<"    outflow="<<nodeIter->second->outFlow().transpose()<<std::endl;
                        std::cout<<"     inflow="<<nodeIter->second->inFlow().transpose()<<std::endl;
                        assert(0 && "NODE IS NOT BALANCED");
                    }
				}
			}
		}
		
		
		/**********************************************************************/
		// energy
		double energy()
        {/*! The total elastic energy of the dislocation network
          */
			double temp(0.0);
			for (typename NetworkLinkContainerType::iterator linkIter=this->linkBegin();linkIter!=this->linkEnd();++linkIter)
            {
				temp+= linkIter->second->energy();
			}
			return 0.5*temp;
		}
		
		/**********************************************************************/
		template <bool useFullField=true>
		MatrixDimD stress(const VectorDimD & Rfield) const
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
        
		/**********************************************************************/
		MatrixDimD stressFromGlidePlane(const Eigen::Matrix<double,dim+1,1>& key, const VectorDimD& Rfield) const
        {
			//GlidePlaneObserver<dim,LinkType> gpObsever;
			MatrixDimD temp(MatrixDimD::Zero());
			typedef typename GlidePlaneObserver<LinkType>::GlidePlaneType GlidePlaneType;
			std::pair<bool, const GlidePlaneType* const> isGp(this->isGlidePlane(key));
			if (isGp.first)
            {
				temp=isGp.second->stress(Rfield);
			}
            if (shared.use_bvp && (shared.boundary_type==1))
            {
                temp+= shared.vbsc.stressFromGlidePlane(key,Rfield);
            }
			return temp;
		}
		
		/**********************************************************************/
		VectorDimD displacement(const VectorDimD & Rfield,const VectorDimD & S) const
        {
			VectorDimD temp(VectorDimD::Zero());
			for (typename NetworkLinkContainerType::const_iterator linkIter=this->linkBegin();linkIter!=this->linkEnd();++linkIter)
            {
				temp+= linkIter->second->displacement(Rfield,S);
			}
			return temp;
		}
		
        
        /**********************************************************************/
		MatrixDimD plasticDistortionRate() const
        {
			MatrixDimD temp(MatrixDimD::Zero());
			for (typename NetworkLinkContainerType::const_iterator linkIter=this->linkBegin();linkIter!=this->linkEnd();++linkIter)
            {
				temp+= linkIter->second->plasticDistortionRate();
			}
			return temp;
		}
        
        
        /**********************************************************************/
		MatrixDimD plasticStrainRate() const
        {
            //			MatrixDimD temp(MatrixDimD::Zero());
            //			for (typename NetworkLinkContainerType::const_iterator linkIter=this->linkBegin();linkIter!=this->linkEnd();++linkIter){
            //				temp+= linkIter->second->plasticStrainRate();
            //			}
            const MatrixDimD temp(plasticDistortionRate());
			return (temp+temp.transpose())*0.5;
		}
		
		/********************************************************/
		MatrixDimD lattice_rotation(const VectorDimD & Rfield) const
        {
			MatrixDimD temp(MatrixDimD::Zero());
			for (typename NetworkLinkContainerType::const_iterator linkIter=this->linkBegin();linkIter!=this->linkEnd();++linkIter)
            {
				temp+= linkIter->second->lattice_rotation_source(Rfield);
			}
			return temp;
		}
		
		/********************************************************/
		MatrixDimD elasticDistortion(const VectorDimD & Rfield) const
        {
			MatrixDimD temp(MatrixDimD::Zero());
			for (typename NetworkLinkContainerType::const_iterator linkIter=this->linkBegin();linkIter!=this->linkEnd();++linkIter)
            {
				temp+= linkIter->second->displacement_gradient_source(Rfield);
			}
			return temp;
		}
		
		/********************************************************/
		// function to return the length of the whole dislocation network
		double network_length() const
        {/*! The length of the dislocation network.
          */
			double temp(0.0);
			for (typename NetworkLinkContainerType::const_iterator linkIter=this->linkBegin();linkIter!=this->linkEnd();++linkIter)
            {
				temp+= linkIter->second->template arcLength<qOrder,QuadratureRule>();
			}
			return temp;
		}
        
        /********************************************************/
        const unsigned int& runningID() const
        {/*! The current simulation step ID.
          */
            return runID;
        }
        
        /********************************************************/
        const double& vMax() const
        {/*! The max vertex velocity.
          */
            return vmax;
        }
		
	};
    
    
    // static data
    template <short unsigned int _dim, short unsigned int corder, typename InterpolationType,
	/*	   */ double & alpha, short unsigned int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule>
    bool DislocationNetwork<_dim,corder,InterpolationType,alpha,qOrder,QuadratureRule>::useImplicitTimeIntegration=false;
    
    
	//////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////
} // namespace model
#endif
