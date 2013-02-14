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
// 2- Zeo node Tangent Problem
// BVP: integration method
// 3 crossSlip
// Motion of BoundarySubnetworks

// TO DO
// -3 MAKE isIsolated and isBalanced data members of NetworkNode that are modified by addToNeighbors and removeFromNeigbors
// implement IndependentPath
// -2: updateQuadraturePoints should be called as part of the topologyChangeActions() and move()

// 0 - Finish DepthFirst class. Don't allow to search/execute if N=0, so that N=1 means that node only, n=2 means first neighbor. Or change SpineNodeBase_CatmullRom::TopologyChangeActions
// 1- Implement operator << in DislocationNode, GlidePlane, SpaceCell
// 2- remove AddressBook, wherever possible, Done in Node chain
// 40 - clean MultiExpand in Network Layer, remove get_r from there
// 35- Simplify Neighborhood structure, remove boost::tuple, (use std::touple? Available in c++11, starting with g++4.7)
// 25- Remove template parameter alpha and make template member function
// 14- RENAME ORIGINAL MESH FOLDER /M
// 16- READ/WRITE IN BINARY FORMAT
// 18- Should define linear=1, quadratic=2, cubic=3 and use polyDegree instead of corder. Put corder in SplineEnums
// 30- Make coreL static
// 37- IS PLANAR SHOULD RETURN 0 IF IS A LINE!!!!! CHANGE ALSO IN SPLINESEGMENTBASE
// 38- IS PLANAR SHOULD RETURN the normal as pair<bool,normal>
// 31- Implement LinearElasticGreensFunction Class
// 12- Finish SimplexMesh library
// 37- NetworkNode, initializations from expansion and contraction
// ALSO: SUBNETWORKS SOULD BE CALLED PARTITIONS
// 38- Finish cleaning STL-Eigen compatibility issue (http://eigen.tuxfamily.org/dox-devel/TopicStlContainers.html) examples include: DislocationSegment::boundaryCollision(),
// IMPLEMENT NEIGHBOR ITERATORS IN NETWORKNODE
// CHANGE CONST VERSION OF EDGEFINDER/VERTEXFINDER PASS this IN MEMBER FUNCTION, REMOVE TEMPLATE SPECIALIZATION AND MAKE STATIC FUNCTIONS
// 9- cellSize should depend on applied load. Or better the number of cell neighbors used in each cell should depend on the applied stress to that cell
// 11- Implement operator << in DislocationNode and DislocationLink to be used by DislocationNetwork::output


// ISSUES
// 3- SplineIntersection::planePlaneType uses the wrong tolerance (10x)
// 10- SplineNetworkBase_common :		assert(iter->second.size()==2);
//     SplineIntersection: //			assert(T0.norm()>tol && "SplineIntersection<3,3,1,2>: T0 too small.");  // NOT RIGHT FOR DIPOLAR LOOPS
//									assert(T1.norm()>tol && "SplineIntersection<3,3,1,2>: T1 too small.");	// NOT RIGHT FOR DIPOLAR LOOPS
// 14- If T0=0 or T1=0, then rl(0)=0/0=NaN !!!!! This is not true since rl still tends to a finite vector. Remove class Parametric curve and implement special case of rl at 0 and 1 for vanishing nodal tangents


#ifndef model_DISLOCATIONNETWORK_H_
#define model_DISLOCATIONNETWORK_H_

//#include <stdlib.h> // random numbers
//#include <time.h> // seed random numbers with time
#include <iostream>
#include <iomanip>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#endif

#define EIGEN_DONT_PARALLELIZE // disable Eigen Internal Parallelization
#include <Eigen/Dense>


#ifdef MPIheaders
#include MPIheaders
#endif

#include <model/Network/Readers/VertexReader.h>
#include <model/Network/Readers/EdgeReader.h>
#include <model/Network/Network.h>

#include <model/Utilities/EigenDataReader.h>
#include <model/Utilities/OutputFile.h>
#include <model/Utilities/SequentialOutputFile.h>
#include <model/Utilities/UniqueOutputFile.h>
#include <model/Utilities/SequentialBinFile.h>

#include <model/DislocationDynamics/DislocationConsts.h>
#include <model/DislocationDynamics/DislocationNetworkTraits.h>
#include <model/DislocationDynamics/DislocationSubNetwork.h>
#include <model/DislocationDynamics/DislocationNode.h>
#include <model/DislocationDynamics/DislocationSegment.h>
#include <model/DislocationDynamics/DislocationSharedObjects.h>
#include <model/DislocationDynamics/GlidePlanes/GlidePlaneObserver.h>
#include <model/DislocationDynamics/Remeshing/DislocationNetworkRemesh.h>
#include <model/DislocationDynamics/Junctions/DislocationJunctionFormation.h>
#include <model/DislocationDynamics/CrossSlip/DislocationCrossSlip.h>
#include <model/DislocationDynamics/Materials/Material.h>

#include <model/DislocationDynamics/IO/DislocationNetworkIO.h>


namespace model {
	
    //	std::string defaultColor    = "\033[0m";	  // the default color for the console
    //	std::string redBoldColor    = "\033[1;31m";   // a bold red color
    //	std::string greenBoldColor  = "\033[1;32m";   // a bold green color
    //	std::string blueBoldColor   = "\033[1;34m";   // a bold blue color
    //	std::string magentaColor    = "\033[0;35m";   // a magenta color
	
	/**************************************************************************/
	/**************************************************************************/
	template <short unsigned int _dim, short unsigned int corder, typename InterpolationType,
	/*	   */ double & alpha, short unsigned int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule>
	class DislocationNetwork : public Network<DislocationNetwork<_dim,corder,InterpolationType,alpha,qOrder,QuadratureRule> >,
	/*                      */ public GlidePlaneObserver<typename TypeTraits<DislocationNetwork<_dim,corder,InterpolationType,alpha,qOrder,QuadratureRule> >::LinkType>{
		
    public:
        
        enum {dim=_dim}; // make dim available outside class
        
		typedef DislocationNetwork<dim,corder,InterpolationType,alpha,qOrder,QuadratureRule> DislocationNetworkType;
		typedef DislocationNetworkType Derived;
#include <model/Network/NetworkTypedefs.h>
		typedef Eigen::Matrix<double,dim,dim>	MatrixDimD;
		typedef Eigen::Matrix<double,dim,1>		VectorDimD;
		typedef typename LinkType::DislocationQuadratureParticleType DislocationQuadratureParticleType;
		typedef DislocationCell<dim,cellSize> SpaceCellType;
		typedef SpaceCellObserver<SpaceCellType,dim,cellSize> SpaceCellObserverType;
		typedef typename SpaceCellObserverType::CellMapType CellMapType;
		typedef GlidePlaneObserver<LinkType> GlidePlaneObserverType;
		enum {NdofXnode=NodeType::NdofXnode};
		
#ifdef UpdateBoundaryConditionsFile
#include UpdateBoundaryConditionsFile
#endif
        
#ifdef DislocationNucleationFile
#include DislocationNucleationFile
//         int nucleationFreq;
#endif
		
#ifdef DislocationNetworkMPI
#include DislocationNetworkMPI
#endif
		
	private:
        
        
        
		
		short unsigned int use_redistribution;
		bool use_junctions;
		
		unsigned int runID;
		
		double Lmax,Lmin,thetaDeg;
		
		bool use_crossSlip;
		double crossSlipDeg;
		double crossSlipLength;
		
		double totalTime;
		
		double dx, dt;
        double vmax;
        
		
		EigenDataReader EDR;
        
		/* readNodes **********************************************************/
		void readNodes(const unsigned int& fileID)
        {/*! Reads file V/V_0.txt and creates DislocationNodes
          */
			typedef VertexReader<'V',11,double> VertexReaderType;
			VertexReaderType  vReader;	// sID,Px,Py,Pz,Tx,Ty,Tz,snID
			assert(vReader.isGood(fileID) && "UNABLE TO READ VERTEX FILE V/V_x (x is the requested file ID).");
			vReader.read(fileID);
			for (VertexReaderType::iterator vIter=vReader.begin();vIter!=vReader.end();++vIter){
				const size_t nodeIDinFile(vIter->first);
				NodeType::set_count(nodeIDinFile);
				const size_t nodeID(this->insert(vIter->second.template segment<NdofXnode>(0).transpose()));
				assert(nodeID==nodeIDinFile);
			}
		}
		
		/* readLinks **********************************************************/
		void readLinks(const unsigned int& fileID)
        {/*! Reads file E/E_0.txt and creates DislocationSegments
          */
			typedef EdgeReader  <'E',11,double>	EdgeReaderType;
			EdgeReaderType    eReader;	// sourceID,sinkID,Bx,By,Bz,Nx,Ny,Nz
			assert(eReader.isGood<true>(fileID) && "Unable to read vertex file E/E_x (x is the requested fileID).");
			eReader.read<true>(fileID);
			for (EdgeReaderType::iterator eIter=eReader.begin();eIter!=eReader.end();++eIter){
				VectorDimD B(eIter->second.template segment<dim>(0  ).transpose()); // Burgers vector
				VectorDimD N(eIter->second.template segment<dim>(dim).transpose()); // Glide plane normal
				const size_t sourceID(eIter->first.first );
				const size_t   sinkID(eIter->first.second);
				assert(this->connect(sourceID,sinkID,B) && "UNABLE TO CREATE CURRENT DISLOCATION SEGMENT.");
			}
		}
        
		/* formJunctions ******************************************************/
		void formJunctions()
        {/*! Performs dislocation junction formation if use_junctions==true
          */
			if (use_junctions)
            {
				double t0=clock();
				std::cout<<"		Forming Junctions: found (";
				DislocationJunctionFormation<DislocationNetworkType>(*this).formJunctions(dx,0.05);
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
			
			for (typename NetworkNodeContainerType::iterator nodeIter=this->nodeBegin();nodeIter!=this->nodeEnd();++nodeIter){
				const double vNorm(nodeIter->second->get_V().norm());
				if (vNorm>vmax){
					vmax=vNorm;
				}
			}
			
			double shearWaveFraction(0.1);
			//short unsigned int shearWaveExp=1;
			if (vmax > Material<Isotropic>::cs*shearWaveFraction){
				dt=dx/vmax;
			}
			else{
                //dt=dx/std::pow(shared.material.cs*shearWaveFraction,shearWaveExp+1)*std::pow(vmax,shearWaveExp);
				//dt=dx/(shared.material.cs*shearWaveFraction)*std::pow(vmax/(shared.material.cs*shearWaveFraction),1);
				dt=dx/(Material<Isotropic>::cs*shearWaveFraction);
			}
			std::cout<<std::setprecision(3)<<std::scientific<<" vmax="<<vmax;
			std::cout<<std::setprecision(3)<<std::scientific<<" dt="<<dt;
			std::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]."<<defaultColor<<std::endl;
		}
		
		
		
		
		/* crossSlip *****************************************************************/
		void crossSlip()
        {/*! Performs dislocation cross-slip if use_crossSlip==true
          */
			if(use_crossSlip){
				double t0=clock();
				std::cout<<"		Performing Cross Slip ... "<<std::flush;
				size_t crossSlipEvents(DislocationCrossSlip<DislocationNetworkType>(*this).crossSlip(crossSlipDeg,crossSlipLength));
				std::cout<<crossSlipEvents<<" cross slip events found ";
				std::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]."<<defaultColor<<std::endl;
			}
		}
		
		
		/* loopInversion *****************************************************************/
		void loopInversion(){
			double t0=clock();
			std::cout<<"		Checking for loop inversions ... "<<std::flush;
			//! 3- Check and remove loop inversions
			std::vector<int> toBeErased;
			for (typename SubNetworkContainerType::iterator snIter=this->ABbegin(); snIter!=this->ABend();++snIter){
				if (snIter->second->loopInversion(dt)){
					std::cout<<"SubNetwork "<<snIter->second->sID<<" containing "<<snIter->second->nodeOrder()<<" is an inverted loop"<<std::endl;
					for (typename SubNetworkNodeContainerType::const_iterator nodeIter=snIter->second->nodeBegin();nodeIter!=snIter->second->nodeEnd();++nodeIter){
						toBeErased.push_back(nodeIter->second->sID);
					}
				}
			}
			std::cout<<" found "<<toBeErased.size()<<" inverted nodes ... ";
			for (unsigned int nn=0;nn<toBeErased.size();++nn){
				this->template remove<true>(toBeErased[nn]);
			}
			std::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]."<<defaultColor<<std::endl;
		}
        
        /***********************************************************/
		void singleStep(const bool& updateUserBC=false){
			//! A simulation step consists of the following:
			std::cout<<blueBoldColor<<"runID="<<runID<<", time="<<totalTime<< ": nodeOrder="<<this->nodeOrder()
            /*                                                           */<< ", linkOrder="<<this->linkOrder()
            /*                                                           */<< ", subNetworks="<<this->Naddresses()
			/*                                                           */<< defaultColor<<std::endl;
            
            
#ifdef DislocationNucleationFile
            if(shared.use_bvp && !(runID%shared.use_bvp)){
                nucleateDislocations(); // needs to be called before updateQuadraturePoints()
                removeBoundarySegments();
            }
#endif
            
			//! 1- Check that all nodes are balanced
			checkBalance();
			
			//! 1 - Update quadrature points
			updateQuadraturePoints();
            
			//! 2- Calculate BVP correction
			if (shared.use_bvp){
				double t0=clock();
				std::cout<<"		Updating bvp stress ... ";
				if(!(runID%shared.use_bvp)){
					shared.domain.update_BVP_Solution(updateUserBC,this);
                    //shared.domain.update_BVP_Solution(updateUserBC,this,*dynamic_cast<const GlidePlaneObserverType*>(this));
                    //#ifdef DislocationNucleationFile
                    //                    nucleateDislocations();
                    //#endif
				}
				for (typename NetworkNodeContainerType::iterator nodeIter=this->nodeBegin();nodeIter!=this->nodeEnd();++nodeIter){ // THIS SHOULD BE PUT IN THE MOVE FUNCTION OF DISLOCATIONNODE
					nodeIter->second->updateBvpStress();
				}
				std::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]."<<defaultColor<<std::endl;
			}
			
			//! 3- Solve the equation of motion
#ifdef DislocationNetworkMPI
			MPIstep();
#endif
			assembleAndSolve();
			
			
			//! 4- Compute time step dt (based on max nodal velocity) and increment totalTime
			make_dt();
			totalTime+=dt;
            
            plasticDistortion+=plasticDistortionRate()*dt;
            
            
            //! 11- Output the current configuration, the corresponding velocities, PK force and dt
            output();
            
            
			
			//! 5- Moves DislocationNodes(s) to their new configuration using stored velocity and dt
			move(dt);
			DislocationNetworkRemesh<DislocationNetworkType>(*this).contract0chordSegments();
			
			//! 6- Moves DislocationNodes(s) to their new configuration using stored velocity and dt
			loopInversion();
			
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
			if (shared.boundary_type==softBoundary){
				double t0=clock();
				std::cout<<"		Removing Segments outside Mesh Boundaries... ";
				typedef bool (LinkType::*link_member_function_pointer_type)(void) const;
				link_member_function_pointer_type boundarySegment_Lmfp;
				boundarySegment_Lmfp=&LinkType::is_boundarySegment;
				this->template disconnect_if<1>(boundarySegment_Lmfp);
				std::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]."<<defaultColor<<std::endl;
			}
		}
        
        
        /***********************************************************/
		void segmentMeshCollision(){
			typedef std::map<double,VectorDimD> MultiIntersectionType; // keeps intersections sorted for increasing parameter
			typedef std::map<LinkIDType,MultiIntersectionType> MultiIntersectionContainerType; // Container of MultiIntersectionType for different links
			MultiIntersectionContainerType multiIntersectionContainer;
			for (typename NetworkLinkContainerType::iterator linkIter=this->linkBegin();linkIter!=this->linkEnd();++linkIter){
				multiIntersectionContainer[linkIter->second->nodeIDPair()]=linkIter->second->boundaryCollision();
			}
			for (typename MultiIntersectionContainerType::const_iterator iter=multiIntersectionContainer.begin();iter!=multiIntersectionContainer.end();++iter){
				this->multiExpand(iter->first.first,iter->first.second,iter->second);
			}
		}
        
		/**********************************************************************************/
		/* PUBLIC SECTION *****************************************************************/
	public:
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
		
		DislocationSharedObjects<LinkType> shared;
		
		//! The number of simulation steps taken by the next call to runByStep()
		int Nsteps;
		double timeWindow;
        //		int outputFrequency;
        //     bool outputGlidePlanes;
        //     bool outputSpaceCells;
        //     bool outputPKforce;
        //     bool outputMeshDisplacement;
        
        MatrixDimD plasticDistortion;
        
		/* Constructor ************************************************************/
        DislocationNetwork(const int& argc, char * const argv[]){
        	if (argc==1){
                // argv[0] is by default the name of the executable so use default name for inputfile
                read("./","DDinput.txt");
            }
            else{
                // argv[1] is assumed to be the filename with working directory the current directory
                read("./",argv[1]);
            }
            
            plasticDistortion.setZero();
        }
        
		
		/* remesh *****************************************************************/
		void remesh(){
			if (use_redistribution){
				if(!(runID%use_redistribution)){
					double t0=clock();
					std::cout<<"		remeshing network... "<<std::flush;
					DislocationNetworkRemesh<DislocationNetworkType>(*this).remesh(Lmax,Lmin,thetaDeg);
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
		void read(const std::string& inputDirectoryName_in, std::string inputFileName){ // TO DO: move this to DislocationNetworkIO.h
			
            std::ostringstream fullName;
            fullName<<inputDirectoryName_in<<inputFileName;
            
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
            EDR.readScalarInFile(fullName.str(),"coreWidthSquared",DislocationQuadratureParticle<dim,cellSize>::a2); // core-width
            assert((DislocationQuadratureParticle<dim,cellSize>::a2)>0.0 && "coreWidthSquared MUST BE > 0.");
            //            EDR.readScalarInFile(fullName.str(),"useMultipoleStress",DislocationQuadratureParticle<dim,cellSize>::useMultipoleStress); // useMultipoleStress
            EDR.readScalarInFile(fullName.str(),"nearCellStressApproximation",DislocationQuadratureParticle<dim,cellSize>::nearCellStressApproximation); // useMultipoleStress
            EDR.readScalarInFile(fullName.str(),"farCellStressApproximation",DislocationQuadratureParticle<dim,cellSize>::farCellStressApproximation); // useMultipoleStress
            assert((DislocationQuadratureParticle<dim,cellSize>::farCellStressApproximation >= DislocationQuadratureParticle<dim,cellSize>::nearCellStressApproximation) && "NEAR FIELD APPROXIMATION IS COARSER THAN FAR FIELD APPROXIMATION");
            
            
            EDR.readMatrixInFile(fullName.str(),"externalStress",shared.externalStress);
			
            
			
			
            EDR.readScalarInFile(fullName.str(),"startAtTimeStep",runID);
			
            EDR.readScalarInFile(fullName.str(),"outputFrequency",DislocationNetworkIO<DislocationNetworkType>::outputFrequency);
            EDR.readScalarInFile(fullName.str(),"outputGlidePlanes",DislocationNetworkIO<DislocationNetworkType>::outputGlidePlanes);
            EDR.readScalarInFile(fullName.str(),"outputSpaceCells",DislocationNetworkIO<DislocationNetworkType>::outputSpaceCells);
            EDR.readScalarInFile(fullName.str(),"outputPKforce",DislocationNetworkIO<DislocationNetworkType>::outputPKforce);
            EDR.readScalarInFile(fullName.str(),"outputMeshDisplacement",DislocationNetworkIO<DislocationNetworkType>::outputMeshDisplacement);
            
            
			
			EDR.readScalarInFile(fullName.str(),"boundary_type",shared.boundary_type);
			if (shared.boundary_type){
				EDR.readScalarInFile(fullName.str(),"use_bvp",shared.use_bvp);
                shared.domain.readMesh();
				if(shared.use_bvp){
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
            
            
            EDR.readScalarInFile(fullName.str(),"use_redistribution",use_redistribution);
            EDR.readScalarInFile(fullName.str(),"Lmin",Lmin);
            assert(Lmin>=0.0);
            assert(Lmin>2.0*dx && "YOU MUST CHOOSE Lmin>2*dx.");
            EDR.readScalarInFile(fullName.str(),"Lmax",Lmax);
            assert(Lmax>Lmin);
            EDR.readScalarInFile(fullName.str(),"thetaDeg",thetaDeg);
            assert(thetaDeg>=0.0);
            assert(thetaDeg<=90.0);
			
            // Cross-Slip inputs
            EDR.readScalarInFile(fullName.str(),"use_crossSlip",use_crossSlip);
            if(use_crossSlip){
                EDR.readScalarInFile(fullName.str(),"crossSlipDeg",crossSlipDeg);
                assert(crossSlipDeg>=0.0 && crossSlipDeg <= 90.0 && "YOU MUST CHOOSE 0.0<= crossSlipDeg <= 90.0");
                EDR.readScalarInFile(fullName.str(),"crossSlipLength",crossSlipLength);
                assert(crossSlipLength>=Lmin && "YOU MUST CHOOSE crossSlipLength>=Lmin.");
				//				assert(crossSlipLength<Lmin && "YOU MUST CHOOSE crossSlipLength<Lmin."); // Because otherwise cross-slip points would go outside segment
            }
			
            // Read Vertex and Edge information
            readNodes(runID);
            readLinks(runID);
            
            
            
            if (shared.use_bvp && (shared.boundary_type==softBoundary)) { // MOVE THIS WITH REST OB BVP STUFF
                //                shared.vbsc.read(runID,&shared);
                shared.vbsc.initializeVirtualSegments(*this);
            }
            
            
//#ifdef DislocationNucleationFile
//            EDR.readScalarInFile(fullName.str(),"nucleationFreq",nucleationFreq);
//#endif

            
			
			// Initializing initial configuration
			std::cout<<redBoldColor<<"runID "<<runID<<" (initial configuration). nodeOrder="<<this->nodeOrder()<<", linkOrder="<<this->linkOrder()<<defaultColor<<std::endl;
			move(0.0);	// initial configuration
			output();	// initial configuration, this overwrites the input file
			if (runID==0){ // not a restart
				remesh();	// expand initial FR sources
			}
			updateQuadraturePoints();
			++runID;     // increment the runID counter
		}
		
		/* solve ****************************************************************/
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
            for (unsigned int k=0;k<this->Naddresses();++k){
                typename SubNetworkContainerType::iterator snIter(this->ABbegin()); //  data within a parallel region is private to each thread
                std::advance(snIter,k);
                //				snIter->second->solve();
                if (snIter->second->nodeOrder()>=shared.minSNorderForSolve){
                    snIter->second->sparseSolve();
                }
            }
#else
			for (typename SubNetworkContainerType::iterator snIter=this->ABbegin(); snIter!=this->ABend();++snIter){
                //				snIter->second->solve();
                //snIter->second->sparseSolve();
                if (snIter->second->nodeOrder()>=shared.minSNorderForSolve){
                    snIter->second->sparseSolve();
                }
			}
#endif
			std::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]."<<defaultColor<<std::endl;
		}
		
		/* output ***************************************************************/
		void output() const
        {/*! Outputs DislocationNetwork data
          */
            double t0=clock();
            if (!(runID%DislocationNetworkIO<DislocationNetworkType>::outputFrequency)){
#ifdef DislocationNetworkMPI
				if(localProc==0){
                    DislocationNetworkIO<DislocationNetworkType>(*this).outputTXT(runID);
				}
#else
                DislocationNetworkIO<DislocationNetworkType>(*this).outputTXT(runID);
#endif
			}
			std::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]."<<defaultColor<<std::endl;
		}
		
		
		/* updateQuadraturePoints ********************************/
		void updateQuadraturePoints()
        {
			double t0=clock();
			std::cout<<"		Updating Quadrature Points... "<<std::flush;
			for (typename NetworkLinkContainerType::iterator linkIter=this->linkBegin();linkIter!=this->linkEnd();++linkIter){
				linkIter->second->quadratureParticleVector.clear(); // first clear all quadrature points for all segments
			}
            for (typename NetworkLinkContainerType::iterator linkIter=this->linkBegin();linkIter!=this->linkEnd();++linkIter){
				Quadrature<1,qOrder>::execute(linkIter->second,&LinkType::updateQuadGeometryKernel); // then update again
			}
            for (typename CellMapType::const_iterator cellIter=SpaceCellObserverType::begin();cellIter!=SpaceCellObserverType::end();++cellIter){
                cellIter->second->computeCenterStress();
            }
			std::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]."<<defaultColor<<std::endl;
		}
		
		
		/***********************************************************/
		void move(const double & dt_in)
        {/*! Moves all nodes in the DislocationNetwork using the stored velocity and current dt
          */
			std::cout<<"		Moving Dislocation Nodes... "<<std::flush;
			double t0=clock();
            //			typedef void (NodeType::*NodeMemberFunctionPointerType)(const double&); // define type of Link member function
            //			NodeMemberFunctionPointerType Nmfp(&NodeType::move); // Lmfp is a member function pointer to Link::assemble
            //			double t0=clock();
            //			this->parallelExecute(Nmfp,dt_in);
            
			for (typename NetworkNodeContainerType::iterator nodeIter=this->nodeBegin();nodeIter!=this->nodeEnd();++nodeIter){
				nodeIter->second->move(dt_in);
			}
			std::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]."<<defaultColor<<std::endl;
		}
		
		/***********************************************************/
		void runByStep(const bool& updateUserBC=false)
        {/*! Runs Nsteps simulation steps
          */
			double ts=clock();
			for (int k=0;k<Nsteps;++k){
				std::cout<<std::endl; // leave a blank line
				std::cout<<blueBoldColor<<"Step "<<k+1<<" of "<<Nsteps<<defaultColor<<std::endl;
				singleStep(updateUserBC);
			}
			std::cout<<greenBoldColor<<std::setprecision(3)<<std::scientific<<Nsteps<< " simulation steps completed in "<<(clock()-ts)/CLOCKS_PER_SEC<<" [sec]"<<defaultColor<<std::endl;
		}
		
		
		/***********************************************************/
		void runByTime(const bool& updateUserBC=false)
        {/*! Runs a number simulation steps corresponding to a total
          * dimensionless time timeWindow
          */
			double ts=clock();
			double elapsedTime(0.0);
			while (elapsedTime<timeWindow){
				std::cout<<std::endl; // leave a blank line
				std::cout<<blueBoldColor<<"Time "<<elapsedTime<<" of "<<timeWindow<<defaultColor<<std::endl;
				singleStep(updateUserBC);
				elapsedTime+=dt;
			}
			std::cout<<greenBoldColor<<std::setprecision(3)<<std::scientific<<timeWindow<< " simulation time completed in "<<(clock()-ts)/CLOCKS_PER_SEC<<" [sec]"<<defaultColor<<std::endl;
		}
		
		
		/**********************************************************************/
		void checkBalance() const
        {/*! Checks that each DislocationNode is balanced, and asserts otherwise
          * Exceptions are nodes with only one neighbors (FR source)
          * and nodes on the boundary.
          */
			for (typename NetworkNodeContainerType::const_iterator nodeIter=this->nodeBegin();nodeIter!=this->nodeEnd();++nodeIter){
				if (nodeIter->second->neighborhood().size()>2){
                    const bool nodeIsBalanced(nodeIter->second->is_balanced());
                    if (!nodeIsBalanced && nodeIter->second->nodeMeshLocation==insideMesh) {
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
			double temp=0.0;
			for (typename NetworkLinkContainerType::iterator linkIter=this->linkBegin();linkIter!=this->linkEnd();++linkIter){
				temp+= linkIter->second->energy();
			}
			return temp;
		}
		
		/**********************************************************************/
		template <bool useFullField=true>
		MatrixDimD stress(const VectorDimD & Rfield) const {
			MatrixDimD temp=MatrixDimD::Zero();
			if(useFullField){
				for (typename NetworkLinkContainerType::const_iterator linkIter=this->linkBegin();linkIter!=this->linkEnd();++linkIter){
					temp+= linkIter->second->stress_source(Rfield);
				}
			}
			else{
                assert(0 && "RE-ENABLE THIS WITH NEW CELL CLASS");
                //				SpaceCellObserverType sCO;
                //				typename SpaceCellObserverType::SharedPtrType sharedCellPointer(sCO.getCellByPosition(Rfield));
                //				for(typename SpaceCellType::ParticleContainerType::const_iterator particleIter =sharedCellPointer->neighborParticleContainer.begin();
                //					/*                                                         */ particleIter!=sharedCellPointer->neighborParticleContainer.end();
                //					/*                                                         */ ++particleIter){
                //					temp+=(*particleIter)->stress_at(Rfield);
                //				}
			}
            if (shared.use_bvp && (shared.boundary_type==1)){
                temp+= shared.vbsc.stress(Rfield);
            }
			return temp;
		}
        
		/********************************************************/
		MatrixDimD stressFromGlidePlane(const Eigen::Matrix<double,dim+1,1>& key, const VectorDimD& Rfield) const {
			//GlidePlaneObserver<dim,LinkType> gpObsever;
			MatrixDimD temp(MatrixDimD::Zero());
			typedef typename GlidePlaneObserver<LinkType>::GlidePlaneType GlidePlaneType;
			std::pair<bool, const GlidePlaneType* const> isGp(this->isGlidePlane(key));
			if (isGp.first){
				temp=isGp.second->stress(Rfield);
			}
            if (shared.use_bvp && (shared.boundary_type==1)) {
                temp+= shared.vbsc.stressFromGlidePlane(key,Rfield);
            }
			return temp;
		}
		
		/********************************************************/
		VectorDimD displacement(const VectorDimD & Rfield,const VectorDimD & S) const {
			VectorDimD temp(VectorDimD::Zero());
			for (typename NetworkLinkContainerType::const_iterator linkIter=this->linkBegin();linkIter!=this->linkEnd();++linkIter){
				temp+= linkIter->second->displacement(Rfield,S);
			}
			return temp;
		}
		
        
        /********************************************************/
		MatrixDimD plasticDistortionRate() const {
			MatrixDimD temp(MatrixDimD::Zero());
			for (typename NetworkLinkContainerType::const_iterator linkIter=this->linkBegin();linkIter!=this->linkEnd();++linkIter){
				temp+= linkIter->second->plasticDistortionRate();
			}
			return temp;
		}
        
        
        /********************************************************/
		MatrixDimD plasticStrainRate() const {
            //			MatrixDimD temp(MatrixDimD::Zero());
            //			for (typename NetworkLinkContainerType::const_iterator linkIter=this->linkBegin();linkIter!=this->linkEnd();++linkIter){
            //				temp+= linkIter->second->plasticStrainRate();
            //			}
            const MatrixDimD temp(plasticDistortionRate());
			return (temp+temp.transpose())*0.5;
		}
		
		/********************************************************/
		MatrixDimD lattice_rotation(const VectorDimD & Rfield) const {
			MatrixDimD temp(MatrixDimD::Zero());
			for (typename NetworkLinkContainerType::const_iterator linkIter=this->linkBegin();linkIter!=this->linkEnd();++linkIter){
				temp+= linkIter->second->lattice_rotation_source(Rfield);
			}
			return temp;
		}
		
		/********************************************************/
		MatrixDimD elasticDistortion(const VectorDimD & Rfield) const {
			MatrixDimD temp(MatrixDimD::Zero());
			for (typename NetworkLinkContainerType::const_iterator linkIter=this->linkBegin();linkIter!=this->linkEnd();++linkIter){
				temp+= linkIter->second->displacement_gradient_source(Rfield);
			}
			return temp;
		}
		
		/********************************************************/
		// function to return the length of the whole dislocation network
		double network_length() const {
			double temp(0.0);
			for (typename NetworkLinkContainerType::const_iterator linkIter=this->linkBegin();linkIter!=this->linkEnd();++linkIter){
				temp+= linkIter->second->template arcLength<qOrder,QuadratureRule>();
			}
			return temp;
		}
        
        /********************************************************/
        const unsigned int& runningID() const {
            return runID;
        }
        
        /********************************************************/
        const double& vMax() const
        {
            return vmax;
        }
		
	};
    
    
	//////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////
} // namespace model
#endif


//		/***********************************************************/
//		void checkPlanarity(){
//			for (typename NetworkLinkContainerType::const_iterator linkIter=this->linkBegin();linkIter!=this->linkEnd();++linkIter){
//				std::cout<<"checking planarity for link "<<linkIter->second->nodeIDPair.first<<"->"<<linkIter->second->nodeIDPair.second<<std::endl;
//				Eigen::Matrix<double,dim,4> H2(linkIter->second->hermiteCoefficients());
//				Eigen::Matrix<double,dim,1> n2(linkIter->second->glidePlaneNormal);
//				const double absN2dotC2(std::fabs(n2.dot(H2.col(2)-H2.col(0))));
//				const double absN2dotT2source(std::fabs(n2.dot(H2.col(1))));
//				const double   absN2dotT2sink(std::fabs(n2.dot(H2.col(3))));
//				if (absN2dotC2>FLT_EPSILON || absN2dotT2source>FLT_EPSILON || absN2dotT2sink>FLT_EPSILON){
//					std::cout<<"n2="<<n2.transpose()<<std::endl;
//					std::cout<<"H2="<<H2<<std::endl;
//					std::cout<<"absN2dotC2="<<absN2dotC2<<"\n";
//					std::cout<<"absN2dotT2source="<<absN2dotT2source<<"\n";
//					std::cout<<"absN2dotT2sink="<<absN2dotT2sink<<"\n";
//					//							std::cout<<"n2"std::fabs(n2.dot(H2.col(3)))<<std::endl;
//				}
//
//				assert(      absN2dotC2<FLT_EPSILON);
//				assert(absN2dotT2source<FLT_EPSILON);
//				assert(  absN2dotT2sink<FLT_EPSILON);
//			}
//		}
