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



// DONE 15- STORE BVP DISPACEMENTS IN /D
// DONE Clean GramSchmidt, derive from std::vector with right allocator for Eigen
// 1- DONEBecause Output is after redistribution, Cell Files are not correct because GaussPoints are destroyed when segments are destroyed
// DONE: CROSS_SLIP: NO, THIS CAUSES THE SAME NORMAL TO BE USED
// DONE: Strategy changed with new BVP. 17- add inifinite line stress fields for segments terminating on the surface
// DONE: STRATEGY CHANGED WITH NEW BVP. -1 - code isBoundarySubnetwork. If true don't solve AND set velocity to 0.

// BEING MODIFIED
// 1- BIN WRITER/READER. DEFINE DISLOCATIONSEGMENT::OUTDATA as Eigen::Matrix<double,1,9>. Then in Network add the function friend void operator<< (SequentialBinFile<template Char,???>, ...)
// 2- Zeo node Tangent Problem
// BVP: integration method
// 3 crossSlip
// Motion of BoundarySubnetworks

// TO DO
// -2: updateQuadraturePoints should be called as part of the topologyChangeActions() and move()

// 0 - Finish DepthFirst class. Don't allow to search/execute if N=0, so that N=1 means that node only, n=2 means first neighbor. Or change SpineNodeBase_CatmullRom::TopologyChangeActions
// 1- Implement operator << in DislocationNode, GlidePlane, SpaceCell
// 2- remove AddressBook, wherever possible, Done in Node chain
// 40 - clean MultiExpand in Network Layer, remove get_r from there
// 39 - CLEAN NEW CATMUll-ROM (ENERGY). REMOVE CRNEGHBORS AND USE NEIGHBORHOOD ITSELF
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
// 12 PlanarSplineImplicitization::intersect with has problems with some ill-conditioned matrices and gives wrong results.
// 14- If T0=0 or T1=0, then rl(0)=0/0=NaN !!!!! This is not true since rl still tends to a finite vector. Remove class Parametric curve and implement special case of rl at 0 and 1 for vanishing nodal tangents


#ifndef model_DISLOCATIONNETWORK_H_
#define model_DISLOCATIONNETWORK_H_


#ifndef VERBOSELEVEL
#define VERBOSELEVEL 3
#endif

#include <iostream>
#include <iomanip>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#endif

#include <Eigen/Dense>


#ifdef MPIheaders
#include MPIheaders
#endif


#include <model/Network/Readers/VertexReader.h>
#include <model/Network/Readers/EdgeReader.h>
#include <model/Network/Network.h>

#include <model/Dislocations/DislocationConsts.h>
#include <model/Dislocations/DislocationNetworkTraits.h>
#include <model/Dislocations/DislocationSubNetwork.h>
#include <model/Dislocations/DislocationNode.h>
#include <model/Dislocations/DislocationSegment.h>
#include <model/Utilities/EigenDataReader.h>
#include <model/Utilities/OutputFile.h>
#include <model/Utilities/SequentialOutputFile.h>
#include <model/Utilities/UniqueOutputFile.h>
#include <model/Dislocations/DislocationSharedObjects.h>
//#include <model/Dislocations/DislocationQuadratureParticle.h>
#include <model/Dislocations/GlidePlanes/GlidePlaneObserver.h>
#include <model/Dislocations/Remeshing/DislocationNetworkRemesh.h>
#include <model/Dislocations/Junctions/DislocationJunctionFormation.h>
#include <model/Dislocations/CrossSlip/DislocationCrossSlip.h>

#include <model/Utilities/SequentialBinFile.h>




namespace model {
	
	std::string defaultColor    = "\033[0m";	    // the default color for the console
	std::string redBoldColor    = "\033[1;31m";   // a bold red color
	std::string greenBoldColor  = "\033[1;32m";   // a bold green color
	std::string blueBoldColor   = "\033[1;34m";   // a bold blue color
	std::string magentaColor    = "\033[0;35m";   // a magenta color
	
	
	//////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////
	template <short unsigned int dim, short unsigned int corder, typename InterpolationType,
	/*	   */ double & alpha, short unsigned int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule, 
	/*	   */ typename MaterialType>
	class DislocationNetwork : public Network<DislocationNetwork<dim,corder,InterpolationType,alpha,qOrder,QuadratureRule,MaterialType> >,
	/*                      */ public GlidePlaneObserver<dim,typename TypeTraits<DislocationNetwork<dim,corder,InterpolationType,alpha,qOrder,QuadratureRule,MaterialType> >::LinkType> {
		
public:

		typedef DislocationNetwork<dim,corder,InterpolationType,alpha,qOrder,QuadratureRule,MaterialType> DislocationNetworkType;
		typedef DislocationNetworkType Derived;
#include <model/Network/NetworkTypedefs.h>
		typedef Eigen::Matrix<double,dim,dim>	MatrixDimD;
		typedef Eigen::Matrix<double,dim,1>		VectorDimD;
		typedef typename LinkType::DislocationQuadratureParticleType DislocationQuadratureParticleType;
		typedef SpaceCell<DislocationQuadratureParticleType,dim,cellSize> SpaceCellType;
		typedef SpaceCellObserver<SpaceCellType,dim,cellSize> SpaceCellObserverType;
		typedef typename SpaceCellObserverType::CellMapType CellMapType;
		typedef GlidePlaneObserver<dim,LinkType> GlidePlaneObserverType;
		enum {NdofXnode=NodeType::NdofXnode};
		


#ifdef UpdateBoundaryConditionsFile
#include UpdateBoundaryConditionsFile
#endif
		
		
#ifdef DislocationNetworkMPI
#include DislocationNetworkMPI	
#endif
		
	private:
		
		short unsigned int use_redistribution; 
		short unsigned int use_junctions;	// 0= no junctions, 1=only annahilation, 2= all junction types
		
		unsigned int runID;
		
		double Lmax,Lmin,thetaDeg;
		
		bool use_crossSlip;
		double crossSlipDeg;
		double crossSlipLength;
		
		double totalTime;
		
		double dx, dt;
		
		
		EigenDataReader EDR;
		
		
		/* readNodes ***************************************************************/
		void readNodes(const unsigned int& fileID/*, std::map<size_t,size_t>& nodeIDrestartMap*/){
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
		
		/* readLinks ***************************************************************/
		void readLinks(const unsigned int& fileID/*, const std::map<size_t,size_t>& nodeIDrestartMap*/){	
			typedef EdgeReader  <'E',11,double>	EdgeReaderType;
			EdgeReaderType    eReader;	// sourceID,sinkID,Bx,By,Bz,Nx,Ny,Nz
			assert(eReader.isGood<true>(fileID) && "Unable to read vertex file E/E_x where x is the requested fileID.");
			eReader.read<true>(fileID);			
			for (EdgeReaderType::iterator eIter=eReader.begin();eIter!=eReader.end();++eIter){
				VectorDimD B(eIter->second.template segment<dim>(0  ).transpose()); // Burgers vector
				VectorDimD N(eIter->second.template segment<dim>(dim).transpose()); // Glide plane normal
				const size_t sourceID(eIter->first.first );
				const size_t   sinkID(eIter->first.second);
				assert(this->connect(sourceID,sinkID,B) && "UNABLE TO CREATE CURRENT DISLOCATION SEGMENT.");
			}
		}
		
		/* formJunctions ************************************************************/
		void formJunctions(){
			if (use_junctions){
				double t0=clock();
				std::cout<<"		Forming Junctions: found (";
				DislocationJunctionFormation<DislocationNetworkType>(*this).formJunctions(dx,0.05);
				std::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]."<<defaultColor<<std::endl;			
			}
		}
		
		/* make_dt *****************************************************************/
		void make_dt(){
			/*! Computes the time step size \f$dt\f$ for the current simulation step, 
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
			
			double vmax(0.0);
			
			for (typename NetworkNodeContainerType::iterator nodeIter=this->nodeBegin();nodeIter!=this->nodeEnd();++nodeIter){
				const double vNorm(nodeIter->second->get_V().norm());
				if (vNorm>vmax){
					vmax=vNorm;
				}
			}
			
			double shearWaveFraction(0.1);
			//short unsigned int shearWaveExp=1;
			if (vmax > shared.material.cs*shearWaveFraction){
				dt=dx/vmax;
			}
			else{
				//	dt=dx/std::pow(shared.material.cs*shearWaveFraction,shearWaveExp+1)*std::pow(vmax,shearWaveExp);
				dt=dx/(shared.material.cs*shearWaveFraction)*std::pow(vmax/(shared.material.cs*shearWaveFraction),1);
			}
			std::cout<<std::setprecision(3)<<std::scientific<<" vmax="<<vmax;
			std::cout<<std::setprecision(3)<<std::scientific<<" dt="<<dt;
			std::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]."<<defaultColor<<std::endl;			
		}
		
		
		
		
		/* crossSlip *****************************************************************/
		void crossSlip(){
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
		
		/**********************************************************************************/
		/* PUBLIC SECTION *****************************************************************/
	public:
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
		
		DislocationSharedObjects<dim,MaterialType> shared;
		
		//! The number of simulation steps taken by the next call to runByStep()
		int Nsteps;
		double timeWindow;
		int outputFrequency;
		
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
		
		/* get_dt *****************************************************************/
		const double& get_dt() const {
			return dt;
		}
		
		const double& get_totalTime() const {
			return totalTime;
		}
		
		
		
		/* read *****************************************************************/
		void read(const std::string& inputDirectoryName_in, std::string inputFileName){
			
			std::ostringstream fullName;
			fullName<<inputDirectoryName_in<<inputFileName;
			
			
			EDR.readMatrixInFile(fullName.str(),"initialStress",shared.externalStress);
			
			
			MatrixDimD C2Gtemp;
			EDR.readMatrixInFile(fullName.str(),"C2G",C2Gtemp);
			shared.material.rotate(C2Gtemp);
			
			EDR.readScalarInFile(fullName.str(),"startAtTimeStep",runID);
			
			EDR.readScalarInFile(fullName.str(),"outputFrequency",outputFrequency);
			
			EDR.readScalarInFile(fullName.str(),"boundary_type",shared.boundary_type);
			if (shared.boundary_type){
				EDR.readScalarInFile(fullName.str(),"use_bvp",shared.use_bvp);
				shared.domain.readMesh();
				if(shared.use_bvp){
					shared.domain.readInputBCs();
					//shared.domain.setBoundaryConditions();		// temporary
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
			
			EDR.readScalarInFile(fullName.str(),"use_redistribution",use_redistribution);
			EDR.readScalarInFile(fullName.str(),"use_junctions",use_junctions);
			
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
			
			// Initializing initial configuration 
			std::cout<<redBoldColor<<"runID "<<runID<<" (initial configuration). nodeOrder="<<this->nodeOrder()<<", linkOrder="<<this->linkOrder()<<defaultColor<<std::endl;
			move(0.0);	// initial configuration
			output();	// initial configuration, this overwrites the input file
			if (runID==0){ // not a restart
				remesh();	// expand initial FR sources
			}
			updateQuadraturePoints();			
			++runID;     // increment the runID counter first, since the original configuration has been outputed already
		}
		
		/* solve ****************************************************************/
		void solve(){
			//! 1- Loop over DislocationSegments and assemble stiffness matrix and force vector
			std::cout<<"		Assembling..."<<std::flush;
			typedef void (LinkType::*LinkMemberFunctionPointerType)(void); // define type of Link member function
			LinkMemberFunctionPointerType Lmfp(&LinkType::assemble); // Lmfp is a member function pointer to Link::assemble
			double t0=clock();
			this->parallelExecute(Lmfp);			
			std::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]."<<defaultColor<<std::endl;			

			//! 2- Loop over DislocationSubNetworks, assemble subnetwork stiffness matrix and force vector, and solve
			std::cout<<"		Solving..."<<std::flush;
			t0=clock();
#ifdef _OPENMP
#pragma omp parallel
#pragma omp single
{
#endif
			for (typename SubNetworkContainerType::iterator snIter=this->ABbegin(); snIter!=this->ABend();++snIter){ 
#ifdef _OPENMP
#pragma omp task firstprivate(snIter)	
#endif					
				snIter->second->solve();
			}
#ifdef _OPENMP
#pragma omp taskwait
}
#endif			
			std::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]."<<defaultColor<<std::endl;			
		}
		
		/* output ***************************************************************/
		void output() const {
			std::cout<<"		Outputing to";
			double t0=clock();
			
			//! 1- Outputs the Edge informations to file E_*.txt where * is the current simulation step
			SequentialOutputFile<'E',1>::set_increment(outputFrequency); // Edges_file;
			SequentialOutputFile<'E',1>::set_count(runID); // Edges_file;
			SequentialOutputFile<'E',1> Edges_file;
			Edges_file << *dynamic_cast<const NetworkLinkContainerType*>(this);
			std::cout<<" E/E_"<<Edges_file.sID;
			
			typedef std::pair<std::pair<int,int>,Eigen::Matrix<double,1,9> > BinEdgeType;
			SequentialBinFile<'E',BinEdgeType>::set_increment(outputFrequency);
			SequentialBinFile<'E',BinEdgeType>::set_count(runID);
			SequentialBinFile<'E',BinEdgeType> binEdgeFile;
			for (typename NetworkLinkContainerType::const_iterator linkIter=this->linkBegin();linkIter!=this->linkEnd();++linkIter){
				Eigen::Matrix<double,1,9> temp;
				temp<< linkIter->second->flow.transpose(),
				/*  */ linkIter->second->glidePlaneNormal.transpose(),
				/*  */ linkIter->second->sourceTfactor,
				/*  */ linkIter->second->sinkTfactor,
				/*  */ linkIter->second->pSN()->sID;
				binEdgeFile.write(std::make_pair(linkIter->first,temp));
			}
			
			//! 2- Outputs the Vertex informations to file V_*.txt where * is the current simulation step
			SequentialOutputFile<'V',1>::set_increment(outputFrequency); // Vertices_file;
			SequentialOutputFile<'V',1>::set_count(runID); // Vertices_file;
			SequentialOutputFile<'V',1> Vertices_file;
			for (typename NetworkNodeContainerType::const_iterator nodeIter=this->nodeBegin();nodeIter!=this->nodeEnd();++nodeIter){				
				Vertices_file << nodeIter->second->sID<<"\t"
				/*         */ << std::setprecision(15)<<std::scientific<<nodeIter->second->get_P().transpose()<<"\t"
				/*         */ << std::setprecision(15)<<std::scientific<<nodeIter->second->get_T().transpose()<<"\t"
				/*         */ << nodeIter->second->pSN()->sID<<"\t";
				if (shared.use_bvp){ //output in deformed configuration
					Vertices_file << std::setprecision(15)<<std::scientific<<nodeIter->second->deformedPosition().transpose()<<"\t";
				}
				else{
					Vertices_file<< VectorDimD::Zero().transpose();
				}
				Vertices_file << std::endl;
				
			}
			std::cout<<", V/V_"<<Vertices_file.sID;
			
			//! 3- Outputs the nearest neighbor Cell structures to file C_*.txt where * is the current simulation step
			SequentialOutputFile<'C',1>::set_increment(outputFrequency); // Cell_file;
			SequentialOutputFile<'C',1>::set_count(runID); // Cell_file;
			SequentialOutputFile<'C',1> Cell_file;
			SpaceCellObserverType SPC;
			int cID(0);
			for (typename CellMapType::const_iterator cellIter=SPC.begin();cellIter!=SPC.end();++cellIter){
				Cell_file<<cID<<"\t"<<cellIter->second->cellID.transpose()<<"\t"<<cellSize<<std::endl;
				++cID;
			}
			std::cout<<", C/C_"<<Cell_file.sID;
			
			//! 4- Outputs the glide planes
			SequentialOutputFile<'G',1>::set_increment(outputFrequency); // GlidePlanes_file;
			SequentialOutputFile<'G',1>::set_count(runID); // GlidePlanes_file;
			SequentialOutputFile<'G',1> glide_file;
			glide_file << *dynamic_cast<const GlidePlaneObserverType*>(this); 
			std::cout<<", G/G_"<<glide_file.sID;
			
#ifdef customUserOutputs
#include customUserOutputs
#endif
			
			std::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]."<<defaultColor<<std::endl;
		}
		
		
		/* updateQuadraturePoints ********************************/
		void updateQuadraturePoints(){
			double t0=clock();
			std::cout<<"		Updating Quadrature Points... "<<std::flush;
			for (typename NetworkLinkContainerType::iterator linkIter=this->linkBegin();linkIter!=this->linkEnd();++linkIter){
				linkIter->second->quadratureParticleVector.clear();
				Quadrature<1,qOrder>::execute(linkIter->second,&LinkType::updateQuadGeometryKernel);
			}
			std::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]."<<defaultColor<<std::endl;						
		}
		
		
		/***********************************************************/
		void move(const double & dt_in){
			//! 1- Moves all nodes in the DislocationNetwork using the stored velocity and current dt
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
		
		
		/***********************************************************/
		void runByStep(const bool& updateUserBC=false){
			/*! A simulation step consists of the following:
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
		void runByTime(const bool& updateUserBC=false){
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
		
		
		void checkBalance() const {
			for (typename NetworkNodeContainerType::const_iterator nodeIter=this->nodeBegin();nodeIter!=this->nodeEnd();++nodeIter){
				if (nodeIter->second->neighborhood().size()>2){
					assert(nodeIter->second->is_balanced() && "NODE IS NOT BALANCED");
				}
				//->move(dt_in);
			}
		}
		
		
		/***********************************************************/
		void singleStep(const bool& updateUserBC=false){
			std::cout<<blueBoldColor<<"runID="<<runID<<", time="<<totalTime<< ": nodeOrder="<<this->nodeOrder()<<", linkOrder="<<this->linkOrder()<<defaultColor<<std::endl;
			//! A simulation step consists of the following:
			
			checkBalance();
			
			//! 1 - Update quadrature points
			updateQuadraturePoints();
			
			//! 2- Calculate BVP correction 				
			if (shared.use_bvp){
				
				double t0=clock();
				std::cout<<"		Updating bvp stress ... ";
				
				if(!(runID%shared.use_bvp)){
					shared.domain.update_BVP_Solution(updateUserBC,this);
					//					shared.domain.setBoundaryConditions();		// temporary					
					//#ifdef UpdateBoundaryConditionsFile
					//					updateBVP_BCs(k);        	// update BCs for the BVP 
					//#endif						
					//					shared.domain.solveBVP(true,this);
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
			solve();
			
			
			//! 4- Compute time step dt (based on max nodal velocity) and increment totalTime
			make_dt();	
			totalTime+=dt;
			
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
			updateQuadraturePoints();
//			updateQuadraturePoints();
			
			
			//! 11- Output the current configuration, the corresponding velocities and dt
			if (!(runID%outputFrequency)){
				
#ifdef DislocationNetworkMPI
				if(localProc==0){
					output();
				}
#else
				output();
#endif
			}
			
			//! 12 - Increment runID counter
			++runID;     // increment the runID counter first, since the original configuration has been outputed already
		}
		
		
		/* remeshByContraction **************************************/
		void removeBoundarySegments(){
			if (shared.boundary_type==softBoundary){
				double t0=clock();
				std::cout<<"		Removing Segments outside Mesh Boundaries... ";
				for (typename NetworkLinkContainerType::iterator linkIter=this->linkBegin();linkIter!=this->linkEnd();++linkIter){
					if(linkIter->second->is_boundarySegment()){
					//	shared.domain.vbsc.push_back(*(linkIter->second));
					}
				}


				typedef bool (LinkType::*link_member_function_pointer_type)(void) const; 
				link_member_function_pointer_type boundarySegment_Lmfp;
				boundarySegment_Lmfp=&LinkType::is_boundarySegment;
				this->template disconnect_if<1>(boundarySegment_Lmfp);
				std::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]."<<defaultColor<<std::endl;			
			}
		}
		
		/********************************************************/
		// energy
		double energy(){
			double temp=0.0;
			for (typename NetworkLinkContainerType::iterator linkIter=this->linkBegin();linkIter!=this->linkEnd();++linkIter){
				temp+= linkIter->second->energy(); 
			}
			return temp;
		}
		
		/********************************************************/
		template <bool useFullField=true>
		MatrixDimD stress(const VectorDimD & Rfield) const {
			MatrixDimD temp=MatrixDimD::Zero();
			if(useFullField){
				for (typename NetworkLinkContainerType::const_iterator linkIter=this->linkBegin();linkIter!=this->linkEnd();++linkIter){
					temp+= linkIter->second->stress_source(Rfield); 
				}				
			}
			else{
				SpaceCellObserverType sCO;
				typename SpaceCellObserverType::SharedPtrType sharedCellPointer(sCO.getCellByPosition(Rfield));
				for(typename SpaceCellType::ParticleContainerType::const_iterator particleIter =sharedCellPointer->neighborParticleContainer.begin();
					/*                                                         */ particleIter!=sharedCellPointer->neighborParticleContainer.end();
					/*                                                         */ ++particleIter){
					temp+=(*particleIter)->stress_at(Rfield);
				}
			}
			return temp;
		}

		/********************************************************/
		MatrixDimD stressFromGlidePlane(const Eigen::Matrix<double,dim+1,1>& key, const VectorDimD& Rfield) const {
			GlidePlaneObserver<dim,LinkType> gpObsever;
			MatrixDimD temp(MatrixDimD::Zero());
			typedef typename GlidePlaneObserver<dim,LinkType>::GlidePlaneType GlidePlaneType;
			std::pair<bool, const GlidePlaneType* const> isGp(gpObsever.isGlidePlane(key));		
			if (isGp.first){
				temp=isGp.second->stress(Rfield);
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
		MatrixDimD plasticStrainRate() const {			
			MatrixDimD temp(MatrixDimD::Zero());
			for (typename NetworkLinkContainerType::const_iterator linkIter=this->linkBegin();linkIter!=this->linkEnd();++linkIter){
				temp+= linkIter->second->plasticStrainRate(); 
			}
			return temp;
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
		
	};
	
	//////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////
} // namespace model
#endif


//			GlidePlaneObserverType gpObsever; // MOVE THIS TO DATA MEMBER
//			GlidePlaneObserver<dim,LinkType> gpObsever;
//			glide_file << gpObsever;



