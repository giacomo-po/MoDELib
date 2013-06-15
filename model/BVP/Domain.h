/* This file is part of finite element solution of BVP attached with model "Mechanics of Defects Evolution  Library".
 *
 * Copyright (C) 2011 by Mamdouh Mohamed <mamdouh.s.mohamed@gmail.com>,
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_domain_H_
#define model_domain_H_

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#define omp_get_num_threads() 1
#endif


#include <iterator>

#include <Eigen/Dense>
#include <Eigen/Sparse>

//#include <Eigen/Sparse>
//#include <model/BVP/UmfPackSupport.h>
//#include "model/BVP/SuperLUSupport.h"

#include <map>

#include "model/Quadrature/Quadrature.h"
#include "model/Utilities/SequentialOutputFile.h"
#include <model/Utilities/UniqueOutputFile.h>
#include <model/DislocationDynamics/GlidePlanes/GlidePlaneObserver.h>

#include "model/BVP/Node.h"
#include "model/BVP/Tetrahedron.h"
#include "model/BVP/Triangle.h"
#include "model/BVP/Face.h"
#include "model/BVP/SearchData.h"


#include <stdio.h>


#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/ptr_container/ptr_map.hpp>


namespace model{
	
    //	template <typename MaterialType>
	class Domain {
		
		enum {dim=3};
        
#include <model/BVP/commonTypeDefs.h>
		
	public:
		typedef boost::ptr_vector<bvpfe::Node<dim> >  NodeContainerType;
		typedef boost::ptr_vector<bvpfe::Tetrahedron> TetContainerType;
		typedef boost::ptr_vector<bvpfe::Triangle>    TriContainerType;
		typedef boost::ptr_map<size_t,bvpfe::Face>    FaceContainerType;
		typedef std::pair<bool,bvpfe::Tetrahedron*>   isTetrahedronType;
		
		typedef std::pair< unsigned int, std::pair<unsigned int, double> > inputBCsType;
		
        //	typedef model::VirtualBoundarySlipContainer<dim> VirtualBoundarySlipContainerType;
		
		
	public:
		NodeContainerType nodeContainer;                  // container for the nodes' pointers
		TetContainerType tetContainer;		         // container for the tetrahedrons' pointers
		std::vector<bvpfe::Triangle*> triContainer;
		
		std::vector <inputBCsType> inputBCsContainer , usrBCsContainer;      // stores the input displacement BCs from the BC_0.txt file
		
		
		std::set<unsigned int> cutTrisSet;     // triangles IDs that has dislocation segments cutting them
		
        //	VirtualBoundarySlipContainerType vbsc;       // container that stores the virtual dislocation segments outside the domain
		
#include <model/BVP/outputVTK.h>
		
		
	private:
		FaceContainerType faceContainer;                 // container for the faces pointers
		
		std::vector<double> Finf;                  // holds the r.h.s. traction (infinite + external) force vector for all the nodes
				
		std::map<bvpfe::Face*,std::vector<double> > facesWtraction;    // container for the faces that have traction, and their traction vector
		
		unsigned int sysDim;                             // dimension of the global linear system
		
		public :
		Domain () {}
		
        /**********************************************************************/
        void readMesh()
        {/*! Reads mesh data from input files mesh.*
          */
			readNodes();
			Finf.resize(nodeContainer.size()*3, 0.0);
			readVolElements();
			readSurfElements();
			outputMesh();
		}
		
        /**********************************************************************/
		void readNodes()
		{/*! Reads mesh vertex data from input files mesh.node
          */
			// int nNodes , d1, d2 , d3, in; // Giacomo 09-30-2011
			unsigned int nNodes , d1, d2 , d3, in;
			VectorDim tempP;
			
			FILE *fp =fopen("mesh/mesh.node", "r");
			
			assert(fscanf (fp, "%u%u%u%u", &nNodes , &d1, &d2, &d3)==4);
			
			if ((d1!=3)&&(d2!=0)&&(d3!=0)) assert(0 &&"Error in .node file format");
			
			float X[3];
			
			for (unsigned int i = 0; i< nNodes; i++){
				assert(fscanf (fp, "%u%f%f%f", &in , &X[0], &X[1], &X[2])==4);
				
				for (unsigned int ii = 0; ii<3; ii++) {tempP(ii) = X[ii];}
				
				std::auto_ptr<bvpfe::Node<dim> > pNode (new bvpfe::Node<dim>(tempP) );
				if(pNode->sID != in) assert(0&&"Error in .node file format. use -z option when creating the mesh to start numbering from 0");
				nodeContainer.push_back(pNode);
			}
			
			fclose(fp);
			
		}
		
        /**********************************************************************/
		void readVolElements()
		{/*! Read & distribute the volume elements data from mesh/mesh.ele and
          *  mesh/mesh.neigh
          */
			
			unsigned int nTets , d1, nn , ti , nt;
			
			FILE *fTet =fopen("mesh/mesh.ele", "r");
			FILE *fneigh =fopen("mesh/mesh.neigh", "r");
			
			assert(fscanf (fTet, "%u%u%u", &nTets , &nn, &d1)==3);
			if ((nn!=4)&&(d1!=0)) assert(0&&"Error in .ele file format");
			
			assert(fscanf (fneigh, "%u%u",  &nt, &nn)==2);
			if ((nn!=4)&&(nt!=nTets)) assert(0&&"Error in .neigh file format");
			
			int neighbor[4];
			int tetNodes[4];
			
			for (unsigned int i=0; i<nTets; i++) {				
				std::auto_ptr<bvpfe::Tetrahedron> pTet (new bvpfe::Tetrahedron () );
                
				assert(fscanf (fneigh, "%u%d%d%d%d",  &ti, &neighbor[0], &neighbor[1], &neighbor[2], &neighbor[3])==5);
				assert(fscanf (fTet, "%u%d%d%d%d",  &ti, &tetNodes[0], &tetNodes[1], &tetNodes[2], &tetNodes[3])==5);
				
				if(pTet->sID != ti) assert(0&&"Error in .ele or .neigh file format. use -z option when creating the mesh to start numbering from 0");
				
				for(unsigned int j = 0; j<nn; j++) {
					pTet->insertNode(&nodeContainer[tetNodes[j]]); // save the node ptr in the element's eleNodes array
				}
				
				//----------- set neighbor elements -----------------
				pTet->addNeighbor(neighbor);
				
				//--------- set the nodes neighborlist ---------
				pTet->setNodesNeighbors();
				
				//assert(tetContainer.insert(std::make_pair(pTet->sID,pTet)).second);
				tetContainer.push_back(pTet);
			}
			
			fclose(fTet);
			fclose(fneigh);
		}
		
		//==================================================================================
		//------------ function to read & distribute the triangular surface elements data----------
		//==================================================================================
		
		void readSurfElements()
		{
			int triNodes[3];
			unsigned int nTris , di , iFc, iTet , ti;
			
			FILE *fTri =fopen("mesh/mesh.face", "r");
			
			assert(fscanf (fTri, "%u%u", &nTris , &di)==2);
			
			if (di!=0) assert(0&&"Error in .face file format");
			
			bvpfe::Triangle* pTri;
			
			for (unsigned int iTri = 0; iTri < nTris; iTri++ ){
				assert(fscanf (fTri, "%u%d%d%d%u%u", &ti , &triNodes[0], &triNodes[1], &triNodes[2], &iTet  , &iFc )==6);
				
				//----------- if the face does NOT exist, creat it ----------
				if(faceContainer.find(iFc) == faceContainer.end()) {
					std::auto_ptr<bvpfe::Face> pFace (new bvpfe::Face);
					assert(faceContainer.insert(iFc,pFace).second);
				}
				
				//---------- add triangles ---------
				pTri = new bvpfe::Triangle;
				
				if(pTri->sID != ti) {
					assert(0&&"Error in .face file format. use -z option when creating the mesh to start numbering from 0");
				}
				
				
				for (unsigned int j = 0; j<3; j++ ) {      // insert triangle nodes
					nodeContainer[triNodes[j]].isBoundaryNode=true;
					pTri->insertNode(&nodeContainer[triNodes[j]]);
					nodeContainer[triNodes[j]].triIDs.push_back(pTri->sID);
				}
				
				triContainer.push_back(pTri);
				
				pTri->outNormal = pTri->triNormal();                       // set the outward normal for the triangle
				
				faceContainer.at(iFc).insertTri(pTri);                     // insert triangle in the face
				
				// ------------ set neighbor tetrahedron for the face triangle ----------------
				
				pTri->neighTetIndx = iTet;         // index of the triangle's neighbor tetrahedron
				tetContainer[iTet].insertSurfTri(pTri);   // set tetrahedron's neighbor surface triangle
			}
			
			setNeighborTriangles();                 // set for each triangle an array of pointers to the 3 neighbor triangles
			
			fclose(fTri);
		}
		
		//===================================================================================
		// function to set the 3 neighbor triangles for each surface triangle
		//===================================================================================
		
		void setNeighborTriangles(){
			
			//Triangle* pTri1, pTri2;
			
			Eigen::Matrix<size_t,2,3>   edgi , edgj;
			
			//unsigned int indxi, indxj;
			
			for (unsigned int i = 0; i<(triContainer.size()-1); i++){
				
				//----------- take here specific order, and take the opposite for the pTri2 to avoid the need to sort
				edgi(0,0)=triContainer[i]->eleNodes[1]->sID;     edgi(1,0)=triContainer[i]->eleNodes[2]->sID;
				edgi(0,1)=triContainer[i]->eleNodes[2]->sID;     edgi(1,1)=triContainer[i]->eleNodes[0]->sID;
				edgi(0,2)=triContainer[i]->eleNodes[0]->sID;     edgi(1,2)=triContainer[i]->eleNodes[1]->sID;
				
				for (unsigned int j = i+1; j<triContainer.size(); j++){
					
					edgj(0,0)=triContainer[j]->eleNodes[2]->sID;     edgj(1,0)=triContainer[j]->eleNodes[1]->sID;
					edgj(0,1)=triContainer[j]->eleNodes[1]->sID;     edgj(1,1)=triContainer[j]->eleNodes[0]->sID;
					edgj(0,2)=triContainer[j]->eleNodes[0]->sID;     edgj(1,2)=triContainer[j]->eleNodes[2]->sID;
					
					for (unsigned int ii = 0; ii<3; ii++){
						for (unsigned int jj = 0; jj<3; jj++){
							
							if(edgi.col(ii) == edgj.col(jj)){
								setTriCouples(i,j,ii,jj); break;
							}
							
						}
					}
					
				}
			}
			
			
			/*for (unsigned int i = 0; i<triContainer.size(); i++){
			 std::cout<< "Triangle: "<< triContainer[i]->sID<< " : " << triContainer[i]->neighbor[0]->sID << " "<< triContainer[i]->neighbor[1]->sID << " "
			 << triContainer[i]->neighbor[2]->sID << std::endl;
			 }*/
			
			
		}
		
		//==================================================================================
		// function to set 2 neighbor triangles as neighbors
		//===================================================================================
		
		void setTriCouples(unsigned int ti,unsigned int tj,unsigned int ni,unsigned int nj){
			
			//-----temporary, and needs modification. Index modification is mainly becasue of different sequence convention for edgj
			if(nj==2) {nj=1;}
			else if (nj==1) {nj=2;}
			
			//std::cout << ti << " " << tj << " " << ni << " " << nj << std::endl;
			
			triContainer[ti]->neighbor[ni] = triContainer[tj];
			triContainer[tj]->neighbor[nj] = triContainer[ti];
			
		}
		
		
		//==================================================================================
		// function to read and store the input traction and displacement boundary conditions for the domain from BC_0.txt
		//===================================================================================
		
		void readInputBCs() {
			
			Eigen::Matrix <float,dim,1>  u , tr;
			Eigen::Matrix <unsigned int,dim,1> isBC;
			unsigned int ni;
			
			inputBCsType u_BC;
			
			FILE *fp =fopen("mesh/BCs_0.txt", "r");
			
			while (!feof(fp)) {
				if(fscanf(fp, "%u%u%u%u%f%f%f%f%f%f", &ni,&isBC(0), &isBC(1), &isBC(2),&u(0),&u(1),&u(2),&tr(0),&tr(1),&tr(2))==10){
					if(nodeContainer[ni].isBoundaryNode) {
						for (unsigned int j = 0; j<3 ; j++) {
							//if(isBC(j))  nodeContainer[ni].setBC(j ,double(u(j)) );
							if(isBC(j))  {
								//std::cout << j << "   " << nodeContainer[ni].P(j) << "  " << double(u(j)) << std::endl;
								u_BC.first = ni;
								u_BC.second.first = j;
								u_BC.second.second = double(u(j)) ;
								inputBCsContainer.push_back(u_BC);
							}
							nodeContainer[ni].traction(j) = double(tr(j));
						}
					}
					else  {assert(0&&"inputting boundary condition for non-boundary node");}
				}
			}
			
			fclose(fp);
		}
		
		//==================================================================================
		// function to set the displacement boundary conditions for the domain
		// 1- Update the infinite medium displacement at domain nodes.
		// 2- set the the displacement BCs that were defined in the BC_0.txt file.
		// 3- Update user defined displacement BCs, if required.
		//===================================================================================
		template<typename T>
		void setBoundaryConditions(bool update_usr_defined_BCs ,const T* const pT) {
			
			//------------------- remove old boundary conditions first "just for safety" ----------
			removeBoundaryConditions();
			
			std::cout<<"======================================================================================="<<std::endl;
			//-------------------- Calculate infinite medium displacement field --------------------
			
			//FILE *finfd =fopen("infDis.txt", "w");
			
            // 			double err = 0.0e00;
            // 			unsigned int nn = 0;
            // 			double errMax = 0.0e00;
            // 			unsigned int inn;
            // 			VectorDim u1, u2;
			
			double t0=clock();
			std::cout<<"Calculate infinite medium displacement field......  ";
            
			for (unsigned int dN=0; dN<nodeContainer.size();++dN){
				if(nodeContainer[dN].isBoundaryNode) {
                    if(nodeContainer[dN].triIDs.size() == 0) assert(0&&"Error: Boundary node without triangle element index array");
                    
                    // 				  VectorDim d1 = pT->displacement        (nodeContainer[dN].P,triContainer[nodeContainer[dN].triIDs[0]]->outNormal);
                    // 				  VectorDim d2 = pT->displacementStraight(nodeContainer[dN].P,triContainer[nodeContainer[dN].triIDs[0]]->outNormal);
                    // 				  //VectorDim d2 = pT->displacement_straight(nodeContainer[dN].P);
                    // 				  VectorDim error =  d1- d2;
                    
                    //std::cout<< d1.transpose() << "  " << d2.transpose() << "  " << error.transpose() << std::endl;
                    /*
                     err+=	error.norm()/d1.norm();
                     if (error.norm()/d1.norm()>errMax) {errMax = error.norm()/d1.norm();  u1 = d1;   u2=d2; inn = dN;}
                     nn++;*/
                    nodeContainer[dN].uInf=pT->displacement(nodeContainer[dN].P,triContainer[nodeContainer[dN].triIDs[0]]->outNormal) + nodeContainer[dN].uVir;
                    // nodeContainer[dN].uInf=pT->displacement_straight(nodeContainer[dN].P);
                    //nodeContainer[dN].uInf= error;
					//  fprintf (finfd, "%u %f %f %f \n",dN , nodeContainer[dN].uInf(0) , nodeContainer[dN].uInf(1), nodeContainer[dN].uInf(2) );
				}
				
				
				
			}
			
			/*
             for (unsigned int dN=0; dN<nodeContainer.size();++dN){
             
             Eigen::Matrix<double,3,3> d1 = pT->stress(nodeContainer[dN].P);
             Eigen::Matrix<double,3,3> d2 = pT->stress_straight_CoorDpndt(nodeContainer[dN].P);
             //Eigen::Matrix<double,3,3> d2 = pT->stress_straight(nodeContainer[dN].P);
             Eigen::Matrix<double,3,3> error =  d1- d2;
             
             //std::cout<< d1.col(0).transpose() << "  " << d1.col(1).transpose() << "  " << d1.col(2).transpose() << std::endl;
             //std::cout<< d2.col(0).transpose() << "  " << d2.col(1).transpose() << "  " << d2.col(2).transpose() << std::endl;
             //std::cout<< std::endl;
             
             err+=	error.norm()/d1.norm();
             
             nn++;
             }
             */
			std::cout<<"DONE in"<<"["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]"<<std::endl;
			
			//std::cout<< "displacement calculation error =  " <<std::setprecision(15) << err/nn << "     " << errMax << "  " << inn << "  " << u1.transpose() << "  " << u2.transpose() << std::endl;
			
			//fclose(finfd);
			//------------------- update input boundary conditions -------------------------------
			
			for (unsigned int i=0; i<inputBCsContainer.size(); i++){
				nodeContainer[inputBCsContainer[i].first].setBC( inputBCsContainer[i].second.first , inputBCsContainer[i].second.second );
			}
			
			//----------------------- update the user defined boundary conditions, if needed -------------------------
            if (update_usr_defined_BCs){
#ifdef UpdateBoundaryConditionsFile
                usrBCsContainer = pT->update_usr_BCs();        // store user defined boundary conditions
#endif
            }
            
#ifdef UpdateBoundaryConditionsFile
            // apply user defined boundary conditions
            for (unsigned int i=0; i<usrBCsContainer.size(); i++){
                nodeContainer[usrBCsContainer[i].first].setBC( usrBCsContainer[i].second.first , usrBCsContainer[i].second.second );
            }
#endif
            //			if (update_usr_defined_BCs){
            //#ifdef UpdateBoundaryConditionsFile
            //			    usrBCsContainer = pT->update_usr_BCs();        	// store user defined boundary conditions
            //
            //			// apply user defined boundary conditions
            //			for (unsigned int i=0; i<usrBCsContainer.size(); i++){
            //				nodeContainer[usrBCsContainer[i].first].setBC( usrBCsContainer[i].second.first , usrBCsContainer[i].second.second );
            //            }
            //
            //#endif
            //            }
			
		}
		
		//========================================================================================
		// Function to set the boundary conditions for the BVP and then solve it
		//=======================================================================================
		template<typename T>
		void update_BVP_Solution(bool update_usr_BCs ,const T* const pT) {
			
            //pT->findSegmentsCuttingTriangles();
            
			//------- Calculate infinte medium displacement, and set displacement boundary conditions ----------
			setBoundaryConditions(update_usr_BCs , pT);
			
			//--------------------- Solve the BVP -----------------------------------------------------------
			solveBVP<T>(true,pT);
		}
		
		//===================================================================================
		// function to assemble the global linear system, solve it, and distribute nodal displacement
		//===================================================================================
		
		template<typename T>
		void solveBVP(bool dislocations_Stress,const T* const pT)
		{
			double t0;
			
			//std::cout<<"======================================================================================="<<std::endl;
			//-------------------- Calculate infinite medium displacement field --------------------
			/*
			 t0=clock();
			 std::cout<<"Calculate infinite medium displacement field......  ";
			 for (unsigned int dN=0; dN<nodeContainer.size();++dN){
			 if(nodeContainer[dN].isBoundaryNode) nodeContainer[dN].uInf=pT->displacement(nodeContainer[dN].P);
			 }
			 std::cout<<"DONE in"<<"["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]"<<std::endl;
			 
			 */
			//---------------------- initiate sparse system -------------------------
			std::vector<int> rr, cc;
			unsigned int  nz;
			
			
			double tS1=clock();
			//------ set the equation number for each dof based on the boundary conditions
			setEquationNumbers();
			
			//-------- initiate the sparse system ------------------
			
			initiateSparseSystem(nz, rr , cc);
			
			//--------------------- Build the sparse matrix ------------------------
			
			Eigen::VectorXi nColNZs(sysDim);           // number of nonzero elements on each column of the sparse matrix
			
			for (unsigned int i=0; i<sysDim; i++) {nColNZs(i)=0;}            // initialize
		    
		    for ( unsigned int i=0; i<nz; i++) {
				++nColNZs (cc[i]);
		    }
		    
		    Eigen::VectorXd U(Eigen::VectorXd::Zero(sysDim));           // Displacement vector
		    Eigen::VectorXd F(Eigen::VectorXd::Zero(sysDim));           // r.h.s. force vector
		    for(unsigned int i=0; i<sysDim; i++) {F(i) = 0.0e+00; U(i) =  0.0e+00;}
		    
		    Eigen::SparseMatrix<double> K(sysDim,sysDim);         // The sparse global stiffness matrix
		    K.reserve(nColNZs);
		    
		    
		    //-------- Assemble the the linear system ---------
		    t0=clock();
		    std::cout<<"Assembling stiffness matrix.......................  "<<std::flush;
		    
		    assemble (K,F);    // assmble stiffness matrix and force vector from displacment BCs
		    
		    K.makeCompressed();
		    std::cout<<"DONE in"<<"["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]"<<std::endl;
		    
		    
            
		    //-------- Assemble the r.h.s force vector from infinite and external traction ---------
            
		    if (dislocations_Stress){
				t0=clock();
				std::cout<<"Assembling infinite traction r.h.s. vector .......  "<<std::flush;
				assembleInfiniteTraction <T> (F,pT);   // assemble force vector from infinite medium surface traction
				std::cout<<"DONE in"<<"["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]"<<std::endl;
		    }
			
            
		    //---------------- Solving linear system ---------------------------
		    t0=clock();
		    std::cout<<"Solving linear system ............................  "<<std::flush;
            if (F.squaredNorm()>0.0){
                Eigen::ConjugateGradient<Eigen::SparseMatrix<double> > solver;
                solver.setTolerance(1.0e-8);
                U = solver.compute(K).solveWithGuess(F,Eigen::VectorXd::Zero(sysDim));
                std::cout<<" "<< solver.info() << " " ;
                if(solver.info()!= Eigen::Success) {
                    std::cout << "FAILED: " << solver.info() << std::endl;
                    assert(0);
                }
                //model::SequentialOutputFile<'F',1>::set_increment(outputFrequency); // Vertices_file;
                //model::SequentialOutputFile<'F',1>::set_count(runID); // Vertices_file;
                //                model::SequentialOutputFile<'F',true> f_file;
                //                for (unsigned int i = 0 ; i < nodeContainer.size(); i++ )
                //                {
                //                    if (nodeContainer[i].P(2)==4000)
                //                    {
                //                        f_file<< nodeContainer[i].P.transpose()<<" ";
                //                        for (int j=0;j<3;++j)
                //                        {
                //                            f_file<<  F(nodeContainer[i].eqnNumber[j])<<" ";
                //                        }
                //                        f_file<<std::endl;
                //                    }
                //                }
                
            }
            
            
		    
		    std::cout<<"DONE in"<<"["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]"<<std::endl;
		    
		    //std::cout<<"Sparse Eigen time= "<<(clock()-tS0)/CLOCKS_PER_SEC<<" [sec]"<<std::endl;
		    
		    
		    
		    //-------------- distribute displacement to nodes -----------
		    for (unsigned int i = 0 ; i < nodeContainer.size(); i++ )
		    {
				for (int j=0; j<dim; j++)
				{
					if (!nodeContainer[i].isBC(j)) {nodeContainer[i].u(j) = U(nodeContainer[i].eqnNumber(j));}
				}
				nodeContainer[i].displaceNode();
				//std::cout<< "node: "<< nodeContainer[i].sID<< ": " <<nodeContainer[i].u.value[0]<< " " <<nodeContainer[i].u.value[1]<< " " <<nodeContainer[i].u.value[2] << std::endl;
		    }
		    
		    std::cout<<"Total BVP solution time=    "<<(clock()-tS1)/CLOCKS_PER_SEC<<" [sec]"<<std::endl;
		    std::cout<<"======================================================================================="<<std::endl;
#ifdef Contact_Loading
		    while(pT->overConstraintNodes()) {
				t0=clock();
				std::cout<<"Resolving BVP .........................................  "<<std::flush;
				reSolveBVP(true);
				std::cout<<"DONE in"<<"["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]"<<std::endl;
		    }
#endif
		    std::cout<<"======================================================================================="<<std::endl;
		    //output();
		    //writeVTK_file (pT);
		}
		
		
		//===================================================================================
		// function to assemble the global linear system, solve it, and distribute nodal displacement
		//===================================================================================
		
		void reSolveBVP(bool dislocations_Stress)
		{
		  	std::vector<int> rr, cc;
			unsigned int  nz;
			
			
			//double tS1=clock();
			//------ set the equation number for each dof based on the boundary conditions
			setEquationNumbers();
			
			//-------- initiate the sparse system ------------------
			
			initiateSparseSystem(nz, rr , cc);
			
			//--------------------- Build the sparse matrix ------------------------
			
			Eigen::VectorXi nColNZs(sysDim);           // number of nonzero elements on each column of the sparse matrix
			
			for (unsigned int i=0; i<sysDim; i++) {nColNZs(i)=0;}            // initialize
			
			for (unsigned int i=0; i<nz; i++) {
				++nColNZs (cc[i]);
			}
			
			Eigen::VectorXd U(sysDim);           // Displacement vector
			Eigen::VectorXd F(sysDim);           // r.h.s. force vector
			for(unsigned int i=0; i<sysDim; i++) {F(i) = 0.0e+00; U(i) =  0.0e+00;}
			
			Eigen::SparseMatrix<double> K(sysDim,sysDim);         // The sparse global stiffness matrix
			K.reserve(nColNZs);
			
			
			//-------- Assemble the the linear system ---------
			
			assemble (K,F);    // assmble stiffness matrix and force vector from displacment BCs
			
			K.makeCompressed();
			
			//assembleTractionForce (F);   // assemble force vector from surface traction
			
			if (dislocations_Stress) copyTractionVector(F);
			
			//double tS0=clock();
			
			Eigen::ConjugateGradient<Eigen::SparseMatrix<double> > solver;
			U = solver.compute(K).solve(F);
			
			//std::cout<<"Sparse Eigen time= "<<(clock()-tS0)/CLOCKS_PER_SEC<<" [sec]"<<std::endl;
			//std::cout<<"Total BVP time=    "<<(clock()-tS1)/CLOCKS_PER_SEC<<" [sec]"<<std::endl;
			
			//-------------- distribute displacement to nodes -----------
			for (unsigned int i = 0 ; i < nodeContainer.size(); i++ )
			{
				for (int j=0; j<dim; j++)
				{
					if (!nodeContainer[i].isBC(j)) {nodeContainer[i].u(j) = U(nodeContainer[i].eqnNumber(j));}
				}
				nodeContainer[i].displaceNode();
				//std::cout<< "node: "<< nodeContainer[i].sID<< ": " <<nodeContainer[i].u.value[0]<< " " <<nodeContainer[i].u.value[1]<< " " <<nodeContainer[i].u.value[2] << std::endl;
			}
			
			
			
		}
		
		//=================================================================================================
		// Copy the r.h.s. traction (infinite and external) force vector that was calculated inside solveBVP() function
		//=================================================================================================
		
		void copyTractionVector(Eigen::VectorXd& Fm) {
			
			int qi;
			
			for (unsigned int ni = 0; ni < nodeContainer.size(); ni ++ ) {
				
				for (unsigned int dofi = 0; dofi<3; dofi++){
					
					qi = nodeContainer[ni].eqnNumber(dofi);
					
					if ( qi != -1) Fm(qi) = Fm(qi) + Finf[(ni*3)+dofi];
					
				}
				
			}
		}
		
		//===================================================================================
		// function to set the equation numbers for all dofs
		//===================================================================================
		
		void setEquationNumbers()
		{
			//unsigned int dim = 3;       // 3 dofs
			int num;
			bvpfe::Node<dim>* pNode;
			
			num = -1;
			
			for (unsigned int i = 0 ; i < nodeContainer.size(); i++ )
			{
				pNode = &nodeContainer[i];
				for (int j=0; j<dim; j++)
				{
					if (!pNode->isBC(j))
					{
						num = num + 1;
						pNode->eqnNumber(j) = num;
					}
				}
			}
			
			sysDim = num + 1;          // set linear system dimension
			
			//std::cout<< "number of dof = " << sysDim << std::endl;
			
			/*for (int i = 0 ; i < nodeContainer.size(); i++ )
			 {
			 pNode = nodeContainer[i];
			 std::cout<< "Node=   "<< pNode->sID << ":  " << pNode->u.isBC[0] << "  " << pNode->u.isBC[1] << "  "<< pNode->u.isBC[2] << "  "<< pNode->u.eqnNumber[0] << "  " << pNode->u.eqnNumber[1] << "  "<< pNode->u.eqnNumber[2] << "  "<< std::endl;
			 }*/
			
			
		}
		
		//===================================================================================
		// function to assemble the global stiffness matrix
		//===================================================================================
		
		void assemble (Eigen::SparseMatrix<double> & Km, Eigen::VectorXd& Fm)
		{
			//bvpfe::Tetrahedron* pTet;
			//std::vector<Vector> tetMat;         // 12x12 stiffness matrix for each tetrahedron
			Eigen::Matrix<double,12,12> tetMat;         // 12x12 stiffness matrix for each tetrahedron
			
			//std::vector<Vector> Np;             // shape functions derivatives mapped to actual element
			Eigen::Matrix<double,dim,bvpfe::Tetrahedron::Nnodes> Np;   // shape functions derivatives mapped to actual element
			
			int qi, qj;
			int ie , je;
			
			double wv;
			
			for (unsigned int i=0; i<tetContainer.size() ; i++) {
				
				//pTet = &tetContainer[i];
				
				tetMat = tetContainer[i].getElementStiffness(Np, wv);
				
				//---------- assemble tetMat to the global matrix ---------------------------
				for ( unsigned int ai = 0 ; ai < 4 ; ai ++)            // no. of nodes per Tet
				{
					ie= ai*3;
					for ( unsigned int in = 0 ; in < 3 ; in ++)    //no. of DOFs per node
					{
						qi = tetContainer[i].eleNodes[ai]->eqnNumber[in];
						if (qi == -1) continue;
						
						for ( unsigned int bi = 0 ; bi < 4 ; bi ++)
						{
							je = bi*3 ;
							for ( unsigned int jn = 0 ; jn < 3 ; jn ++)
							{
								qj = tetContainer[i].eleNodes[bi]->eqnNumber(jn);
								if (qj == -1)
								{
									Fm(qi)=Fm(qi)-(tetMat(ie+in,je+jn)*tetContainer[i].eleNodes[bi]->bcValue[jn]);
									//std::cout<< "out111=== " << tetMat(ie+in,je+jn) << "  " << tetContainer[i].eleNodes[bi]->u.bcValue[jn] << std::endl;
									continue;
								}
								
								//Km.set(qi,qj) = Km(qi,qj) + (tetMat(ie+in,je+jn));
								Km.coeffRef(qi,qj) += tetMat(ie+in,je+jn);
							}
						}
					}
				}
				
				
				//------------- assemble the dislocation force vector, if required ----------
				
				/*if(disStress)
				 {
				 
				 Eigen::Matrix<double,dim,Tetrahedron::Nnodes> temp=tetContainer[i].get_dislocationForceVector<4,T>(pT,Np);
				 
				 
				 //----------- assemble in the global force vector ------------------
				 for ( unsigned int in = 0 ; in < Tetrahedron::Nnodes ; in ++)            // no. of nodes per Tet
				 {
				 for ( unsigned int id = 0 ; id < dim ; id ++)    //no. of DOFs per node
				 {
				 qi = tetContainer[i].eleNodes[in]->u.eqnNumber[id];
				 if (qi == -1) continue;
				 
				 Fm(qi)=Fm(qi)-temp(id,in);
				 
				 }
				 }
				 }*/
				
			}
		}
		
		
		//===================================================================================
		// function to initiate the sparse linear system
		//===================================================================================
		
		void initiateSparseSystem(unsigned int& nz, std::vector<int>& rr, std::vector<int>& cc )
		{
			int qi, qj;
			
			bvpfe::Node<dim>* pNodei;
			bvpfe::Node<dim>* pNodej;
			
			// ---- sort the neighbor vector for each node. This makes it easy to initiate sparse system ---
			
			for (unsigned int i=0; i<nodeContainer.size() ; i++) {
				nodeContainer[i].sortNeighbors();
				//				pNodei = &nodeContainer[i];
				//				pNodei->sortNeighbors();
				
				/*std::cout<< pNodei->sID<< " : " ;
				 for (int j=0; j<pNodei->neighbor.size() ; j++)
				 {
				 std::cout<< pNodei->neighbor[j]->sID<< "  " ;
				 }
				 std::cout<< std::endl;*/
			}
			
			//---- loop over all nodes to get the positions of nonzero elements------
			
			for (unsigned int i=0; i<nodeContainer.size() ; i++) {
				
				pNodei = &nodeContainer[i];
				
				for (int ii = 0; ii<dim; ii++)
				{
					qi = pNodei->eqnNumber(ii);
					if (qi == -1) continue;
					
					for (unsigned int j = 0; j<pNodei->neighbor.size(); j++)
					{
						pNodej = pNodei->neighbor[j];
						
						for (int jj = 0; jj<dim; jj++)
						{
							qj = pNodej->eqnNumber(jj);
							if (qj == -1) continue;
							
							rr.push_back(qi);     cc.push_back(qj);
						}
					}
				}
				
			}
			
			nz = rr.size();
			
			//std::cout << "nz = " << nz << std::endl;
			
			/*for (int i = 0; i<rr.size() ; i++)
			 {
			 std::cout << rr[i] << "  "<< cc[i] << std::endl;
			 }*/
			
			
		}
		
		//===================================================================================
		// function to assemble the force vector from infinite medium surface traction
		//===================================================================================
		template<typename T >
		void assembleInfiniteTraction (Eigen::VectorXd& Fm , const T* const pT) {
			
			for (unsigned int i=0 ; i< Finf.size() ; i++ )  Finf[i]= 0.0;
			
            //		Eigen::Matrix<double,dim,3> TriVec;       // 3 is the number of nodes per triangle
			
			unsigned int nID;
			bvpfe::Triangle* pTri;
			std::vector<bvpfe::Triangle*>::iterator itt;
			
			
			// Loop and compute infinite force and store in each triangle
			
			
			
            //				for (itt= triContainer.begin() ; itt != triContainer.end(); itt++ ) {
            
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (unsigned int kkk=0; kkk<triContainer.size(); kkk++ ) {
                //	std::vector<Triangle*>::iterator itt(triContainer.begin());
                //	std::advance(itt,kkk);
                //	pTri = *itt;
                //	pTri->TriVec = pTri->getTriInfiniteForce_gp<T,LinkType>(pT, gpObsever);
                triContainer[kkk]->TriVec = triContainer[kkk]->getTriInfiniteForce_gp<T>(pT);
			}
            
            
            
            //            model::SequentialOutputFile<'H',true> h_file;
            //            for (unsigned int kkk=0; kkk<triContainer.size(); kkk++ ) {
            //                for (int nnn=0;nnn<3;++nnn){
            //                    h_file<< triContainer[kkk]->eleNodes[nnn]->P.transpose()<<" "<< triContainer[kkk]->TriVec.col(nnn).transpose()<<" ";
            //                }
            //                h_file<<std::endl;
            //			}
			
			
			
			for (itt= triContainer.begin() ; itt != triContainer.end(); itt++ ) {
				
				pTri = *itt;
				
				
				//------------- assemble the infinite traction vector --------------
				int  qi;
				for (int ai = 0 ; ai < 3; ai++ )             // loop over the 3 surface element's nodes
				{
					nID = pTri->eleNodes[ai]->sID;        // node global index
					
					for  (int in = 0 ; in < 3; in++ )   // loop over the 3 dofs
					{
						Finf[(nID*3)+in] = Finf[(nID*3)+in] + pTri->TriVec(in,ai);  // -ve sigen for the infinite field is considered in Triangle::dislocationStressKernel
						
						qi = pTri->eleNodes[ai]->eqnNumber(in);
						if (qi != -1) Fm(qi) = Fm(qi) + pTri->TriVec(in,ai);  // -ve sigen for the infinite field is considered in Triangle::dislocationStressKernel
					}
				}
				
			}
			
			
		}
		
        
        
        /**********************************************************************/
		Eigen::Matrix<double,dim,dim> stressAt(const VectorDim& P)
		{/*! @param[in] P the field point
          * The stress tensor at the field point.
          */
			Eigen::Matrix<double,dim,dim> stress(Eigen::Matrix<double,dim,dim>::Zero());
			const isTetrahedronType isT(findIncludingTet(P));
			if (isT.first)
            {
                stress = isT.second->getStress();
            }
			return stress;
		}
		
		//===================================================================================
		// function to search for the including tetrahedron given the new position for the point,
		// old point's tetrahedron index, and search direction vector
		//===================================================================================
		
		void SearchMovingNode (model::SearchData<dim> & data){
			unsigned int nn ;
			// --------------NOT boundary node (coming from inside or outside) -----------------------
			if(data.nodeMeshLocation == 0 || data.nodeMeshLocation == 1) {
				
				data.newMeshID = data.currentMeshID;	    // initialize
				int dI, sI;
				sI = 0;  dI = 0;
				Eigen::Matrix<double,4,1>::Index ii;
				// , jj;
				VectorDim interP;                              // should save the velocity-surface intersection point
				bvpfe::Triangle* pTri;
				bvpfe::Tetrahedron* pTet;
				
				Eigen::Matrix<double,4,1>  Bary;
				
				while (!data.found){
					bool insideTet =  tetContainer[data.newMeshID].isInsideTet(data, sI, Bary);
					
					if(insideTet) break;
					
					//--------if node is found outside, find intersection with boundary, and set it as boundary node ----
					if(data.nodeMeshLocation == 0){
						pTet = &tetContainer[data.newMeshID];
						nn = 0;
						while(!data.found){
							
							//---- find the vector-surface intersection point ---------
							data.newMeshID = pTet->sID;
							
							interP = getSurfaceIntersection(pTet->TetSurfTris.find(sI)->second->outNormal,
															pTet->TetSurfTris.find(sI)->second->eleNodes[0]->P, -data.normalizedDir, data.P);
							
							//------ calculate the barycentric coordinates for this intersection point -------
							Eigen::Matrix<double,4,1> bary = pTet->getBarycentric(interP);
							double baryMin = bary.minCoeff(&ii);
							//double baryMin = pTet->getBarycentric(interP).minCoeff(&ii);
							
							//------ check if the point is inside, or just outside (but too close) the tetrahedron surface --> DONE ---------
							//------ if (ii==sI) so the point is on the tetrahedron surface, with small error from calculations ---
							//if(ii==sI) std::cout<< "minBary = " << baryMin << std::endl;
							if( baryMin>=-1.0e-12 && baryMin<=1.0e-12   /*||(ii==sI)*/){
                                if(pTet->neighbor[ii]==-1){
                                    data.found=true;
                                    data.projectedP = interP;
                                    data.outwardFaceNormal= pTet->TetSurfTris.find(ii)->second->outNormal;
                                    data.nodeMeshLocation = 2;      // boundary node
                                    data.triIndex = pTet->TetSurfTris.find(ii)->second->sID;
                                    break;
                                }
                                else{
                                    for ( int ib=0; ib<4 ; ib++){
                                        if (ib==ii) continue;
                                        if(bary(ib)>=-1.0e-12 && bary(ib)<=1.0e-12 && pTet->neighbor[ib]==-1) {
                                            data.found=true;
                                            data.projectedP = interP;
                                            data.outwardFaceNormal= pTet->TetSurfTris.find(ib)->second->outNormal;
                                            data.nodeMeshLocation = 2;      // boundary node
                                            data.triIndex = pTet->TetSurfTris.find(ib)->second->sID;
                                            break;
                                        }
                                    }
                                }
							}
							if (data.found) break;
							//-------- if the intersection point is not in this tetrahedron, move to the next -----
							
							// --- find the index (from 0 to 2) for the point with lowest bary for the surface triangle ----
							
							for(unsigned int k=0; k<3; k++){
								if(pTet->TetSurfTris.find(sI)->second->eleNodes[k]->sID == pTet->eleNodes[ii]->sID) {dI = k; break;}
							}
							
							pTri = pTet->TetSurfTris.find(sI)->second->neighbor[dI];    // next neighbor triangle
							pTet = &tetContainer[pTri->neighTetIndx];                   // new tetrahedron to move to
							
							for (std::map<size_t,bvpfe::Triangle*>::iterator it = pTet->TetSurfTris.begin(); it!= pTet->TetSurfTris.end(); it++){
								if (it->second->sID == pTri->sID){sI = it->first; break;}
							}
							
							nn++;
							if (nn>100003) {
                                std::cout<< "minBary = " <<std::setprecision(15) << baryMin<< " " << pTet->neighbor[ii] << std::endl;
                                assert(0&&"The code went inside an infinite loop: Failed to find node-surface intersection point");
							}
						}
					}
				}
			}
			
			// ----------------- boundary node ---------------------
			
			else if (data.nodeMeshLocation == 2) {
				
				VectorDim org = data.P - data.Dir;      // old position of the node
				unsigned int iTri;
				
				
				if(triContainer[data.triIndex]->intersectWithLine(org,data.P,data.projectedP,iTri)) {
                    //VectorDim avgNormal = data.outwardFaceNormal + triContainer[iTri]->outNormal;   // average of old & new Triangle normals to avoid vibrations when segments has node that is commuting between two neighbor triangles
                    //std::cout<<std::setprecision(15)<< data.outwardFaceNormal.transpose()<< std::endl;
                    //std::cout<<std::setprecision(15)<< triContainer[iTri]->outNormal.transpose()<< std::endl;
                    //std::cout<<std::setprecision(15)<< avgNormal.normalized().transpose()<< std::endl;
                    //std::cout<< std::endl;
                    //data.outwardFaceNormal = avgNormal.normalized();
                    
                    data.outwardFaceNormal= triContainer[iTri]->outNormal;
                    data.newMeshID = triContainer[iTri]->neighTetIndx;
                    data.triIndex = iTri;
                    //std::cout<< "new mesh ID ================== "<<data.newMeshID  << std::endl;
				}
				else{
                    data.newMeshID = data.currentMeshID;
                    data.outwardFaceNormal= triContainer[data.triIndex]->outNormal;
				}
				data.found=true;
			}
			
			
		}
		//===================================================================================
		// function to calculate the intersection point between a line and a plane, given
		// the plane normal, point on the plane, line direction, point on the line
		//==================================================================================
		
		VectorDim getSurfaceIntersection (VectorDim pD, VectorDim pP,const VectorDim lD,const VectorDim lP){
			
			double dm = lD.dot(pD);
			
			if(dm == 0.0) assert(0&&"Dislocation Node is moving parallel to the boundary: unable to get intersection point");
			
			double d = ((pP-lP).dot(pD))/dm;
			
			return (lP+d*lD);
			
		}
		
		//===================================================================================
		// function to return a pointer to the including tetrahedron for any given point,
		// or return null if the point is outside the domain
		//==================================================================================
		isTetrahedronType findIncludingTet(VectorDim P){
			bool found = false;
			int ci; // index of neighboor tet. -1 means outside.
			int ni; //
			
			isTetrahedronType temp=std::make_pair(true, (bvpfe::Tetrahedron*) NULL);
			
			//----- starting Tet, random start point ------------
			ci = int(tetContainer.size()/2);
			
			while ((!found)&&(ci>=0))
			{
				temp.second = &tetContainer[ci];
				ci = temp.second->nextNeighbor(P,found,ni);    // returns -10 if the point is inside the Tet
			}
			
			temp.first=found;
			
			//			if (!found)   // this means it is outside the domain
			//			{
			//				temp.first=false;
			//			//	temp.second = NULL;
			//				//std::cout << "============= ERROR: Searching for point outside the domain ==============" << std ::endl;
			//				//assert(0);
			//			}
			
			return temp;
		}
		
		//====================================================================================
		// Function to initally find the dislocation node position in the FE mesh
		//=====================================================================================
		void findIncludingTet (model::SearchData<dim> & data){
            
            //---- initialize ---------------
            //data.currentMeshID =
            data.newMeshID = int(tetContainer.size()/2);
            
            bool found = false;
            
            while (!found){
                found = tetContainer[data.newMeshID].isInsideTet(data) || data.nodeMeshLocation != -1;
            }
            
		}
		
		
		//====================================================================================
		// Function to initally find the dislocation node position in the FE mesh
		//=====================================================================================
		void findIncludingTet (model::SearchData<dim> & data, const unsigned int& tetID){
            
            //---- initialize ---------------
            data.newMeshID = tetID;
            
            bool found = false;
            
            while (!found){
                found = tetContainer[data.newMeshID].isInsideTet(data) || data.nodeMeshLocation != -1;
            }
            
		}
		
		//===================================================================================
		// function to detect the intersection line between a slip plane, and the mesh surface
		//==================================================================================
		//void get_planeMeshIntersection(const VectorDim& x0, const VectorDim& n, std::map< unsigned int, std::pair<VectorDim,VectorDim> >& collisionContainer, const unsigned int gpID){
		void get_planeMeshIntersection(const VectorDim& x0, const VectorDim& n,std::map< unsigned int, std::pair<VectorDim,VectorDim>,std::less<unsigned int>, Eigen::aligned_allocator<std::pair<const unsigned int,std::pair<VectorDim,VectorDim> > > >& collisionContainer){
			
			double tol = 1.0e-7;
            
			//FILE *ft =fopen("triangles_gp.txt", "a");
			
			for(unsigned int i=0; i<triContainer.size();i++){
				std::vector<VectorDim> intersectionPoints;
				
				for(unsigned int j = 0; j< 3; j++){
					VectorDim v0 = triContainer[i]->eleNodes[j]->P;
					VectorDim v1 = triContainer[i]->eleNodes[(j+1)%3]->P;
					double u = getSurfaceIntersectionPar (x0,n,v0,v1);
					
					if(u>=0 && u<=1.0){
						VectorDim P = v0 + u*(v1-v0);
						
						
						bool isDifferent=true;
						for (unsigned int k=0;k<intersectionPoints.size();k++){
							isDifferent*= (P-intersectionPoints[k]).squaredNorm()>tol;
						}
						
						if (isDifferent){
							intersectionPoints.push_back(P);
						}
					}
				}
				
				assert(intersectionPoints.size()<3);
				
				
				
				if (intersectionPoints.size()==2){
					//collisionContainer.push_back(std::make_pair(intersectionPoints[0],intersectionPoints[1]));
					collisionContainer.insert(std::make_pair(triContainer[i]->sID,std::make_pair(intersectionPoints[0],intersectionPoints[1])));
					
					
					//---------------- generate custom Quadrature points around the triangle-glide plane intersection line------
					//triContainer[i]->makeLocalQuadPoints(intersectionPoints[0],intersectionPoints[1] , gpID);
					
                    // 					fprintf (ft, "%22.15e %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e \n",
                    // 						 triContainer[i]->eleNodes[0]->P(0),triContainer[i]->eleNodes[0]->P(1), triContainer[i]->eleNodes[0]->P(2),
                    // 						 triContainer[i]->eleNodes[1]->P(0),triContainer[i]->eleNodes[1]->P(1), triContainer[i]->eleNodes[1]->P(2),
                    // 						 triContainer[i]->eleNodes[2]->P(0),triContainer[i]->eleNodes[2]->P(1), triContainer[i]->eleNodes[2]->P(2));
                    
					
                    // 				std::cout<<"Tri ID  " << triContainer[i]->sID << std::endl;
                    // 				std::cout<< std::setprecision(15) <<triContainer[i]->eleNodes[0]->P.transpose() << std::endl;
                    // 				std::cout<< std::setprecision(15) <<triContainer[i]->eleNodes[1]->P.transpose() << std::endl;
                    // 				std::cout<< std::setprecision(15) <<triContainer[i]->eleNodes[2]->P.transpose() << std::endl;
                    // 				std::cout<< std::setprecision(15) << intersectionPoints[0].transpose() << "   " << intersectionPoints[1].transpose() << std::endl;
                    // 				std::cout<<  std::endl;
					
				}
				
			}
			
			//fclose(ft);
            
			
		}
		
		
		//===================================================================================
		// function to calculate the intersection point between a line and a plane, given
		// the plane normal, point on the plane, line direction, point on the line
		//==================================================================================
		
		double getSurfaceIntersectionPar (VectorDim x0, VectorDim n,const VectorDim v0,const VectorDim v1){
			
			return ((x0-v0).dot(n)) / (v1-v0).dot(n);
		}
		
		//===================================================================================
		// function to output Tetrahedron mesh
		//==================================================================================
		//		void outputMesh(const size_t& runID, const size_t& outputFrequency){
		void outputMesh() const {
			
			//model::SequentialOutputFile<'N',1>::set_increment(outputFrequency); // Vertices_file;
			//model::SequentialOutputFile<'N',1>::set_count(runID); // Vertices_file;
			model::SequentialOutputFile<'N',true> nodesFile;
			//			model::UniqueOutputFile<'N'> nodesFile;
			for (NodeContainerType::const_iterator iter=nodeContainer.begin();iter!=nodeContainer.end();++iter){
				nodesFile<< iter->sID<<"	" << iter->P.transpose()<<"	"<<iter->isBoundaryNode<<std::endl;
			}
			
			//output_T_File();
			
			//		}
			
			//======================================================================================================
			// Function to output the T files that contains the mesh connectivity data
			//=====================================================================================================
			//		void output_T_File(){
			model::SequentialOutputFile<'T',true> tetFile;
			//		model::UniqueOutputFile<'T'> tetFile;
			
			for (TetContainerType::const_iterator iter=tetContainer.begin();iter!=tetContainer.end();++iter){
				for (unsigned int k=0;k<iter->eleNodes.size()-1;++k){
					for (unsigned int j=k+1;j<iter->eleNodes.size();++j){
						tetFile<< iter->eleNodes[k]->sID<<"	"<<iter->eleNodes[j]->sID<<" 0"<<std::endl;
					}
				}
			}
		}
		
		//==================================================================================
		// function to remove all previous traction and displacement boundary conditions
		// called before setting new ones
		//==================================================================================
		
		void removeBoundaryConditions(){
			
			for (unsigned int i =0; i<nodeContainer.size(); i++){
				if(nodeContainer[i].isBoundaryNode){
					//nodeContainer[i].traction = VectorDim::Zero();         // set traction to be zero
					nodeContainer[i].remove_BCs();                        // set displacement BC to be false
				}
			}
			
		}
		
		
		
		//================================================================================
		//  function to set the dof of a face (the whole face)
		//================================================================================
		//	template <typename T>
		/*
		 void setFaceDof(int iface, int idof, double val)
		 {
		 bvpfe::Face* pFace;
		 std::map<size_t,bvpfe::Triangle*>::iterator it;
		 
		 idof = idof -1;
		 pFace = faceContainer.find(iface)->second;    // pointer to the targetted face
		 
		 //----- loop over the face triangles ---
		 
		 for (it= pFace->triContainer.begin() ; it != pFace->triContainer.end(); it++ )
		 {
		 //----- loop over the triangle nodes ---
		 for(unsigned int i = 0; i<it->second->Nnodes; i++)
		 {
		 //std::cout<< it->second->eleNodes[i]->sID;
		 //VectorDim uInf = pT->displacement(it->second->eleNodes[i]->P);
		 it->second->eleNodes[i]->setBC (idof , val-it->second->eleNodes[i]->uInf(idof)*0);
		 }
		 }
		 
		 
		 }*/
		//===========================================================================================
	};
	
	
	
}  //  namespace bvpfe
#endif




//===================================================================================
// function to assemble the force vector from surface traction
// we assume here uniform traction over the surface, so one quadrature point is used
//===================================================================================
/*
 void assembleTractionForce (VECTOR_double& Fm )
 {
 bvpfe::Face* pFace;
 std::map<bvpfe::Face*,std::vector<double> >::iterator itf;
 std::vector<double> tr;
 
 bvpfe::Triangle* pTri;
 std::map<size_t,bvpfe::Triangle*>::iterator itt;
 
 std::vector<double> triVec;
 triVec.resize(9,0.0e+00);     //The force vector for each triangle
 
 int ie , qi;
 
 //--------- loop over all faces subjected to traction ----------
 for (itf= facesWtraction.begin() ; itf != facesWtraction.end(); itf++ )
 {
 pFace = itf->first;
 tr = itf->second;
 
 //std::cout<< "traction " << tr[0] << " " << tr[1] << " " << tr[2] << std::endl;
 
 // ---------- loop over all triangles of the face ------------
 for (itt= pFace->triContainer.begin() ; itt != pFace->triContainer.end(); itt++ )
 {
 pTri = itt->second;
 
 triVec = pTri->getForceVec(tr);
 
 for (int ai = 0 ; ai < 3; ai++ )             // loop over the 3 surface element's nodes
 {
 ie = ai*3;
 
 for  (int in = 0 ; in < 3; in++ )   // loop over the 3 dofs
 {
 qi = pTri->eleNodes[ai]->eqnNumber(in);
 if (qi == -1) continue;
 Fm(qi) = Fm(qi) +  triVec[ie+in];
 }
 }
 }
 }
 
 
 }
 */
