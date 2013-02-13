/* This file is part of finite element solution of BVP attached with model "the Mechanics of Defects Evolution Library".
 *
 * Copyright (C) 2011 by Mamdouh Mohamed <mamdouh.s.mohamed@gmail.com>, 
 * Copyright (C) 2011 by Giacomo Po <giacomopo@gmail.com>.
 * 
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef bvpfe_tetrahedron_H_
#define bvpfe_tetrahedron_H_

#include "model/BVP/Triangle.h"
#include <vector>
#include <set>
#include <map>
#include <utility>
#include <memory>
#include <model/BVP/SearchData.h>
#include <model/DislocationDynamics/Materials/Material.h>


namespace bvpfe{
	
	/*template<typename S>
	struct PointerWrapper {
		const S* const pt;
		PointerWrapper(const S* const ptin):pt(ptin){}
	};*/
	
//	enum{dim=3};
	
	
	class Tetrahedron : public Element<3>, 
	/*               */ public model::StaticID<Tetrahedron>   {
		
		
		enum{dim=3};				
#include <model/BVP/commonTypeDefs.h>
		
		public :
		//static Coor Nprime[4][3];                       // definition of that as std::vector<std::vector> made a problem
		static Eigen::Matrix<double,dim,Tetrahedron::Nnodes> Nprime;
		static double mu , nu , lambda , c11 , c12 , c44;
		
		std::vector<int> neighbor;                     // set of neighbor Tets
		std::map<size_t,Triangle*> TetSurfTris;        // pointers to external surface triangles of the tetrahedron. The index is for the node opposite to the surf 
		//std::map<unsigned int,std::auto_ptr<Triangle> > TetSurfTris;
		
		typedef std::pair<bool,Tetrahedron*> isTetrahedronType;
		
        /* Constructor ******************************************************/
//		template<typename SharedType>
//		Tetrahedron(const SharedType* sharedPtr)
        		Tetrahedron()
		{
			//--------- set the basis functions (linear) derivative values  // do it only once
			if  (this->sID == 0)
			{
				setBasisDerivative();
				setElasticConstants();
			}
		}
		
		///////////////////////////////////////////////////////////////////
		/*Tetrahedron()
		{
			//--------- set the basis functions (linear) derivative values  // do it only once
			if  (this->sID == 0)
			{
				setBasisDerivative();
				setElasticConstants();
			}
		}
		*/
		//===================================================================
		// function to calculate the basis functions derivatives. Linear basis 
		// are used so their derivative is constant over the element
		//===================================================================
		
		void setBasisDerivative ()
		{
			Nprime.col(0) << -1.0e+00,
			                 -1.0e+00,
					 -1.0e+00;
					 
			Nprime.col(1) <<  1.0e+00,
			                  0.0e+00,
					  0.0e+00;
					  
			Nprime.col(2) <<  0.0e+00,
			                  1.0e+00,
					  0.0e+00;
					  
			Nprime.col(3) <<  0.0e+00,
			                  0.0e+00,
					  1.0e+00;

		}
		
		//===================================================================
		// function to insert new surface triangle pointer to the  TetSurfTris
		//===================================================================
		void insertSurfTri(Triangle* pTri)
		//void insertSurfTri(std::auto_ptr<Triangle> pTri)
		{
		  bool Included;
		  unsigned int indx;
		  
		  for (unsigned int i =0; i<Tetrahedron::Nnodes; i++){
		    Included = false;
		    for (unsigned int j =0; j<Triangle::Nnodes; j++){
		      if (this->eleNodes[i]->sID == pTri->eleNodes[j]->sID) {Included = true;}
		    }
		    if(!Included) {indx = i;   break;}
		  }
		  
		  TetSurfTris.insert (std::pair<size_t,Triangle*>(indx,pTri) );
		  //TetSurfTris.insert (std::pair<unsigned int,std::auto_ptr<Triangle> >(indx,pTri));
		}
		
		//===================================================================
		// function to calculate the elastic constants
		//===================================================================
		//template<typename SharedType>
		void setElasticConstants ()
		{
			//mu =this->material.mu;
			//nu =this->material.nu;
//			mu =sharedPtr->material.mu;
//			nu =sharedPtr->material.nu;
            mu =model::Material<model::Isotropic>::mu;
			nu =model::Material<model::Isotropic>::nu;
						
			lambda = (2.0e+00*mu * nu)/(1.0e+00-(2.0e+00*nu)) ;
			
			c11 = lambda + (2.0e+00*mu);
			c12 = lambda;
			c44 = mu;
		}
		
		//===================================================================
		// function to add neighbor elements to the neighbor list
		//==================================================================
		void addNeighbor(int nbr[4])
		{
			for (int i = 0; i<4; i++) {neighbor.push_back(nbr[i]);}
			
			//std::cout<<"neighbor "<< this->sID<< " : "<<nbr[0]<<" "<<nbr[1]<<" "<<nbr[2]<<" "<<nbr[3]<<" "<<std::endl;
		}
		
		//===================================================================
		// function to calculate the stiffness matric for the tetrahedron
		// [calculations are done only at one Gaussian point]
		//===================================================================
		Eigen::Matrix<double,12,12> getElementStiffness(Eigen::Matrix<double,dim,Tetrahedron::Nnodes>& Np, double & wv)
		{
			MatrixDim  Jm = computeJacobianMatrix();                // the Jacobian matrix transpose
			
			wv = Jm.determinant() / 6.0e+00;       // must mutiply by 1/6 which is the reference volume
			
			MatrixDim JmInv = Jm.inverse();                       // inverse the matrix
			
			//------------- map shape function derivatives to actual element --
						
			for (unsigned int i=0;i<4;i++){
			  Np.col(i) = JmInv*Nprime.col(i);
			}
			
			//------------- loop over nodes and form loop matrices ---------------------------
			
			Eigen::Matrix<double,12,12> tetMat = Eigen::Matrix<double,12,12>::Zero();         // 12x12 stiffness matrix for each tetrahedron
			
			Eigen::Matrix<double,dim,dim> Kn;                  // stiffness matrix for one node
			unsigned int ie , je;
			
			for ( unsigned int ai = 0 ; ai < 4 ; ai ++)
			{
				ie= ai*3;                  // index in the Ke array
				for ( unsigned int bi = 0 ; bi < 4 ; bi ++)
				{
					je = bi*3 ;       // index in the Ke array
					
					Kn = getNodeMatrix (ai,bi,Np.transpose(),wv);
					
					for ( unsigned int in = 0 ; in < 3 ; in ++)
					{
						for ( unsigned int jn = 0 ; jn < 3 ; jn ++)
						{
							tetMat(ie+in,je+jn) += Kn(in,jn);
						}
					}
				}
			}
			
			//std::cout << tetMat << std::endl;
			
			return tetMat;
		}
		
		//===================================================================
		// function to calculate the stiffness matric for the tetrahedron
		// [calculations are done only at one Gaussian point]
		//===================================================================
		Eigen::Matrix<double,12,12> getElementStiffness()
		{
			Eigen::Matrix<double,dim,Tetrahedron::Nnodes> Np;
			double  wv;
		  
			MatrixDim  Jm = computeJacobianMatrix();                // the Jacobian matrix transpose
			
			wv = Jm.determinant() / 6.0e+00;       // must mutiply by 1/6 which is the reference volume
			
			MatrixDim JmInv = Jm.inverse();                       // inverse the matrix
			
			//------------- map shape function derivatives to actual element --
						
			for (unsigned int i=0;i<4;i++){
			  Np.col(i) = JmInv*Nprime.col(i);
			}
			
			//------------- loop over nodes and form loop matrices ---------------------------
			
			Eigen::Matrix<double,12,12> tetMat = Eigen::Matrix<double,12,12>::Zero();         // 12x12 stiffness matrix for each tetrahedron
			
			Eigen::Matrix<double,dim,dim> Kn;                  // stiffness matrix for one node
			unsigned int ie , je;
			
			for ( unsigned int ai = 0 ; ai < 4 ; ai ++)
			{
				ie= ai*3;                  // index in the Ke array
				for ( unsigned int bi = 0 ; bi < 4 ; bi ++)
				{
					je = bi*3 ;       // index in the Ke array
					
					Kn = getNodeMatrix (ai,bi,Np.transpose(),wv);
					
					for ( unsigned int in = 0 ; in < 3 ; in ++)
					{
						for ( unsigned int jn = 0 ; jn < 3 ; jn ++)
						{
							tetMat(ie+in,je+jn) += Kn(in,jn);
						}
					}
				}
			}
			
			//std::cout << tetMat << std::endl;
			
			return tetMat;
		}
		
		//===================================================================
		// function to calculate Jacobian matrix
		//===================================================================
		MatrixDim computeJacobianMatrix() const
        {
			MatrixDim J;
			
			for (int c=0;c<dim;++c){
				J.col(c)= this->eleNodes[c+1]->P - this->eleNodes[0]->P;
			}
			
			return J.transpose();
		}
		
		//===================================================================
		// function to calculate Jacobian matrix
		//===================================================================
		inline Eigen::Matrix<double,dim,dim> getNodeMatrix (unsigned int ai,unsigned int bi, Eigen::Matrix<double,Tetrahedron::Nnodes,dim> phi, double wv)
		{
			Eigen::Matrix<double,dim,dim> Kn;
						
			Kn(0,0)= wv *((phi(ai,0)*c11*phi(bi,0)) + (phi(ai,1)*c44*phi(bi,1)) +  
						   (phi(ai,2)*c44*phi(bi,2)));
			Kn(0,1)= wv *((phi(ai,0)*c12*phi(bi,1)) + (phi(ai,1)*c44*phi(bi,0)) );
			Kn(0,2)= wv *((phi(ai,0)*c12*phi(bi,2)) + (phi(ai,2)*c44*phi(bi,0)) );
			
			Kn(1,0)= wv *((phi(ai,1)*c12*phi(bi,0)) + (phi(ai,0)*c44*phi(bi,1)) );
			Kn(1,1)= wv *((phi(ai,1)*c11*phi(bi,1)) + (phi(ai,0)*c44*phi(bi,0)) + 
						   (phi(ai,2)*c44*phi(bi,2)) );
			Kn(1,2)= wv *((phi(ai,1)*c12*phi(bi,2)) + (phi(ai,2)*c44*phi(bi,1)) );
			
			Kn(2,0)= wv *((phi(ai,2)*c12*phi(bi,0)) + (phi(ai,0)*c44*phi(bi,2)) );
			Kn(2,1)= wv *((phi(ai,2)*c12*phi(bi,1)) + (phi(ai,1)*c44*phi(bi,2)) );
			Kn(2,2)= wv *((phi(ai,2)*c11*phi(bi,2)) + (phi(ai,1)*c44*phi(bi,1)) +  
						   (phi(ai,0)*c44*phi(bi,0)) );
			
			return Kn;
		}
		
		//===================================================================
		// function to calculate stress inside the tetrahedron. Also One Gauss point is used
		//===================================================================
		
		Eigen::Matrix<double,dim,dim> getStress() const 
		{						
			Eigen::Matrix<double,dim,dim> uprim = getUprim();
						
			//------- The stress components ---------------
			Eigen::Matrix<double,dim,dim> stress;	
			
			double em = lambda*uprim.trace();
			
			stress(0,0) = em + (2.0e0*mu*uprim(0,0));
			stress(1,1) = em + (2.0e0*mu*uprim(1,1));
			stress(2,2) = em + (2.0e0*mu*uprim(2,2));
			
			stress(0,1) = mu*(uprim(0,1)+uprim(1,0));               stress(1,0) = stress(0,1);
			stress(0,2) = mu*(uprim(2,0)+uprim(0,2));               stress(2,0) = stress(0,2);
			stress(1,2) = mu*(uprim(1,2)+uprim(2,1));               stress(2,1) = stress(1,2);
			
			
			return stress;
			
		}
		
		//===================================================================
		// function to calculate elastic lattice rotation inside the tetrahedron. Also One Gauss point is used
		//===================================================================
		
		Eigen::Matrix<double,dim,dim> getLatticeRotation() const 
		{						
			Eigen::Matrix<double,dim,dim> uprim = getUprim();
						
			return 0.5e00*(uprim.transpose() - uprim);
			
		}
		
		//===================================================================
		// function to calculate derivative of the displacement inside the tet
		//===================================================================
		
		Eigen::Matrix<double,dim,dim> getUprim() const
		{
			
			Eigen::Matrix<double,dim,dim> uprim = Eigen::Matrix<double,dim,dim>::Zero();
							
			MatrixDim  Jm = computeJacobianMatrix();                // the Jacobian matrix transpose
			
			MatrixDim  JmInv = Jm.inverse();                       // inverse the matrix
			
			//------------- map shape function derivatives to actual element --
			Eigen::Matrix<double,dim,Tetrahedron::Nnodes> Np;
			
			for (unsigned int i=0;i<4;i++){
			  Np.col(i) = JmInv*Nprime.col(i);
			}
			
			for(unsigned int  iu = 0 ; iu< 3; iu++)
			{
				for(unsigned int ix = 0 ; ix< 3; ix++)
				{
					for(unsigned int i = 0 ; i< 4; i++)
					{						
						uprim (iu,ix) += (Np(ix,i) * this->eleNodes[i]->u(iu) );
					}
				}
			}
			
			return uprim;
		}
		
		//=====================================================================
		// function to return check if any point P is inside the Tet (return -10 in this case).
		// otherwise return the index of the neighbor Tet in the direction of the point
		//========================================================================
		
		int nextNeighbor(VectorDim P, bool & found , int & ni)
		{	
			Eigen::Matrix<double,4,1>::Index ii;
			double baryMin = getBarycentric(P).minCoeff(&ii);
			ni = ii;
			//if (baryMin>=0.0) found = true; 
			if ((baryMin >= -1.0e-8 && neighbor[ii]==-1)||(baryMin>= 0.0)) found = true; 
			
			return neighbor[ii];
		}
		
		
		//=====================================================================
		// function to return the position of the dislocation node in the FE mesh
		// or return the ID for the next search tetrahedron
		//========================================================================
		/*
		void nextNeighbor(VectorDim P, bool & found , int & ni)
		{	
			Eigen::Matrix<double,4,1>::Index ii;
			double baryMin = getBarycentric(P).minCoeff(&ii);
			ni = ii;
			//if (baryMin>=0.0) found = true; 
			if ((baryMin >= -1.0e-3 && neighbor[ii]!=-1)||(baryMin> 0.0)) found = true; 
			
			return neighbor[ii];
		}
		*/
//		//=======================================================================
//		// function to return an approximation for the position of a point (that
//		// was found outside the mesh) to be repositioned on the mesh surface corresponding to node in
//		//========================================================================
//		
//		VectorDim findApproxSurfPoint(const VectorDim& P_in, const int& ni) const {
//		  
//		  Eigen::Matrix<double,4,1> Bary = Eigen::Matrix<double,4,1>::Zero();
//		  
////		  for(unsigned int i=0; i<4; i++){ Giacomo 09-30-2011
//			  for(int i=0; i<4; i++){
//		    if (i != ni) {
//		      Bary(i) = 1.0/((this->eleNodes[i]->P - P_in).norm());
//		    }
//		  }
//		  
//		  Eigen::Matrix<double,4,1> B = Bary/Bary.sum();   // correct the summation to be one
//		  
//		  VectorDim P = VectorDim::Zero();
//		  
//		  for(unsigned int i=0; i<4; i++) {
//		    P += (B(i)*this->eleNodes[i]->P);
//		  } 
//		  
//		  return P;
//		}
		
		//=======================================================================
		// function to return a point in actual domain given its barycentric coordinates
		//========================================================================
		
		VectorDim mapBary2Actual (Eigen::Matrix<double,4,1>  Bary){
		  
		  Eigen::Matrix<double,4,1> B = Bary/Bary.sum();   // insure the summation to be one
		  
		  VectorDim P = VectorDim::Zero();
		  
		  for(unsigned int i=0; i<4; i++) {
		    P += (B(i)*this->eleNodes[i]->P);
		  } 
		  
		  return P;
		}
		
		
		//=====================================================================
		// function to check if a point is inside a tetrahedron
		//=====================================================================
		
		bool isInsideTet(model::SearchData<dim> & data, int & sI, Eigen::Matrix<double,4,1> & Bary) {

			// Find minimum baricentric coordinate and corresponding face			
			Eigen::Matrix<double,4,1>::Index ii;
			Bary = 	getBarycentric(data.P);		
			
			double baryMin = Bary.minCoeff(&ii);	
			bool isInsideT = (baryMin >= 0.0e00 && neighbor[ii]!=-1)||(baryMin> 0.0e00);    //completely inside, or on the edge but must be interior edge
			
			//bool foundOnDomainBoundary = baryMin == 0.0e00 && neighbor[ii]==-1;
			bool foundOnDomainBoundary = baryMin <= 1.0e-12 && baryMin >= -1.0e-12 && neighbor[ii]==-1;
			
			//std::cout << baryMin << " " << neighbor[ii] << " " << isInsideT << " " << foundOnDomainBoundary << std::endl;
			//bool isInsideT = baryMin>=0.0;

			
			bool isOutsideD = neighbor[ii]==-1;
			
			if(isInsideT){            // point is inside a tetrahedron
			  data.found=true;
			  //data.nodeMeshLocation = 1;     // NOT NEEDED
			  data.outwardFaceNormal=VectorDim::Zero();     // no surface normal  (NEEDS RETHINKING FOR A MOVING SURFACE NODE)
			}
			
			else if(isOutsideD && !foundOnDomainBoundary){
			  data.nodeMeshLocation = 0;       // point is currently outside the domain, this has to be changed later to be boundary
			  sI = ii;
			}
			
			else if (foundOnDomainBoundary) {
			  data.nodeMeshLocation = 2;
			  data.found=true;
			  data.outwardFaceNormal= TetSurfTris.find(ii)->second->outNormal;
			  data.triIndex = TetSurfTris.find(ii)->second->sID;
			  data.projectedP = data.P;
			}
			
			else{                          // continue looking in another tetrahedron
			  data.found=false;
			  data.nodeMeshLocation = 1;      
			  data.outwardFaceNormal=VectorDim::Zero();     // NO NEED FOR THAT
			  data.newMeshID = neighbor[ii];
			}
			
			return isInsideT;		  
		}
		
		
		//=====================================================================
		// function to check if a point is inside a tetrahedron
		//=====================================================================
		
		bool isInsideTet(model::SearchData<dim> & data ) {

			// Find minimum baricentric coordinate and corresponding face			
			Eigen::Matrix<double,4,1>::Index ii;
			Eigen::Matrix<double,4,1> Bary = 	getBarycentric(data.P);		
			
			double baryMin = Bary.minCoeff(&ii);	
			bool isInsideT = (baryMin >= 0.0e00 && neighbor[ii]!=-1)||(baryMin> 0.0e00);    //completely inside, or on the edge but must be interior edge
			//bool foundOnDomainBoundary = baryMin == 0.0e00 && neighbor[ii]==-1;
			bool foundOnDomainBoundary = (baryMin <= 1.0e-7) && (baryMin >= -1.0e-7) && (neighbor[ii]==-1);
			//std::cout << baryMin << " " << neighbor[ii] << " " << isInsideT << " " << foundOnDomainBoundary << std::endl;
			//bool isInsideT = baryMin>=0.0;

			
			bool isOutsideD = neighbor[ii]==-1;
			
			if(isInsideT && !foundOnDomainBoundary){            // point is inside a tetrahedron
			  data.found=true;
			  data.nodeMeshLocation = 1;     // NOT NEEDED
			  data.outwardFaceNormal=VectorDim::Zero();     // no surface normal  (NEEDS RETHINKING FOR A MOVING SURFACE NODE)
			}
			
			else if(isOutsideD && !foundOnDomainBoundary){
			  data.nodeMeshLocation = 0;       // point is currently outside the domain, this has to be changed later to be boundary
			  //int sI = ii;
			}
			
			else if (foundOnDomainBoundary) {
			  data.nodeMeshLocation = 2;
			  data.found=true;
			  data.outwardFaceNormal= TetSurfTris.find(ii)->second->outNormal;
			  data.triIndex = TetSurfTris.find(ii)->second->sID;
			  data.projectedP = data.P;
			}
			
			else{                          // continue looking in another tetrahedron
			  data.found=false;
			  //data.nodeMeshLocation = 1;      
			  data.outwardFaceNormal=VectorDim::Zero();     // NO NEED FOR THAT
			  data.newMeshID = neighbor[ii];
			}
			
			return isInsideT;		  
		}
		
		
		
		//===================================================================================
		// function to calculate the unit normal of the tetrahedron face facing node number i
		//===================================================================================
		
		//----- This function became USELESS. get instead TetSurfTris.find(ii)->second->outNormal
		
		VectorDim getElementOutNormal(unsigned int faceID) {
		  
		  /* nodes are ordered in a way that if you move over the first 3 nodes 
		   with the right hand rule, you will point to the 4th one */
		  
		  /* this way of implementation does not depend on the nodes order  */

		Eigen::Matrix<double,dim,1> org = this->eleNodes[faceID]->P;
		Eigen::Matrix<double,dim,1> v1 = this->eleNodes[(faceID+2)%4]->P - this->eleNodes[(faceID+1)%4]->P;
		Eigen::Matrix<double,dim,1> v2 = this->eleNodes[(faceID+3)%4]->P - this->eleNodes[(faceID+1)%4]->P;
		Eigen::Matrix<double,dim,1> v3 = org - this->eleNodes[(faceID+1)%4]->P;
		  
		Eigen::Matrix<double,dim,1> Nv = v1.cross(v2).normalized();    //Nv = Nv.normalized();
		  
		if(Nv.dot(v3) > 0.0) Nv = -Nv;
			
//			std::cout<<"Tet "<<"I'm here 6"<<std::endl;
		  return Nv;
		}

		//=====================================================================
		// function to return the barycentric coordinates of any point P
		//========================================================================
		
		Eigen::Matrix<double,4,1>  getBarycentric (VectorDim P)
		{
			//Vector bary;            bary.resize(4);
			Eigen::Matrix<double,4,1> bary;
			double V0;
			
			V0 = getVol(this->eleNodes[0]->P,this->eleNodes[1]->P,this->eleNodes[2]->P,this->eleNodes[3]->P);
			
			bary(0) = getVol(                   P,this->eleNodes[1]->P,this->eleNodes[2]->P,this->eleNodes[3]->P)/V0;
			bary(1) = getVol(this->eleNodes[0]->P,                   P,this->eleNodes[2]->P,this->eleNodes[3]->P)/V0;
			bary(2) = getVol(this->eleNodes[0]->P,this->eleNodes[1]->P,                   P,this->eleNodes[3]->P)/V0;
			bary(3) = getVol(this->eleNodes[0]->P,this->eleNodes[1]->P,this->eleNodes[2]->P,                   P)/V0;
			
			return bary;
		}
		
		//======================================================================
		// return the volume of a tetradedron, given its points coordinates in order
		//======================================================================
		double getVol(VectorDim a, VectorDim b, VectorDim c, VectorDim d)
		{

			MatrixDim temp;
			
			temp.col(0) = a - d;
			temp.col(1) = b - d;
			temp.col(2) = c - d;
			
			return temp.determinant()/6.0;
		}
		
		//======================================================================
		// return the volume of the current mesh tetradedron
		//======================================================================
		double getTetVolume(){
		  return getVol(this->eleNodes[0]->P,this->eleNodes[1]->P,this->eleNodes[2]->P,this->eleNodes[3]->P);
		}
		
		//=====================================================================
		// function to map a given point (i) to the reference coordinate system
		//=====================================================================
		
		template <short unsigned int qOrder>				
		Eigen::Matrix<double,dim,qOrder> mapGaussPoints(){
			
			//-----------------mapping constants ------------------------
			MatrixDim A;
			VectorDim b = this->eleNodes[0]->P;

			for (int i = 0; i<dim; i++){
				A.col(i) = this->eleNodes[i+1]->P - this->eleNodes[0]->P;
			}
			//-----------------map the Gauss points ------------------------
			// gp = A*P + b
					
			Eigen::Matrix<double,dim,qOrder> gp = A * model::Quadrature<dim,qOrder>::abcsissas_  ; // FIX THIS !!!! add vector to matrix

			for (int i=0; i < qOrder; i++){
				gp.col(i)+=b;
			}
			
		    return gp;
		}
		

		
		//==============================================================================================
		// Function to integrate the infinite medium stress field over Tetrahedron element
		//==============================================================================================
		
		template <short unsigned int qOrder, typename T>
		Eigen::Matrix<double,dim,3> getTetInfiniteStress (const T* const pt) {
		  
		  Eigen::Matrix<double,dim,3> stressInt=Eigen::Matrix<double,dim,3>::Zero();
		  
		  PointerWrapper<T> pts(pt);
		  		  
		  model::Quadrature<3,qOrder>::integrate(this,stressInt,&Tetrahedron::dislocationStressKernel<PointerWrapper<T> > , pts);
		  
		  return stressInt;
		}
		
		//============================================================================
		// function to return the integration kernel for the infinite medium stress field
		//=============================================================================
		
		template<typename T> 
		MatrixDim dislocationStressKernel(const VectorDim& Rstd, const T& pts) const {
			double J;
			VectorDim R=mapStdPoint(Rstd,J);
			return pts.pt->stress(R)*6.0e00;           // this is averaging over the volume, NOT integration over the volume
			                                           // and the sum of intergration weights = (1/6)
		}
		
		
		
		
		//==============================================================================================
		// Function to integrate the infinite medium lattice rotation field over Tetrahedron element
		//==============================================================================================
		
		template <short unsigned int qOrder, typename T>
		Eigen::Matrix<double,dim,3> getTetInfiniteRotation (const T* const pt) {
		  
		  Eigen::Matrix<double,dim,3> rotationInt=Eigen::Matrix<double,dim,3>::Zero();
		  
		  PointerWrapper<T> pts(pt);
		  		  
		  model::Quadrature<3,qOrder>::integrate(this,rotationInt,&Tetrahedron::dislocationRotationKernel<PointerWrapper<T> > , pts);
		  
		  return rotationInt;
		}
		
		//============================================================================
		// function to return the integration kernel for the infinite medium stress field
		//=============================================================================
		
		template<typename T> 
		MatrixDim dislocationRotationKernel(const VectorDim& Rstd, const T& pts) const {
			double J;
			VectorDim R=mapStdPoint(Rstd,J);
			return pts.pt->lattice_rotation(R)*6.0e00;           // this is averaging over the volume, NOT integration over the volume
			                                           // and the sum of intergration weights = (1/6)
		}
		
		//==============================================================================================
		// Function to integrate the infinite medium displacement gradient field over Tetrahedron element
		//==============================================================================================
		
		template <short unsigned int qOrder, typename T>
		Eigen::Matrix<double,dim,3> getTetInfiniteUprim (const T* const pt) {
		  
		  Eigen::Matrix<double,dim,3> uprimInt=Eigen::Matrix<double,dim,3>::Zero();
		  
		  PointerWrapper<T> pts(pt);
		  		  
		  model::Quadrature<3,qOrder>::integrate(this,uprimInt,&Tetrahedron::dislocationUprimKernel<PointerWrapper<T> > , pts);
		  
		  return uprimInt;
		}
		
		//============================================================================
		// function to return the integration kernel for the infinite medium stress field
		//=============================================================================
		
		template<typename T> 
		MatrixDim dislocationUprimKernel(const VectorDim& Rstd, const T& pts) const {
			double J;
			VectorDim R=mapStdPoint(Rstd,J);
			return pts.pt->displacement_gradient(R)*6.0e00;           // this is averaging over the volume, NOT integration over the volume
			                                           // and the sum of intergration weights = (1/6)
		}
		
		/***************************************************************************/
		VectorDim mapStdPoint(const VectorDim& Rstd, double& J) const {
			MatrixDim A;
			VectorDim b = this->eleNodes[0]->P;
		    for (int i = 0; i<dim; i++){
				A.col(i) = this->eleNodes[i+1]->P - this->eleNodes[0]->P;
		    }
			J=A.determinant();
			assert(J>0.0);
			return  A * Rstd +b ; 
		}
		
		//double elementJacobian
		
		
	};
	
	
	
	//==========================================================================
	//=========================================================================
	
	//------------ redefine the static variable [required] ------------
	Eigen::Matrix<double,3,Tetrahedron::Nnodes> Tetrahedron::Nprime;
	
	double Tetrahedron::mu ;
	double Tetrahedron::nu ;
	double Tetrahedron::lambda;
	
	double Tetrahedron::c11 ;
	double Tetrahedron::c12 ;
	double Tetrahedron::c44 ;
	
	//	bvpfe::QuadratureTet<4> Tetrahedron::Quad;
	
}  //  namespace bvpfe
#endif


/*********************************************************************/

/*Curently USELESS*/
/*
 template <short unsigned int qOrder, typename T>
 Eigen::Matrix<double,dim,Element<3>::Nnodes> get_dislocationForceVector(const T* const pt, std::vector<std::vector<double> > const Np) const{
 
 Eigen::Matrix<double,dim,dim> stressInt=Eigen::Matrix<double,dim,dim>::Zero();
 
 PointerWrapper<T> pts(pt);
 //			pts.pt=pt;
 
 model::Quadrature<dim,qOrder>::integrate(this,stressInt,&Tetrahedron::dislocationStressKernel<PointerWrapper<T> > ,pts);
 
 //-------multiply stress by the shape function derivative --------------------------
 
 Eigen::Matrix<double,dim,Element<3>::Nnodes> temp=Eigen::Matrix<double,dim,Element<3>::Nnodes>::Zero();
 
 for (unsigned int i = 0; i< Tetrahedron::Nnodes; ++i){
 temp(0,i) = (Np[i][0]*stressInt(0,0))+(Np[i][1]*stressInt(1,1))+(Np[i][2]*stressInt(2,2));
 temp(1,i) = (Np[i][1]*stressInt(0,1))+(Np[i][0]*stressInt(1,1))+(Np[i][2]*stressInt(1,2));
 temp(2,i) = (Np[i][2]*stressInt(0,2))+(Np[i][1]*stressInt(1,2))+(Np[i][0]*stressInt(2,2));
 }
 
 return temp;
 }*/

/*********************************************************************/
