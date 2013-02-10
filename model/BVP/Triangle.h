/* This file is part of finite element solution of BVP attached with model "the Mechanics Of Defects Evolution Library".
 *
 * Copyright (C) 2011 by Mamdouh Mohamed <mamdouh.s.mohamed@gmail.com>, 
 * Copyright (C) 2011 by Giacomo Po <giacomopo@gmail.com>.
 * 
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef bvpfe_triangle_H_
#define bvpfe_triangle_H_

#include "model/Utilities/StaticID.h"
#include <model/DislocationDynamics/GlidePlanes/GlidePlaneObserver.h>
#include "model/BVP/Element.h"
#include <model/Utilities/CompareVectorsByComponent.h>
#include <cmath>
#include <vector>
#include<Eigen/StdVector>
//#include "model/BVP/Tetrahedron.h"

//#include "model/Quadrature/Quadrature.h"

namespace bvpfe {
  
	template<typename S>
	struct PointerWrapper {
		const S* const pt;
		PointerWrapper(const S* const ptin):pt(ptin){}
	};
  

	class Triangle    : public bvpfe::Element<2>, 
		            public model::StaticID<Triangle>  

	{	 

		public :

		  enum{dim=3};
		  
#include <model/BVP/commonTypeDefs.h>
		  
		  Eigen::Matrix<double,dim,3> TriVec; 

		//--- the value of the 3 shape functions at three Gauss points (1/3, 1/3)
		static double N[3][3];          // definition of that as std::vector<double> made a problem
		size_t neighTetIndx;        // index to the neighbor tetrahedron for this surface triangle
		
		std::vector<Triangle*> neighbor;   // pointers to the 3 neighbor triangles to this triangle
		
		std::vector<std::pair<unsigned int,unsigned int> > cuttingSegments;   // pair nodes for dislocation segments that cut the triangle
		
		typedef std::vector< VectorDim , Eigen::aligned_allocator<VectorDim> > vectorDimVectorType;		
		
		typedef std::map<Eigen::Matrix<double,dim+1,1>, vectorDimVectorType ,  
		                 model::CompareVectorsByComponent<double,dim+1,float>,
				 Eigen::aligned_allocator<std::pair<const Eigen::Matrix<double,dim+1,1>,vectorDimVectorType > > > localQuadraturePointsContainerType;
		
		
		localQuadraturePointsContainerType localQuadPnts;     // stores the customly generated Quadrature points around each glide plane that intersects the triangle
		                                                                   // , where the key is the glide plane sID
		                                                                   
		//std::map<Eigen::Matrix<double,dim+1,1>, std::vector<VectorDim> , model::CompareVectorsByComponent<double,dim+1,float> > localQuadPnts;
		                                                                   
		                                                                   
		//std::map<Eigen::Matrix<double,dim+1,1>, std::vector<VectorDim> > localQuadPnts;     // stores the customly generated Quadrature points around each glide plane that intersects the triangle
		                                                                                    // , where the key is vector that contains the glide plane normal and hight
										   
		//std::map<unsigned int, double > localQuadPnts_w;                    // the associated weight with each quadrature points set. Weights are equal for all points, and they have to sum to 0.5
		
		Eigen::Matrix<unsigned int,2,1> projPlaneIndx;                        // the index for the projection plane for this triangle
		
		VectorDim outNormal;       // The outward normal to the surface triangle
		
		static double edgeTol;

		Triangle()
		{
			//--------- set the basis functions (linear) derivative values  // do it only once
			if  (this->sID == 0)   setBasis();
			neighbor.resize(3);
			
			//---------- set the projection plane for this triangle -----
			//projPlaneIndx= find2DProjectionPlane();
		}		

		//===================================================================
		// function to calculate the basis functions at 3 Gauss point (1/6, 2/3) . Linear basisare used 
		//===================================================================

		void setBasis ()
		{
			double twothird , onesixth;
			
			twothird = 2.0e+00/ 3.0e+00;
			onesixth = 1.0e+00/ 6.0e+00;

			N[0][0] = 1.0e+00-twothird-onesixth ;  N[0][1] = twothird;   N[0][2] = onesixth;
			N[1][0] = 1.0e+00-twothird-onesixth ;  N[1][1] = onesixth;   N[1][2] = twothird;
			N[2][0] = 1.0e+00-onesixth-onesixth ;  N[2][1] = onesixth;   N[2][2] = onesixth;

			//std::cout<< N[0] <<" " << N[1]<< " " << N [2] << std::endl;
 			
		}
				
		//============================================================================
		// function to return area of the actual triangle (used in Gauss integration)
		//=============================================================================

		double area() const
		{
			double ar;
			double a[3] , b[3];

			a[0] = this->eleNodes[1]->P(0) - this->eleNodes[0]->P(0);
			a[1] = this->eleNodes[1]->P(1) - this->eleNodes[0]->P(1);
			a[2] = this->eleNodes[1]->P(2) - this->eleNodes[0]->P(2);

			b[0] = this->eleNodes[2]->P(0) - this->eleNodes[0]->P(0);
			b[1] = this->eleNodes[2]->P(1) - this->eleNodes[0]->P(1);
			b[2] = this->eleNodes[2]->P(2) - this->eleNodes[0]->P(2);
		
			ar = pow((a[1]*b[2])-(a[2]*b[1]),2) + pow((a[2]*b[0])-(a[0]*b[2]),2) +  pow((a[0]*b[1])-(a[1]*b[0]),2);

			ar = sqrt(ar)/ 2.0e+00;

			//std::cout<< "area " << ar << std::endl;
			return ar;
		}
		
		//============================================================================
		// function to return the current (deformed) area of the triangle element (used in indentation force calsulation)
		//=============================================================================

		double deformed_area() const
		{
			double ar;
			double a[3] , b[3];

			a[0] = this->eleNodes[1]->Pc(0) - this->eleNodes[0]->Pc(0);
			a[1] = this->eleNodes[1]->Pc(1) - this->eleNodes[0]->Pc(1);
			a[2] = this->eleNodes[1]->Pc(2) - this->eleNodes[0]->Pc(2);

			b[0] = this->eleNodes[2]->Pc(0) - this->eleNodes[0]->Pc(0);
			b[1] = this->eleNodes[2]->Pc(1) - this->eleNodes[0]->Pc(1);
			b[2] = this->eleNodes[2]->Pc(2) - this->eleNodes[0]->Pc(2);
		
			ar = pow((a[1]*b[2])-(a[2]*b[1]),2) + pow((a[2]*b[0])-(a[0]*b[2]),2) +  pow((a[0]*b[1])-(a[1]*b[0]),2);

			ar = sqrt(ar)/ 2.0e+00;

			//std::cout<< "area " << ar << std::endl;
			return ar;
		}


		//============================================================================
		// function to return the triangle force vector resulted from surface traction
		//=============================================================================

		std::vector<double> getForceVec(std::vector<double> tr)
		{
			std::vector<double> triVec;      triVec.resize(9,0.0e+00);
			double wa;
			int ie;

			//wa = (1.0e+00/3.0e+00)*area()/0.5e+00;    // this works as the det of the Jacobian  = area of actual/ area of reference
			wa = (1.0e+00/3.0e+00)*area();   //this is the right one because the integration wieght = (1/3) * (1/2) "reference area"
			//std::cout<< "wa ====" << wa << std::endl;

			
			for (int ai = 0 ; ai < 3; ai++ )             // loop over the 3 triangle nodes
			{
				ie = ai*3;

				for  (int jj = 0 ; jj < 3; jj++ )   // loop over the 3 dofs
				{
					for (int ii = 0 ; ii < 3; ii++ )             // loop over 3 Gauss points
					{		
					triVec[ie+jj] = triVec[ie+jj] + ( wa * N[ii][ai] * tr[jj]); 
					}
				}
			}			
			
			return triVec;
		}
		
		//============================================================================
		// function to return the triangle force vector resulted from infinite medium surface traction
		//=============================================================================
		template <short unsigned int qOrder, typename T>
		Eigen::Matrix<double,dim,3> getTriInfiniteForce (const T* const pt) {
		  
		  Eigen::Matrix<double,dim,3> tractionInt=Eigen::Matrix<double,dim,3>::Zero();
		  
		  PointerWrapper<T> pts(pt);
		  		  
		  Eigen::Matrix<double,3,3> tractionMatrix;
		  for (unsigned int i = 0; i<3; i++){tractionMatrix.col(i) = this->eleNodes[i]->traction;}

		  model::Quadrature<2,qOrder>::integrate(this,tractionInt,&Triangle::dislocationStressKernel<PointerWrapper<T> > , pts, tractionMatrix);
		  
		  return tractionInt;
		}
		
		
		//============================================================================
		// infinite medium stress integration kernel
		//=============================================================================
		template<typename T> 
		Eigen::Matrix<double,dim,3> dislocationStressKernel(const Eigen::Matrix<double,2,1>& Rstd, const T& pts, const Eigen::Matrix<double,3,3>& tractionMatrix) const {
			double J ;
			VectorDim R=mapStdPoint<2>(Rstd,J);
			
			//---------- interpolate externally applied traction between nodes ---------------
			VectorDim externalTraction = VectorDim::Zero();
			Eigen::Matrix<double,1,3> shapeFunc = Shapefunc(Rstd);
			for (unsigned int i = 0; i<3; i++){externalTraction+=shapeFunc(i)*tractionMatrix.col(i);}
			
			//Eigen::Matrix<double,dim,dim> triStress = pts.pt->stress(R,this->sID);
			Eigen::Matrix<double,dim,dim> triStress = pts.pt->stress(R);
			
			return (externalTraction-(triStress*outNormal))*shapeFunc*J ;
			//return (externalTraction-(pts.pt->stress(R,this->sID)*triN))*shapeFunc*J ;
		
		}
		
		//==============================================================================
		// function to map the Gauss point to the actual domain
		//==============================================================================
		template<short unsigned int OtherDim>
		VectorDim mapStdPoint(const Eigen::Matrix<double,2,1>& Rstd, double& J) const {
			
		  Eigen::Matrix<double,dim,OtherDim> A;
		  VectorDim b = this->eleNodes[0]->P;
		  
		  Eigen::Matrix<double,dim,1> temp;
		  
		    for (int i = 0; i<OtherDim; i++){
				A.col(i) = this->eleNodes[i+1]->P - this->eleNodes[0]->P;
		    }
		
		    temp = A.col(0).cross(A.col(1));
		    J= temp.norm();    // this is = area of actual triangle * 2 , but the sum of integration weights for triangle is actually in this code (1/2)               
		    
		    assert(J>0.0);
		    return  A * Rstd +b ; 
		}
		
		//================================================================================
		// function to return the unit vector normal to the triangle
		//=================================================================================
		
		VectorDim triNormal ()
		{
		  VectorDim n, v1 , v2;
		  
		  v1 = this->eleNodes[1]->P - this->eleNodes[0]->P;
		  v2 = this->eleNodes[2]->P - this->eleNodes[0]->P;
		  
		  n = v2.cross(v1);
		  
		  return n.normalized();
		}
		
		//================================================================================
		// function to return the unit vector normal to the triangle in deformed configuration
		//=================================================================================
		
		VectorDim triNormalDeformed() const
		{
		  VectorDim n, v1 , v2;
		  
		  v1 = this->eleNodes[1]->Pc - this->eleNodes[0]->Pc;
		  v2 = this->eleNodes[2]->Pc - this->eleNodes[0]->Pc;
		  
		  n = v2.cross(v1);
		  
		  return n.normalized();
		}
		//================================================================================
		// function to return the 3 shape functions at any reference point
		//=================================================================================
		
		Eigen::Matrix<double,1,3> Shapefunc (Eigen::Matrix<double,2,1> Rstd) const {
		  Eigen::Matrix<double,1,3> sf;
		  
		      sf(0) = 1.0 - Rstd(0) - Rstd(1);
		      sf(1) = Rstd(0);
		      sf(2) = Rstd(1);
		        
		  return sf;
		}
		
		//================================================================================
		// function to return the shape function for node (in) at any reference point
		//=================================================================================
		
		double Shapefunc_i (Eigen::Matrix<double,2,1> Rstd , unsigned int in)
		{
		  double sf;
		  
		  switch (in){
		    case 0:
		      sf = 1.0 - Rstd(0) - Rstd(1);    break;
		    case 1:
		      sf = Rstd(0);                    break;
		    case 2:
		      sf = Rstd(1);                    break;
		  }		  
		  return sf;
		}
		
		//=====================================================================================================
		// function to check if a given line lies on one of the triangle edges or not
		//=====================================================================================================
		
		bool onTriEdge (std::pair<VectorDim,VectorDim> lineP) {
		  VectorDim p1 = lineP.first;
		  VectorDim p2 = lineP.second;
		  
		  double tol  = edgeTol;
		  
		  return ( ((p1-this->eleNodes[0]->P).norm() < tol || (p1-this->eleNodes[1]->P).norm() < tol || (p1-this->eleNodes[2]->P).norm() < tol) && 
		           ((p2-this->eleNodes[0]->P).norm() < tol || (p2-this->eleNodes[1]->P).norm() < tol || (p2-this->eleNodes[2]->P).norm() < tol)    );
		  
		}
		
		
		//======================================================================================================
		// function to find the 2D projection plane for the triangle. On this plane, the Barycentric coordinates
		// in 2D will be calculated. If the returned values ix = 2 , iy = 3 so the projection plane is y-z plane
		//=======================================================================================================
		Eigen::Matrix<unsigned int,2,1> find2DProjectionPlane() const{
		  
		  VectorDim uV;
		  std::vector<int> normal , not_normal;
		  Eigen::Matrix<unsigned int,2,1> ppindx;
		  
		  for(int i=0; i<3; i++){
		    uV = VectorDim::Zero();     uV(i) = 1.0;
		    if(std::abs(uV.dot(outNormal)) < 0.15e00){     // the plane normal makes angle of > 80 degrees with direction x_i
		      normal.push_back(i);
		    }
		    else{not_normal.push_back(i);}
		  }
		  
		  switch (normal.size()){
		    case 2:
		      ppindx(0) = normal[0];
		      ppindx(1) = normal[1];
		      break;
		    case 1:
		      ppindx(0) = normal[0];
		      ppindx(1) = not_normal[0];
		      break;
		    case 0:
		      ppindx(0) = not_normal[0];
		      ppindx(1) = not_normal[1];
		      break;
		  }
		  return ppindx;
		}
		
		//================================================================================
		// function to find the intersection between triangle edges and a line (if any)
		//================================================================================
		
		bool intersectWithLine (VectorDim p0 , VectorDim p1 , VectorDim & P ,unsigned int & indx ){
		  
		  //------------- find plane of projection during intersection calculation -------
		  
		  //VectorDim uV;
		  //double tol = 1.0e-08;
		  //int ix , iy;
		  //ix = 0;   iy = 0;
		  
		  Eigen::Matrix<unsigned int,2,1> ppindx= find2DProjectionPlane();
		  unsigned int ix = ppindx(0);
		  unsigned int iy = ppindx(1);
		  /*
		  std::vector<int> normal , not_normal;
		  
		  for(int i=0; i<3; i++){
		    uV = VectorDim::Zero();     uV(i) = 1.0;
		    if(std::abs(uV.dot(outNormal)) < 0.15e00){     // the plane normal makes angle of > 80 degrees with direction x_i
		      normal.push_back(i);
		    }
		    else{not_normal.push_back(i);}
		  }
		  
		  switch (normal.size()){
		    case 2:
		      ix = normal[0];
		      iy = normal[1];
		      break;
		    case 1:
		      ix = normal[0];
		      iy = not_normal[0];
		      break;
		    case 0:
		      ix = not_normal[0];
		      iy = not_normal[1];
		      break;
		  }
		  */
		  //--------------- calculate barycentric coordinates of the node -----------
		  
		  bool found= false;;
		  
		  Eigen::Matrix<double,3,1>::Index ii;
		  Eigen::Matrix<double,3,1> Bary = getBarycentric(p1 , ix , iy);
		  double baryMin = Bary.minCoeff(&ii);
		  
		  if(baryMin > 0.0 && baryMin < 1.0 ){
		    //-- NO intersection, new point is on the same triangle ----------
		    //found = false;
		    P = p1;
		  }
		  else{
		    //------------------ implement intersection  with the 3 triangle edges -----------
		    VectorDim x0 , x1 , pP1;
		    Eigen::Matrix<double,2,2> mat;
		    Eigen::Matrix<double,2,1> u,b;
		    
		    std::vector<VectorDim> intersectionPoints;
		    std::vector<unsigned int> sideID;
		    VectorDim temp;
		    
		    //pP1 = p0 + 4.0e00*(p1-p0);             // extend the line to avoid the effect of any errors on the calculation
		    pP1 = p1;
		    for(unsigned int i=0; i<3; i++){               // loop over the 3 triangle edges
		      x0 = this->eleNodes[(i+1)%3]->P;
		      x1 = this->eleNodes[(i+2)%3]->P;
		    
		      mat(0,0) = pP1(ix) - p0(ix);        mat(0,1) = x0(ix) - x1(ix);    
		      mat(1,0) = pP1(iy) - p0(iy);        mat(1,1) = x0(iy) - x1(iy);
		    
		      b(0) = x0(ix) - p0(ix);
		      b(1) = x0(iy) - p0(iy);
		    
		      u = mat.inverse()*b;
		    
		      if ( /*u(0) > 1.0e-2 && u(0) <=1.0e00 &&*/ u(1) >= 0.0e00  && u(1) <=1.0e00) {       // intersection point is within the triangle line segments
			temp = x0 + u(1)*(x1-x0);
			intersectionPoints.push_back(temp);
			sideID.push_back(i);
			//P = p0 + u(0)*(pP1-p0);
			/*
			P = x0 + u(1)*(x1-x0);
			found = true;
			indx = neighbor[i]->sID;
			break;
			*/
		      }
		    }
		      
		    if(intersectionPoints.size()==0) {
		      std::cout << "Old position " <<std::setprecision(15) <<  p0.transpose() << std::endl;
		      std::cout << "New position " <<std::setprecision(15) <<  p1.transpose() << std::endl;
		      std::cout << "Traingel ID  " << this->sID << std::endl;
		      assert(0 && "Failed to calculate new position for a moving boundary node");
		    }
		    
		    else if(intersectionPoints.size()==1){
		      P = intersectionPoints[0];
		      found = true;
		      indx = neighbor[sideID[0]]->sID;
		    }
		    
		    else {    // found 2 or 3 intersection points
		      
		      //----------- check which of them is closure to the point p1
		      double minDis = 1.0e20;
		      int ii = -1;
		      for (unsigned int i1=0; i1<intersectionPoints.size();i1++){
			if ((intersectionPoints[i1]-p1).norm() < minDis){
			  ii = i1;
			  minDis = (intersectionPoints[i1]-p1).norm();
			}
		      }
		      
		      P = intersectionPoints[ii];
		      found = true;
		      indx = neighbor[sideID[ii]]->sID;
		      
		      /*
		      int repeatedPoint = -1; 
		      for (unsigned int i1=0; i1<intersectionPoints.size();i1++){
			for (unsigned int i2=i1+1; i2<intersectionPoints.size(); i2++){
			  if (std::abs(intersectionPoints[i1].dot(intersectionPoints[i2])) < 1.0e-12 ) repeatedPoint = i1; 
			}
		      }
		
		      if (repeatedPoint != -1) intersectionPoints.erase(intersectionPoints.begin()+repeatedPoint);
		      */
		      
		      
		    }
		      
		    if (!found) {
		      std::cout << "Old position " <<std::setprecision(15) <<  p0.transpose() << std::endl;
		      std::cout << "New position " <<std::setprecision(15) <<  p1.transpose() << std::endl;
		      std::cout << "Traingel ID  " << this->sID << std::endl;
		      assert(0 && "Failed to calculate new position for a moving boundary node");
		    }
		  }
		  
		  return found;
		}
		
		//==================================================================================
		// Function to split this triangle into 3 subtriangles, given the intersection line between the triangle and a segment's glide plane
		//===================================================================================
		
		std::vector<Eigen::Matrix<double,dim,dim> > spliteTriangle (std::pair<VectorDim,VectorDim> glidePlaneLine) {
		  
		  std::vector<Eigen::Matrix<double,dim,dim> > trisPointsSet;
		  Eigen::Matrix<double,dim,dim> trisPoint;
		  
		  std::vector<std::pair<VectorDim,unsigned int> >edgePoints;
		  std::vector<std::pair<VectorDim,unsigned int> >vertexPoints;
		  unsigned int in;
		  
		  if (isPointOnVertex(glidePlaneLine.first, in)) vertexPoints.push_back(std::make_pair(glidePlaneLine.first, in));
		  else edgePoints.push_back(std::make_pair(glidePlaneLine.first, in));
		  
		  if (isPointOnVertex(glidePlaneLine.second, in)) vertexPoints.push_back(std::make_pair(glidePlaneLine.second, in));
		  else edgePoints.push_back(std::make_pair(glidePlaneLine.second, in));
		    
		  // if one point on a vertex and the other on an edge that is connected to the same vertex, so the line is almost on an edge -> no splitting
		  if (edgePoints.size()==1 && vertexPoints.size()==1 && edgePoints[0].second != vertexPoints[0].second ) edgePoints.clear();   
		    
		  
		  std::vector<unsigned int> otherVertexes;
		  unsigned int otherVertex = 5;       // dummy initialization, just to avoid getting warning
		  
		  switch (edgePoints.size()){
		    case 0:          // both points are on vertexes, no spliting
		      for (unsigned int i=0;i<3;i++) trisPoint.col(i) = this->eleNodes[i]->P;
		      trisPointsSet.push_back(trisPoint);
		      break;
		      
		    case 1:         // one point on a vertex and the other point on an edge, split into only two triangles
		      
		      for (unsigned int i=0;i<3;i++) {if (i!=vertexPoints[0].second) otherVertexes.push_back(i);}
		      //---- triangle 1 ----
		      trisPoint.col(0) = this->eleNodes[vertexPoints[0].second]->P;
		      trisPoint.col(1) = edgePoints[0].first;
		      trisPoint.col(2) = this->eleNodes[otherVertexes[0]]->P;
		      trisPointsSet.push_back(trisPoint);
		      
		      //---- triangle 2 ----
		      trisPoint.col(0) = this->eleNodes[vertexPoints[0].second]->P;
		      trisPoint.col(1) = edgePoints[0].first;
		      trisPoint.col(2) = this->eleNodes[otherVertexes[1]]->P;
		      trisPointsSet.push_back(trisPoint);
		      break;
		      
		    case 2:   // two points on edges -> split into 3 triangles
		      
		      for (unsigned int i=0;i<3;i++) {if (i!=edgePoints[0].second && i!=edgePoints[1].second) otherVertex=i;}
		      
		      //---- triangle 1 ----
		      trisPoint.col(0) = edgePoints[0].first;
		      trisPoint.col(1) = edgePoints[1].first;
		      trisPoint.col(2) = this->eleNodes[otherVertex]->P;
		      trisPointsSet.push_back(trisPoint);
		      
		      //---- triangle 2 ----
		      trisPoint.col(0) = edgePoints[0].first;
		      trisPoint.col(1) = edgePoints[1].first;
		      trisPoint.col(2) = this->eleNodes[edgePoints[0].second]->P;
		      trisPointsSet.push_back(trisPoint);
		      
		      //---- triangle 2 ----
		      trisPoint.col(0) = edgePoints[0].first;
		      trisPoint.col(1) = this->eleNodes[edgePoints[0].second]->P;
		      trisPoint.col(2) = this->eleNodes[edgePoints[1].second]->P;
		      trisPointsSet.push_back(trisPoint);
		      break;
		      
		    default:
		      assert(0 && " Error in Trianlge.h::spliteTriangle function " );
		      
		  }
		 		  
		  return trisPointsSet;
		}
		
		//======================================================================================================
		// function to check if a point is on one of the triangle vertexes, and updates the index of this vertex,
		// or if it is on the edge, it update the edge index (which is the index of the corresponding vertex)
		//======================================================================================================
		
		bool isPointOnVertex (VectorDim P, unsigned int & in) {
		  bool onVertex = false;
		  
		  //----- check if the intersection point is on a vertex ----
		  for(unsigned int i=0; i<3; i++){               // loop over the 3 triangle edges
		    if ((P-this->eleNodes[i]->P).norm() < edgeTol) {
		      in= i;
		      onVertex = true;
		      break;
		    }
		  }
		  
		  if (!onVertex) in = findEdgeIndex(P);
		  
		  return onVertex;
		}
		
		//==============================================================================
		// function to return the index of the triangle edge on which the intersection point is on.
		// if the line cuts the triangle in vertex, the function returns a number > 2 (10 for example)
		//==============================================================================
		unsigned int findEdgeIndex (VectorDim P) {
		  
		  unsigned int ii;
		  VectorDim x0 , x1;

		  VectorDim crossProduct;
		  
		  for(unsigned int i=0; i<3; i++){               // loop over the 3 triangle edges
		    x0 = this->eleNodes[(i+1)%3]->P;
		    x1 = this->eleNodes[(i+2)%3]->P;
		    crossProduct(i) = ( (P-x0).cross(P-x1) ).norm();
		  }

		  Eigen::Matrix<double,3,1>::Index in;
		  double cP_min = crossProduct.minCoeff(&in);
		  assert(cP_min<0.01&&"error on getting the triangle edge on which an intersection point is on:: Triangle.h::findEdgeIndex (VectorDim P)");

		  ii = in;

		  return ii;
		}
		
		/*
		//============================================================================
		// function to return the triangle force vector resulted from infinite medium surface traction of a single dislocation segment
		//=============================================================================
		template < typename T>
		Eigen::Matrix<double,dim,3> getTriInfiniteForce_Seg (const T* const pt, const unsigned int iSeg) {
		  
		  Eigen::Matrix<double,dim,3> tractionInt=Eigen::Matrix<double,dim,3>::Zero();
		  		  		  
		  Eigen::Matrix<double,3,3> tractionMatrix;
		  for (unsigned int i = 0; i<3; i++){tractionMatrix.col(i) = this->eleNodes[i]->traction;}

		  integrate_Seg<T> (pt,iSeg ,tractionMatrix, tractionInt);
		  
		  return tractionInt;
		}
		
		//=================================================================================================
		// function to integrate the force coming from the infinite medium stress field of a given dislocation segment 
		// (over the customly definied Gauss points over the triangle)
		//=================================================================================================
		template < typename T>
		void integrate_Seg (const T* const pt, const unsigned int iSeg, Eigen::Matrix<double,3,3> tractionMatrix , Eigen::Matrix<double,dim,3>& tractionInt ) const {
		  
		  unsigned int iGlide = pt->link(cuttingSegments[iSeg].first,cuttingSegments[iSeg].second).second->pGlidePlane->sID;
		  
		  std::vector<VectorDim> abscissas = localQuadPnts.find(iGlide)->second;
		  double weight = localQuadPnts_w.find(iGlide)->second;
		  
		  for (unsigned int i=0; i<abscissas.size(); i++) {
		    tractionInt += dislocationStressKernel_Seg<T> (abscissas[i],pt,tractionMatrix,iSeg);
		  }
		  
		  tractionInt = weight*tractionInt;
		  
		}
		//============================================================================
		// infinite medium stress integration kernel
		//=============================================================================
		template<typename T> 
		Eigen::Matrix<double,dim,3> dislocationStressKernel_Seg(const VectorDim& R, const T* const pt, const Eigen::Matrix<double,3,3>& tractionMatrix,const unsigned int iSeg) const {
		  
			double J = 2.0*area() ;
			
			//---------- interpolate externally applied traction between nodes ---------------
			//Eigen::Matrix<double,1,3> shapeFunc = Shapefunc(Rstd);
			Eigen::Matrix<unsigned int,2,1> ixy =  find2DProjectionPlane();
			Eigen::Matrix<double,1,3> shapeFunc = getBarycentric (R , ixy(0) , ixy(1));       // shape functions w.r.t. the parent triangle
			
			VectorDim externalTraction = VectorDim::Zero();
			for (unsigned int i = 0; i<3; i++){externalTraction+=shapeFunc(i)*tractionMatrix.col(i);}
			
			//Eigen::Matrix<double,dim,dim> triStress = pts.pt->stress(R,this->sID);
			Eigen::Matrix<double,dim,dim> triStress = pt->link(cuttingSegments[iSeg].first,cuttingSegments[iSeg].second).second->stress_source(R);
			
			return (externalTraction-(triStress*outNormal))*shapeFunc*J ;
			//return (externalTraction-(pts.pt->stress(R,this->sID)*triN))*shapeFunc*J ;
		
		}
		*/
		
		//============================================================================
		// function to return the triangle force vector resulted from infinite medium surface traction of a single dislocation segment
		//=============================================================================
		template < typename T>
		Eigen::Matrix<double,dim,3> getTriInfiniteForce_gp (const T* const pt) const {
		  
		  PointerWrapper<T> pts(pt);
		  
		  Eigen::Matrix<double,dim,3> tractionInt=Eigen::Matrix<double,dim,3>::Zero();
		  		  		  
		  Eigen::Matrix<double,3,3> tractionMatrix;
		  for (unsigned int i = 0; i<3; i++){tractionMatrix.col(i) = this->eleNodes[i]->traction;}
		  
		  model::GlidePlaneObserver<typename T::LinkType> gpObsever;
//		  Eigen::Matrix<double,dim+1,1> gpKey;
		  
		  //std::cout << "Local points container size  " << localQuadPnts.size() << std::endl;
		  
		  //------------- loop over all glide planes, if localQuadPnts is found for it so integrate over them, 
		  //------------- otherwise integrate normally over the standard gauss points  -----------
		  for (typename model::GlidePlaneObserver<typename T::LinkType>::const_iterator gpIter=gpObsever.begin(); gpIter!=gpObsever.end(); ++gpIter){
		    
            const Eigen::Matrix<double,dim+1,1> gpKey((Eigen::Matrix<double,dim+1,1>()<<gpIter->second->planeNormal.normalized() , gpIter->second->height).finished());

		    //gpKey << gpIter->second->planeNormal.normalized() , gpIter->second->height;
		    
		    if (localQuadPnts.find(gpKey) != localQuadPnts.end()) integrate_gp<T> (pt ,tractionMatrix, tractionInt, gpKey);
		    
		    else model::Quadrature<2,3>::integrate(this,tractionInt,&Triangle::dislocationStressKernel_gp<PointerWrapper<T> > , pts, tractionMatrix , gpKey);
		    
		  }

		  return tractionInt;
		}
		
		//=================================================================================================
		// function to integrate the force coming from the infinite medium stress field of dislocation segments included on glide plane "gpKey" 
		// (over the customly definied Gauss points over the triangle)
		//=================================================================================================
		template < typename T>
		void integrate_gp (const T* const pt, Eigen::Matrix<double,3,3> tractionMatrix , Eigen::Matrix<double,dim,3>& tractionInt, const  Eigen::Matrix<double,dim+1,1> GlidePlaneKey ) const {
		  
		  //typename std::map<Eigen::Matrix<double,dim+1,1>, std::vector<VectorDim> , model::CompareVectorsByComponent<double,dim+1,float> >::const_iterator itt = localQuadPnts.find(GlidePlaneKey);
		  typename localQuadraturePointsContainerType::const_iterator itt = localQuadPnts.find(GlidePlaneKey);
		  
		  
		  //std::vector<VectorDim> abscissas = (*itt).second;
		  
		  vectorDimVectorType abscissas = (*itt).second;
		  
		  double weight = 0.5e00 / double (abscissas.size());
		  
		  double J = 2.0*area() ;
		  Eigen::Matrix<unsigned int,2,1> ixy =  find2DProjectionPlane();
		 
		  //Eigen::Matrix<double,1,3>::Index iBary;
		  //double baryMin;

		  for (unsigned int i=0; i<abscissas.size(); i++) {
		    
		    Eigen::Matrix<double,1,3> shapeFunc = getBarycentric (abscissas[i] , ixy(0) , ixy(1));       // shape functions at field point R
		    
		    //baryMin = shapeFunc.minCoeff(&iBary);
		    //assert ((baryMin>=0 && baryMin<=1.0) && "Integration point is outside triangle " );
		    
		    //std::cout << shapeFunc << std::endl;
		    
		    //---------- interpolate externally applied traction between nodes ---------------
		    VectorDim externalTraction = VectorDim::Zero();
		    for (unsigned int j = 0; j<3; j++) {externalTraction+=shapeFunc(j)*tractionMatrix.col(j);}
		    
		    tractionInt += (externalTraction-(pt->stressFromGlidePlane(GlidePlaneKey,abscissas[i]) *outNormal))*shapeFunc*J;
		    
		  }
		  
		  tractionInt = weight*tractionInt;
		  
		}
				
		//============================================================================
		// infinite medium stress integration kernel (used when integrating the stress field of dislocations belong to glide plane "GlidePlaneKey" over standard Gauss points)
		//=============================================================================
		template<typename T> 
		Eigen::Matrix<double,dim,3> dislocationStressKernel_gp(const Eigen::Matrix<double,2,1>& Rstd, const T& pts, const Eigen::Matrix<double,3,3>& tractionMatrix,
								    const Eigen::Matrix<double,dim+1,1>& GlidePlaneKey) const {
			double J ;
			VectorDim R=mapStdPoint<2>(Rstd,J);
			Eigen::Matrix<double,1,3> shapeFunc = Shapefunc(Rstd);
			
			//---------- interpolate externally applied traction between nodes ---------------
			VectorDim externalTraction = VectorDim::Zero();
			for (unsigned int i = 0; i<3; i++){externalTraction+=shapeFunc(i)*tractionMatrix.col(i);}
						
			return (externalTraction -( pts.pt->stressFromGlidePlane(GlidePlaneKey,R) * outNormal) )*shapeFunc*J ;
			//return (externalTraction -( pts.pt->stress(R) * outNormal) )*shapeFunc*J ;

            
		}
		/*
		//============================================================================
		// infinite medium stress integration kernel (used when integrating the stress field of dislocations belong to glide plane "GlidePlaneKey" over localQaudPnts)
		//=============================================================================
		template<typename T> 
		Eigen::Matrix<double,dim,3> dislocationStressKernel_gp(const VectorDim& R, const T* const pt, const Eigen::Matrix<double,3,3>& tractionMatrix,
								       const Eigen::Matrix<double,dim+1,1>& GlidePlaneKey) const {
		  
			double J = 2.0*area() ;
			
			Eigen::Matrix<unsigned int,2,1> ixy =  find2DProjectionPlane();
			Eigen::Matrix<double,1,3> shapeFunc = getBarycentric (R , ixy(0) , ixy(1));       // shape functions at field point R
			
			//---------- interpolate externally applied traction between nodes ---------------
			VectorDim externalTraction = VectorDim::Zero();
			for (unsigned int i = 0; i<3; i++){externalTraction+=shapeFunc(i)*tractionMatrix.col(i);}
						
			return (externalTraction-(pt->stressFromGlidePlane(GlidePlaneKey,R) *outNormal))*shapeFunc*J ;		
		}
		*/
		
		
		/*
		//============================================================================
		// function to return the sub-triangle force vector resulted from infinite medium surface traction
		//=============================================================================
		template <short unsigned int qOrder, typename T>
		Eigen::Matrix<double,dim,3> subTri_getTriInfiniteForce (const T* const pt , Eigen::Matrix<double,dim,dim> subTriPoints, unsigned int iSeg) {
		  
		  Eigen::Matrix<double,dim,3> tractionInt=Eigen::Matrix<double,dim,3>::Zero();
		  
		  PointerWrapper<T> pts(pt);
		  
		  VectorDim triN = triNormal();
		  
		  Eigen::Matrix<double,3,3> tractionMatrix;
		  for (unsigned int i = 0; i<3; i++){tractionMatrix.col(i) = this->eleNodes[i]->traction;}

		  model::Quadrature<2,qOrder>::integrate(this,tractionInt,&Triangle::subTri_dislocationStressKernel<PointerWrapper<T> > , pts, triN, tractionMatrix, subTriPoints , iSeg );
		  
		  return tractionInt;
		}
		
		
		//============================================================================
		// infinite medium stress integration kernel on subtriangle
		//=============================================================================
		template<typename T> 
		Eigen::Matrix<double,dim,3> subTri_dislocationStressKernel(const Eigen::Matrix<double,2,1>& Rstd, const T& pts, const VectorDim& triN, 
									   const Eigen::Matrix<double,3,3>& tractionMatrix,const Eigen::Matrix<double,dim,dim>& subTriPoints,
									   const unsigned int& iSeg) const {
			double J ;
			VectorDim R=subTri_mapStdPoint<2>(subTriPoints,Rstd,J);
			//std::cout<< R.transpose() << std::endl;
			//VectorDim R=mapStdPoint<2>(Rstd,J);
			
			//---------- interpolate externally applied traction between nodes ---------------
			VectorDim externalTraction = VectorDim::Zero();
			
			//Eigen::Matrix<double,1,3> shapeFunc = Shapefunc(Rstd);
			
			Eigen::Matrix<unsigned int,2,1> ixy =  find2DProjectionPlane();
			Eigen::Matrix<double,1,3> shapeFunc = getBarycentric (R , ixy(0) , ixy(1));       // shape functions w.r.t. the parent triangle
			
			for (unsigned int i = 0; i<3; i++){externalTraction+=shapeFunc(i)*tractionMatrix.col(i);}
			
			//Eigen::Matrix<double,dim,dim> triStress = pts.pt->stress(R,this->sID);
			Eigen::Matrix<double,dim,dim> triStress = pts.pt->link(cuttingSegments[iSeg].first,cuttingSegments[iSeg].second).second->stress_source(R);
			
			//std::cout<< this->sID << " : "<< triStress.col(0).transpose() << " " <<triStress.col(1).transpose() << " " <<triStress.col(2).transpose() <<std::endl;    
			//if (this->sID==428||this->sID==471||this->sID==1895||this->sID==1899||this->sID==1898||this->sID==442||this->sID==805) {
// 			if (this->sID==26||this->sID==457||this->sID==927||this->sID==456||this->sID==935 
// 			    || this->sID==428||this->sID==471||this->sID==1895||this->sID==1899||this->sID==1898||this->sID==442||this->sID==805) {
// 			  std::cout<< this->sID << " : "<< triStress.col(0).transpose() << " " <<triStress.col(1).transpose() << " " <<triStress.col(2).transpose() <<std::endl;  
// 			}
			return (externalTraction-(triStress*triN))*shapeFunc*J ;
			//return (externalTraction-(pts.pt->stress(R,this->sID)*triN))*shapeFunc*J ;
		
		}
		
		//==============================================================================
		// function to map the Gauss point to the actual domain in a subtriangle
		//==============================================================================
		template<short unsigned int OtherDim>
		VectorDim subTri_mapStdPoint(const Eigen::Matrix<double,dim,dim> subTriPoints, const Eigen::Matrix<double,2,1>& Rstd, double& J) const {
			
		  Eigen::Matrix<double,dim,OtherDim> A;
		  //VectorDim b = this->eleNodes[0]->P;
		  VectorDim b = subTriPoints.col(0);
		  
		  Eigen::Matrix<double,dim,1> temp;
		  
		    for (int i = 0; i<OtherDim; i++){
				//A.col(i) = this->eleNodes[i+1]->P - this->eleNodes[0]->P;
				A.col(i) = subTriPoints.col(i+1) - subTriPoints.col(0);
		    }
		
		    temp = A.col(0).cross(A.col(1));
		    J= temp.norm();    // this is = area of actual triangle * 2 , but the sum of integration weights for triangle is actually in this code (1/2)               
		    
		    assert(J>0.0);
		    return  A * Rstd +b ; 
		}
		*/
		
		
		/*
		//==========================================================================================
		// function to generate a custom set of quadrature points oriented with the glide plane - triangle intersection line
		//==================================================================================================
		
		void makeLocalQuadPoints (VectorDim P0, VectorDim P1,  const unsigned int gpID) {
		  
		  double dl = 20.0e00;                      // the separation distance between the quadrature points rows
		  
		  std::vector<VectorDim> QuadPointsSet;      // vector accumulates the generated points
		  std::vector<VectorDim> rowPointsSet;
		  
		  projPlaneIndx= find2DProjectionPlane();
		  
		  VectorDim cp = 0.5*(P0 + P1);  
		  VectorDim uv = (P1-P0).normalized();
		  VectorDim nd = outNormal.cross(uv);           // vector perpendicular to the intersection line, on the triangle plane
		  
		  bool done;
		  unsigned int iRow;
		  //--------- populate first in the +ve side of the nd ---------
		  done = false;
		  iRow = 0;
		  while (!done) {
		    iRow ++;
		    //rowPointsSet  = makeOneRow(dl,cp,uv,nd,iRow);
		    rowPointsSet  = makeOneRow_onLine(dl,cp,uv,nd,iRow);
		    
		    if(rowPointsSet.size()>0) QuadPointsSet.insert( QuadPointsSet.begin(), rowPointsSet.begin() , rowPointsSet.end() );
		    else done = true;
		  }
		  //--------- populate first in the -ve side of the nd ---------
		  done = false;
		  iRow = 0;
		  while (!done) {
		    iRow ++;
		    //rowPointsSet  = makeOneRow(dl,cp,uv,-nd,iRow);
		    rowPointsSet  = makeOneRow_onLine(dl,cp,uv,-nd,iRow);
		    
		    if(rowPointsSet.size()>0) QuadPointsSet.insert( QuadPointsSet.begin(), rowPointsSet.begin() , rowPointsSet.end() );
		    else done = true;
		  }
		  
		  assert(QuadPointsSet.size() > 0 && " failed to generate quadrature points around the glide plane intersection line on mesh triangles " );
		  
		  double w = 0.5e00 / double (QuadPointsSet.size());      // the sum of all weights should be 0.5
		  
		  //for (unsigned int i=0; i<QuadPointsSet.size(); i++) std::cout<< QuadPointsSet[i].transpose()<<std::endl; 
		  
		  localQuadPnts.insert  ( std::make_pair( gpID , QuadPointsSet ) );
		  localQuadPnts_w.insert( std::make_pair( gpID , w ) );
		  
		}
		*/
		//==================================================================================
		// function to add one row to the custom made quadrature points for each triangle
		//===================================================================================
		
		std::vector<VectorDim> makeOneRow(const double dl, const VectorDim cp, const VectorDim uv, const VectorDim nd, unsigned int iRow ) {
		  std::vector<VectorDim> rowPointsSet;
		  
		  VectorDim P0, P ;
		  Eigen::Matrix<double,3,1>::Index iBary;
		  Eigen::Matrix<double,3,1> Bary; 
		  double baryMin;
		  
		  bool endOfSide;
		  
		  //------- add a point at the center ------------
		  P0 = cp + (iRow*dl*nd);
		  Bary = getBarycentric(P0,projPlaneIndx(0),projPlaneIndx(1));
		  baryMin = Bary.minCoeff(&iBary);
		  //std::cout << Bary.transpose() << std::endl;
		  if (baryMin>0.0) rowPointsSet.push_back(P0);

		  //--------- move in +ve side of uv -----------
		  endOfSide = false;
		  P = P0;
		  
		  while (!endOfSide) {
		    P = P + (dl*uv);
		    Bary = getBarycentric(P,projPlaneIndx(0),projPlaneIndx(1));
		    //std::cout << Bary.transpose() << std::endl;
		    endOfSide = Bary.minCoeff(&iBary) <= 0.0;
		    
		    if (!endOfSide) rowPointsSet.push_back(P);
		  }
		  
		  //--------- move in -ve side of uv -----------
		  endOfSide = false;
		  P = P0;
		  
		  while (!endOfSide) {
		    P = P - (dl*uv);
		    Bary = getBarycentric(P,projPlaneIndx(0),projPlaneIndx(1));
		    //std::cout << Bary.transpose() << std::endl;
		    endOfSide = Bary.minCoeff(&iBary) <= 0.0;
		    
		    if (!endOfSide) rowPointsSet.push_back(P);
		  }
		  
		  return rowPointsSet;
		}
		
		
		//==================================================================================
		// function to add one row to the custom made quadrature points for each triangle
		//===================================================================================
		
		std::vector<VectorDim> makeOneRow_onLine(const double dl, const VectorDim cp, const VectorDim uv, const VectorDim nd, unsigned int iRow ) {
		  std::vector<VectorDim> rowPointsSet;
		  
		  VectorDim c0, p0 , p1 , x0 , x1 , P , v ;
		  unsigned int nPnts;
		  double dx , l;
		  
		  //------- get end points of the new line ------------
		  c0 = cp + (iRow*dl*nd);
		  //std::cout<< c0.transpose() << std::endl;
		  p0 = c0 - (1.0e4*uv);
		  p1 = c0 + (1.0e4*uv);
		  
		  //------------------ implement intersection  with the 3 triangle edges -----------
		  if (intersectLineWithEdges(p0,p1,x0,x1)) {
		    l = (x1-x0).norm();
		    v = (x1-x0).normalized();
		    nPnts = int( l / dl) + 1;              // number of points on the line
		    dx = 0.5e00 * ( l - ((nPnts-1)*dl) );
		    
		    for (unsigned int i=0; i<nPnts; i++){
		      P = x0 + (dx+(i*dl))*v;
		      rowPointsSet.push_back(P);
		    }
		  }    
		  return rowPointsSet;
		}
		
		//===========================================================================================
		// function to intersect a line with the 3 triangle edges, and update the 2 intersection points
		//===========================================================================================
		
		bool intersectLineWithEdges (VectorDim p0, VectorDim p1, VectorDim& x0, VectorDim& x1) {
		  //VectorDim x0 , x1 ;
		  bool intersectionFound = false;
		  
		  Eigen::Matrix<double,2,2> mat;
		  Eigen::Matrix<double,2,1> u,b;
		  
		  std::vector<VectorDim> intersectionPoints;
		  VectorDim temp;
		  
		  Eigen::Matrix<unsigned int,2,1> ppindx= find2DProjectionPlane();
		  unsigned int ix = ppindx(0);
		  unsigned int iy = ppindx(1);
		  
		  for(unsigned int i=0; i<3; i++){               // loop over the 3 triangle edges
		    x0 = this->eleNodes[(i+1)%3]->P;
		    x1 = this->eleNodes[(i+2)%3]->P;
		    
		    mat(0,0) = p1(ix) - p0(ix);        mat(0,1) = x0(ix) - x1(ix);    
		    mat(1,0) = p1(iy) - p0(iy);        mat(1,1) = x0(iy) - x1(iy);
		    
		    b(0) = x0(ix) - p0(ix);
		    b(1) = x0(iy) - p0(iy);
		    
		    u = mat.inverse()*b;
		    
		    if ( u(1) >= 0.0e00  && u(1) <=1.0e00) {       // intersection point is within the triangle line segments
		      temp = x0 + u(1)*(x1-x0);
		      intersectionPoints.push_back(temp);
		    }
		  }
		  
		  if(intersectionPoints.size()==3) {      // when one of the intersection points is a vertex
		    int repeatedPoint = -1; 
		    for (unsigned int i1=0; i1<intersectionPoints.size();i1++){
		      for (unsigned int i2=i1+1; i2<intersectionPoints.size(); i2++){
			if ((intersectionPoints[i1]-intersectionPoints[i2]).norm() < 1.0e-8 ) repeatedPoint = i1; 
		      }
		    }
		    
		    if (repeatedPoint != -1) intersectionPoints.erase(intersectionPoints.begin()+repeatedPoint);
		  }
		  
		  if(intersectionPoints.size()==2) {     // found 2 intersection points
		    x0 = intersectionPoints[0];
		    x1 = intersectionPoints[1];
		    intersectionFound = true;
		  }
		  else intersectionFound = false;
		  
		  return intersectionFound;
		}
		
		//==============================================================================
		Eigen::Matrix<double,3,1>  getBarycentric (VectorDim P , unsigned int ix  , unsigned int iy  ) const
		{
			Eigen::Matrix<double,3,1> bary;
			double V0;
			
			V0 = getVol(this->eleNodes[0]->P,this->eleNodes[1]->P,this->eleNodes[2]->P ,ix, iy);
			
			bary(0) = getVol(                   P,this->eleNodes[1]->P,this->eleNodes[2]->P,  ix,  iy)/V0;
			bary(1) = getVol(this->eleNodes[0]->P,                   P,this->eleNodes[2]->P,  ix,  iy)/V0;
			bary(2) = getVol(this->eleNodes[0]->P,this->eleNodes[1]->P,                   P,  ix,  iy)/V0;
			
			return bary;
		}
		
		//======================================================================
		// return the volume of the tetradedron, given its points coordinates in order
		//======================================================================
		double getVol(VectorDim a, VectorDim b, VectorDim c , int ix, int iy) const
		{

			MatrixDim temp;
			
			temp << a(ix), b(ix), c(ix),
			        a(iy), b(iy), c(iy),
				1   , 1   , 1; 
				
			return temp.determinant()/2.0;
		}

		
		/////////////////////////////////////////////////// DEBUGGING FUNCTIONS ////////////////////////////////////////////

		
		
		//============================================================================
		// function to return the triangle infinite medium surface traction vector resulted from infinite medium dislocations field
		//=============================================================================
		template <short unsigned int qOrder, bool deformed = false , typename T>
		Eigen::Matrix<double,dim,1> getTriInfiniteTraction (const T* const pt) const {
		  
		  Eigen::Matrix<double,dim,1> tractionInt=Eigen::Matrix<double,dim,1>::Zero();
		  
		  PointerWrapper<T> pts(pt);
		  VectorDim triN;     
		  if (deformed) triN = triNormalDeformed();     else triN = outNormal;
		  
		  model::GlidePlaneObserver<typename T::LinkType> gpObsever;
		  Eigen::Matrix<double,dim+1,1> gpKey;
		  
		  for (typename model::GlidePlaneObserver<typename T::LinkType>::const_iterator gpIter=gpObsever.begin(); gpIter!=gpObsever.end(); ++gpIter){
		    
		    gpKey << gpIter->second->planeNormal.normalized() , gpIter->second->height;
		    
		    if (localQuadPnts.find(gpKey) != localQuadPnts.end()) integrate_gp<T> (pt, tractionInt, gpKey, triN );
			    
		    else model::Quadrature<2,qOrder>::integrate(this,tractionInt,&Triangle::dislocationStressKernel<PointerWrapper<T> > , pts , gpKey, triN );
		    
		  }

		  return tractionInt;
		}

		//=================================================================================================
		// function to integrate the force coming from the infinite medium stress field of a given dislocation segment 
		// (over the customly definied Gauss points over the triangle)
		//=================================================================================================
		template < typename T>
		void integrate_gp (const T* const pt , Eigen::Matrix<double,dim,1>& tractionInt, const  Eigen::Matrix<double,dim+1,1> GlidePlaneKey, const VectorDim triN ) const {

		  typename localQuadraturePointsContainerType::const_iterator itt = localQuadPnts.find(GlidePlaneKey);
		  vectorDimVectorType abscissas = (*itt).second;
		  
		  double weight = 1.0e00 / double (abscissas.size());
		  
		  Eigen::Matrix<double,3,3> temp = Eigen::Matrix<double,3,3>::Zero();
		  
		  for (unsigned int i=0; i<abscissas.size(); i++) {
		    temp+= pt->stressFromGlidePlane(GlidePlaneKey,abscissas[i]);
		  }		  
		  tractionInt = weight*(temp*triN);  
		}

		//============================================================================
		// infinite medium surface traction integration kernel
		//=============================================================================
		template<typename T> 
		Eigen::Matrix<double,dim,1> dislocationStressKernel(const Eigen::Matrix<double,2,1>& Rstd, const T& pts,
								    const Eigen::Matrix<double,dim+1,1>& GlidePlaneKey, const VectorDim& triN) const {
			double J ;
			VectorDim R=mapStdPoint<2>(Rstd,J);
			return (pts.pt->stressFromGlidePlane(GlidePlaneKey,R)*triN*2.0);  // multiply by 2 because the sum of the integration weights = 0.5;
			//return (pts.pt->stress(R)*triN)*J ;
		}
	
		//============================================================================
		// function to return the triangle infinite medium force vector resulted from infinite medium dislocations field
		//=============================================================================
		template <short unsigned int qOrder, bool deformed = false , typename T>
		Eigen::Matrix<double,dim,1> forceInfinite (const T* const pt) const {
		  
		  double triArea;
		  
		  if (deformed) triArea = deformed_area();
		  else triArea = area();
		  
		  return triArea*getTriInfiniteTraction<qOrder,deformed,T>(pt);
		}
		
		
	};
	
//------------ redefine the static variable [required] ------------
	double Triangle::N[3][3] ;
	double Triangle::edgeTol = 1.0e-02;
	
}  //  namespace bvpfe
#endif