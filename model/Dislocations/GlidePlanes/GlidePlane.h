/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2011 by Mamdouh Mohamed <mamdouh.s.mohamed@gmail.com>, 
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_GLIDEPLANE_H
#define model_GLIDEPLANE_H
#include <math.h>
#include <set>
#include <assert.h>
#include <boost/shared_ptr.hpp>
#include <Eigen/Core>
#include<Eigen/StdVector>
#include <model/Utilities/StaticID.h>
#include <model/Dislocations/DislocationNetworkTraits.h>
#include <model/Dislocations/GlidePlanes/GlidePlaneObserver.h>
#include <model/Dislocations/DislocationSharedObjects.h>
#include <model/BVP/VirtualBoundarySlipSurface.h>


namespace model {
	
	/*************************************************************/
	/*************************************************************/
	template <typename SegmentType>
	class GlidePlane : public  StaticID<GlidePlane<SegmentType> >,
	/*              */ private std::set<const SegmentType*>,
	/*              */ private std::set<const VirtualBoundarySlipSurface<SegmentType>*>,
	/*              */ private GlidePlaneObserver<SegmentType>{
		
public:
	typedef GlidePlane<SegmentType> GlidePlaneType;
	typedef boost::shared_ptr<GlidePlaneType> GlidePlaneSharedPtrType;

        enum{dim=TypeTraits<SegmentType>::dim};
        
		typedef std::set<const SegmentType*> 				         SegmentContainerType;
		typedef std::set<const VirtualBoundarySlipSurface<SegmentType>*> BoundarySegmentContainerType;
			
		//typedef GlidePlane<SegmentType> GlidePlaneType;
		//typedef typename SegmentType::TempMaterialType MaterialType;
		typedef Eigen::Matrix<double,dim,dim> MatrixDimD;
		typedef Eigen::Matrix<double,dim,1> VectorDimD;
		typedef Eigen::Matrix<double,dim+1,1> VectorDimPlusOneD;
		typedef std::pair<VectorDimD,VectorDimD> segmentMeshCollisionPair;
		//typedef std::vector<segmentMeshCollisionPair> segmentMeshCollisionPairContainerType;
		//typedef std::map<unsigned int , segmentMeshCollisionPair> segmentMeshCollisionPairContainerType;
		typedef std::map<unsigned int, segmentMeshCollisionPair, std::less<unsigned int>,
		/*            */ Eigen::aligned_allocator<std::pair<const unsigned int, segmentMeshCollisionPair> > > segmentMeshCollisionPairContainerType;

		typedef  std::pair<unsigned int, segmentMeshCollisionPair> planeTraingleIntersection;
//		typedef  std::vector<planeTraingleIntersection> planeMeshIntersectionType;
		typedef  std::vector<planeTraingleIntersection, Eigen::aligned_allocator<planeTraingleIntersection> > planeMeshIntersectionType;

		private:
		static DislocationSharedObjects<SegmentType> shared;
		
	public:
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW

		
		//! Make this const and initialize at constructor time
		segmentMeshCollisionPairContainerType segmentMeshCollisionPairContainer;
		
		std::vector <planeMeshIntersectionType> planesMeshIntersectionContainer;         // stores the intersection lines of planes (parallel to this glide plane) with the mesh
		
		
		//! The non-unit vector defining the plane in space
		const VectorDimD planeNormal;
		const double height;
		
		/* Constructor with non-unit plane normal ******************************/
		GlidePlane(const VectorDimD& planeNormal_in, const double& height_in) : planeNormal(planeNormal_in),
		/* init list                                                         */ height(height_in){
//std::cout<<"planeNormal.norm()="<<planeNormal.norm()<<std::endl;
			assert(std::fabs(planeNormal.norm()-1.0)<=DBL_EPSILON && "GLIDE PLANE NORMAL IS NOT UNIT");	
//			assert(planeNormal.norm()==1.0 && "GLIDE PLANE NORMAL IS NOT UNIT");
			std::cout<<"Creating GlidePlane "<<this->sID<<" with unit plane normal: "<<planeNormal.transpose()<<" and height "<<height<<std::endl;
			assert(this->glidePlaneMap.insert(std::make_pair((VectorDimPlusOneD()<< planeNormal, height).finished(),this)).second && "CANNOT INSERT GLIDE PLANE  IN STATIC glidePlaneMap.");
	if (shared.boundary_type){			
				shared.domain.get_planeMeshIntersection(planeNormal*height,planeNormal,segmentMeshCollisionPairContainer);
				//----------- intersect the mesh with parallel planes to this glide plane -------------
				//intersectMeshWithParallelPlanes (planeNormal,height);
				intersectMeshWithParallelPlanes ();
			}
		}
		
		/* Desctructor  **********************************************/
		~GlidePlane(){
			std::cout<<"Deleting GlidePlane "<<this->sID<<std::endl;
			assert(this->glidePlaneMap.erase((VectorDimPlusOneD()<< planeNormal, height).finished())==1 && "CANNOT ERASE GLIDE PLANE  FROM STATIC glidePlaneMap.");
			assert(        SegmentContainerType::empty() && "DELETING NON-EMPTY GLIDE PLANE.");
			assert(BoundarySegmentContainerType::empty() && "DELETING NON-EMPTY GLIDE PLANE.");
		}

		/* getSharedPointer ******************************************/
		GlidePlaneSharedPtrType getSharedPointer() const{
			GlidePlaneSharedPtrType temp;
			if (!SegmentContainerType::empty()){
				temp=(*SegmentContainerType::begin())->pGlidePlane;
			}
			else{ // no segments in the container
				if (!BoundarySegmentContainerType::empty()){
					temp=(*BoundarySegmentContainerType::begin())->pGlidePlane;
				}
				else{
					assert(0 && "GLIDE PLANE CONTAINS NO SEGMENTS AND NO BOUNDARY SEGMENTS");
				}
			}
		return temp;
		}
		
		/* addToGLidePlane  ******************************************/
		void addToGLidePlane(SegmentType* const pS){
			assert(SegmentContainerType::insert(pS).second && "COULD NOT INSERT SEGMENT POINTER IN SEGMENT CONTAINER.");
		}
		
		/* removeFromGlidePlane  *************************************/
		void removeFromGlidePlane(SegmentType* const pS){
			assert(SegmentContainerType::erase(pS)==1 && "COULD NOT ERASE SEGMENT POINTER FROM SEGMENT CONTAINER.");
		}

		/* addToGLidePlane  ******************************************/
		void addToGLidePlane(VirtualBoundarySlipSurface<SegmentType>* const pS){
			assert(BoundarySegmentContainerType::insert(pS).second && "COULD NOT INSERT BOUNDARY SEGMENT POINTER IN SEGMENT CONTAINER.");
		}
		
		/* removeFromGlidePlane  *************************************/
		void removeFromGlidePlane(VirtualBoundarySlipSurface<SegmentType>* const pS){
			assert(BoundarySegmentContainerType::erase(pS)==1 && "COULD NOT ERASE BOUNDARY SEGMENT POINTER FROM SEGMENT CONTAINER.");
		}
		
		/* begin  **************************************************/
		typename std::set<const SegmentType*>::const_iterator begin() const {
			return SegmentContainerType::begin();
		}
		
		/* end  **************************************************/
		typename std::set<const SegmentType*>::const_iterator end() const {
			return SegmentContainerType::end();
		}
		
		
		/* friend T& operator << **********************************************/
		template <class T>
		friend T& operator << (T& os, const GlidePlaneType& gp){
		  unsigned int ii = 0 ;
		  typename segmentMeshCollisionPairContainerType::const_iterator itt;
		  for (itt = gp.segmentMeshCollisionPairContainer.begin(); itt != gp.segmentMeshCollisionPairContainer.end() ;++itt){
		    os << gp.sID<< " "<< ii <<" "
		    << (*itt).second.first.transpose()<<" "
		    << (*itt).second.second.transpose()<<"\n";
		    ii++;
		  }
			/*for (unsigned int k=0; k<gp.segmentMeshCollisionPairContainer.size();++k){
				os << gp.sID<< " "<< k <<" "
				   << gp.segmentMeshCollisionPairContainer[k].first.transpose()<<" "
				   << gp.segmentMeshCollisionPairContainer[k].second.transpose();
			}*/
			return os;
		}
		
		/* stress from all segments on this glide plane ************************************************/
		MatrixDimD stress(const VectorDimD& Rfield) const {
		  MatrixDimD temp(MatrixDimD::Zero());
		  for (typename std::set<const SegmentType*>::const_iterator sIter = SegmentContainerType::begin(); sIter != SegmentContainerType::end() ;++sIter){
		    temp+=(*sIter)->stress_source(Rfield);
		  }
		  return temp;
		}
		
		//===========================================================================
		// function to intersect the mesh with planes parallel to the current slip plane
		//============================================================================
		//void intersectMeshWithParallelPlanes (const VectorDimD planeNormal,const double height) {
		void intersectMeshWithParallelPlanes () {
		  
		  unsigned int nPlanes = 3;            // number of parallel planes on each side of the glise plane
		  double separation = 50.0e00;         // the separation distance between each plane
		  
		  bool isLast = false;
		  
		  //segmentMeshCollisionPairContainerType temp;
		  
		  for (unsigned int i=1; i<= nPlanes; i++ ) {
		    
		    if (i==nPlanes) isLast = true;
		    
		    //addParallelPlane(planeNormal*(height+(i*separation)) , planeNormal , separation );
		    //addParallelPlane(planeNormal*(height-(i*separation)) , planeNormal , separation );
		    
		    //addParallelPlane(height , planeNormal, separation, i ,  1 , isLast );
		    //addParallelPlane(height , planeNormal, separation, i , -1 , isLast );
		    addParallelPlane(separation,  i  , isLast );
		    addParallelPlane(separation, -i ,  isLast );
		  }
		  
		}
		
		//=================================================================================
		// function to intersect a plane parallel to this glide plane with the mesh surface
		//=================================================================================
		
		//void addParallelPlane(VectorDimD x0, VectorDimD planeNormal, double separation) {
		//  void addParallelPlane(const double height, const VectorDimD planeNormal, double separation, unsigned int iP, int sign, const bool isLast) {
		void addParallelPlane(const double separation, const int iP, const bool isLast) {
	
		    VectorDimD x0 = (height+(iP*separation))*planeNormal; 
   
		    segmentMeshCollisionPairContainerType collContainer, temp ;
		    
		    shared.domain.get_planeMeshIntersection(x0,planeNormal,collContainer);
		    
		    temp = collContainer;
	   
		    if (temp.size()>0) {
		      planeMeshIntersectionType contourPointsVector = sortIntersectionLines(temp);
		      generateIntegrationPoints (contourPointsVector, separation);
		      
		      //------------ if this is the last generated parallel plane, add additional points so that triangles will not have just few points
		      // ----------- on one just on side of it -----------------------------------------------------------------------------------------
		      //if (isLast) completeTrianglesPopulation(collContainer,height,planeNormal,separation,iP,sign);
		      if (isLast) completeTrianglesPopulation(collContainer,separation,iP);
		    }
		    
		  }
		  
		//===================================================================================================
		// function to add more points to triangles that has been cut with planes parallel to the glide plane
		//===================================================================================================
		
		void completeTrianglesPopulation(const segmentMeshCollisionPairContainerType collContainer,const double separation,const int iP){
		  
		  bool allDone = false;
		  int ii = std::abs(iP);
		  int sign = iP/ii;
		  
		  //std::cout << iP << " " << ii << " " << sign << std::endl;
		  
		  VectorDimD x0;
		  
		  while (!allDone) {
		    ii++;
		    
		    //std::cout << height << " " << sign << " " << ii << " " << separation << " " << sign*ii << " " << height+(sign*ii*separation) << std::endl;
		    x0 = planeNormal*(height+(sign*ii*separation));
		    
		    segmentMeshCollisionPairContainerType temp ;
		    shared.domain.get_planeMeshIntersection(x0,planeNormal,temp);
		    
		    //std::cout<< temp.size() << " ";
		    
		    if (temp.size()>0) {
		      bool foundTriangles = false;
		      for (typename segmentMeshCollisionPairContainerType::iterator itt=temp.begin(); itt!= temp.end(); itt++ ) {
			if (collContainer.find(itt->first)!= collContainer.end()) {
			  
			  addPointToTriangle(itt->first, itt->second, separation);
			  foundTriangles = true;
			}
		      }
		      
		      if (!foundTriangles) allDone = true;
		      
		    }
		    
		    else allDone = true;
		    
		    //std::cout<< allDone<< std::endl;
		    
		  }
		}
		
		//============================================================================
		// function to add more integration points to a specific triangle
		//===========================================================================
		
		void addPointToTriangle (const unsigned int triID, const segmentMeshCollisionPair pointsPair, const double separation) {
		  
		  double l = (pointsPair.second - pointsPair.first).norm();
		  int nPnts = int(l/ separation) + 1;
		  
		  VectorDimD dir = (pointsPair.second - pointsPair.first).normalized();
		  
		  double dx = 0.5e00*(l - ((nPnts-1)*separation));
		  //std::cout<<l << " " << separation << " " << nPnts << " " << dx << " : ";
		  //std::cout << pointsPair.first.transpose() << std::endl;
		  //std::cout << pointsPair.second.transpose() << std::endl;
		  
		  for (int i=0; i<nPnts; i++) {
		    
		    //std::cout << dx << " ";
		    //VectorDimD P = pointsPair.first + (dx+(i*separation))*dir;
		    VectorDimD P = pointsPair.first + dx*dir;
		    
		    if(shared.domain.triContainer[triID]->localQuadPnts.find( (VectorDimPlusOneD()<<planeNormal.normalized(),height).finished() ) != 
		      shared.domain.triContainer[triID]->localQuadPnts.end() ) {
		      //std::cout << P.transpose()<<std::endl;
		      shared.domain.triContainer[triID]->localQuadPnts.find( (VectorDimPlusOneD()<<planeNormal.normalized(),height).finished() )->second.push_back(P);
		      }
		      
		      dx+=separation;
		  }
		  //std::cout << std::endl;
		  
		}
		
		//===========================================================================
		// function to sort the plane-mesh intersection line to be in a sequence
		//===========================================================================
		
		planeMeshIntersectionType sortIntersectionLines(segmentMeshCollisionPairContainerType & temp) {
		  
		  double tol = 1.0e-7;
		  
		  unsigned int input_size = temp.size();
		  		  		  
		  planeMeshIntersectionType linesVector;
		  VectorDimD P;
		  
		  unsigned int iTriBegin = temp.begin()->first;              // the index of the triangle on  which search begin
		  unsigned int iTriNext  = iTriBegin;                      
		  
		  unsigned int triID;
		  bool found;
		  
		  //------------- insert the first line -------------------
		  linesVector.push_back(std::make_pair(iTriBegin,std::make_pair(temp.begin()->second.first , temp.begin()->second.second)));
		  P = temp.begin()->second.second;
		  temp.erase(iTriBegin);         // already sorted, so remove it from the list
		  
		  unsigned int remaining = temp.size();
		  		  
		  while(remaining != 0) {
		    		    
		    found = false;
		    
		    for (unsigned int in=0; in<3; in++) {
		      
		      triID = shared.domain.triContainer[iTriNext]->neighbor[in]->sID;
		      
		      if (temp.find(triID) != temp.end()) {
			
			if ((P-(*temp.find(triID)).second.first).norm()<tol) {
			  P = temp.find(triID)->second.second;
			  linesVector.push_back(std::make_pair(triID,std::make_pair(temp.find(triID)->second.first , temp.find(triID)->second.second)));
			  iTriNext = triID;
			  found = true;
			  temp.erase(iTriNext);
			  break;
			}
			
			else if ((P-(*temp.find(triID)).second.second).norm()<tol) {
			  P = temp.find(triID)->second.first;
			  linesVector.push_back(std::make_pair(triID,std::make_pair(temp.find(triID)->second.second , temp.find(triID)->second.first)));
			  iTriNext = triID;
			  found = true;
			  temp.erase(iTriNext);
			  break;
			}
		      }
		    }
		    
		    double minDis = 1.0e08;
		    
		    if (!found) {
		      
		      typename segmentMeshCollisionPairContainerType::iterator itt , itt_min;
		      
		      for (itt = temp.begin(); itt != temp.end() ;++itt){
			if      ((P-(*itt).second.first).norm() < minDis)       { minDis = (P-(*itt).second.first).norm();    itt_min=itt;}
			if      ((P-(*itt).second.second).norm()< minDis)       { minDis = (P-(*itt).second.second).norm();   itt_min=itt;} 			
		      }
		      
		      
		      if ((P-(*itt_min).second.first).norm() < (P-(*itt_min).second.second).norm()) {
			P = (*itt_min).second.second;
			linesVector.push_back(std::make_pair((*itt_min).first,std::make_pair((*itt_min).second.first , (*itt_min).second.second)));
			iTriNext = (*itt_min).first;
			found = true;
		      }
		      else {
			P = (*itt_min).second.first;
			linesVector.push_back(std::make_pair((*itt_min).first,std::make_pair((*itt_min).second.second , (*itt_min).second.first)));
			iTriNext = (*itt_min).first;
			found = true;
		      }
		      
		      if (found) temp.erase(iTriNext);
		      
		      if (minDis > 1.0e-2) std::cout<< "Warning: distance between two consecutive lines is " << minDis << std::endl;
		      
		    }
		  
//  		  if (!found) {
//  		    std::cout <<"==================== Glide plane-mesh intersection lines ========================" << std::endl;
// 		    
// 		    for (unsigned int i=0; i<linesVector.size(); i++) std::cout<< linesVector[i].second.first.transpose() << "    " << linesVector[i].second.second.transpose() << std::endl;
// 		    std::cout << std::endl;
//  		    for (typename segmentMeshCollisionPairContainerType::iterator itt = temp.begin(); itt != temp.end() ;++itt) 
// 		           std::cout << (*itt).second.first.transpose() << "  " << (*itt).second.second.transpose()<< std::endl; ;
//  		    
//  		    std::cout<< "===================================================================================== " << std::endl;
//  		  }
		  		  
		  assert(found && "unable to find next segment on which the glide plane intersects the mesh boundary");
		  
		  remaining = temp.size();
		}
		
		assert (linesVector.size() == input_size && "missing line segments during the glide plane intersection lines sorting process " );
		planesMeshIntersectionContainer.push_back(linesVector);
		
		return linesVector;
		
		}
		
		//=======================================================================================================
		// function to generate integration points over the line segments, and distribute them over the triangles
		//=======================================================================================================
		void generateIntegrationPoints (planeMeshIntersectionType contourPointsVector, double separation){
		  
		  double length = 0.0e00;
		  for (unsigned int i=0; i<contourPointsVector.size(); i++) {
		    length += (contourPointsVector[i].second.second - contourPointsVector[i].second.first).norm(); 
		  }
		  
		  unsigned int nPnts = int (length / separation);
		  double dx = length / double(nPnts);                 // the actual separation distance to be used
		  
		  VectorDimD P;
		  
		  double dl = 0.5e00*dx;           // just initial value for the first point 
		  
		  for (unsigned int i=0; i<contourPointsVector.size(); i++) {
		    
		    double l = (contourPointsVector[i].second.second - contourPointsVector[i].second.first).norm();
		    VectorDimD uv = (contourPointsVector[i].second.second - contourPointsVector[i].second.first).normalized();
		    
		    while (dl <= l) {
		      P = contourPointsVector[i].second.first + dl*uv;
		      //std::cout << P.transpose()<<std::endl;
		      //-------- check if there is a container for this glide plane already exists or not -------------
		      if(shared.domain.triContainer[contourPointsVector[i].first]->localQuadPnts.find( (VectorDimPlusOneD()<<planeNormal.normalized(),height).finished() ) != 
			 shared.domain.triContainer[contourPointsVector[i].first]->localQuadPnts.end() ) {
			shared.domain.triContainer[contourPointsVector[i].first]->localQuadPnts.find( (VectorDimPlusOneD()<<planeNormal.normalized(),height).finished() )->second.push_back(P);
		      }
		      else {
			//std::vector<VectorDimD> QuadPointsSet (1,P);
			std::vector< VectorDimD , Eigen::aligned_allocator<VectorDimD> > QuadPointsSet (1,P);
			shared.domain.triContainer[contourPointsVector[i].first]->localQuadPnts.insert( std::make_pair( (VectorDimPlusOneD()<<planeNormal.normalized(),height).finished() , QuadPointsSet ) );
		      }
		      
		      dl+= dx;
		    }
		    
		    dl -= l; 
		  } 
		  
		}
		
		
		
		
	};
	
	/*************************************************************/
	/*************************************************************/
} // namespace model
#endif

//std::set<const SegmentType*> segmentContainer;
//assert(segmentContainer.empty() && "DELETING NON EMPTY GLIDE PLANE.");
//assert(segmentContainer.insert(pS).second && "COULD NOT INSERT SEGMENT POINTER IN SEGMENT CONTAINER.");
//			assert(segmentContainer.erase(pS)==1 && "COULD NOT ERASE SEGMENT POINTER IN SEGMENT CONTAINER.");

