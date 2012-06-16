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
#include <Eigen/Core>
#include <model/Utilities/StaticID.h>
#include <model/Dislocations/GlidePlanes/GlidePlaneObserver.h>
#include <model/Dislocations/DislocationSharedObjects.h>

namespace model {
	
	/*************************************************************/
	/*************************************************************/
	template <short unsigned int dim, typename SegmentType>
	class GlidePlane : public  StaticID<GlidePlane<dim,SegmentType> >,
	/*              */ private std::set<const SegmentType*>,
	/*              */ private GlidePlaneObserver<dim,SegmentType>{
		
public:

		typedef GlidePlane<dim,SegmentType> GlidePlaneType;
		typedef typename SegmentType::TempMaterialType MaterialType;
		typedef Eigen::Matrix<double,dim,dim> MatrixDimD;
		typedef Eigen::Matrix<double,dim,1> VectorDimD;
		typedef Eigen::Matrix<double,dim+1,1> VectorDimPlusOneD;
		typedef std::pair<VectorDimD,VectorDimD> segmentMeshCollisionPair;
		//typedef std::vector<segmentMeshCollisionPair> segmentMeshCollisionPairContainerType;
//		typedef std::map<unsigned int , segmentMeshCollisionPair> segmentMeshCollisionPairContainerType;
		typedef std::map<unsigned int, segmentMeshCollisionPair, std::less<unsigned int>,
		/*            */ Eigen::aligned_allocator<std::pair<const unsigned int, segmentMeshCollisionPair> > > segmentMeshCollisionPairContainerType;

		typedef  std::pair<unsigned int, segmentMeshCollisionPair> planeTraingleIntersection;
//		typedef  std::vector<planeTraingleIntersection> planeMeshIntersectionType;
		typedef  std::vector<planeTraingleIntersection, Eigen::aligned_allocator<planeTraingleIntersection> > planeMeshIntersectionType;

		private:
		static DislocationSharedObjects<dim,MaterialType> shared;
		
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
				shared.domain.get_planeMeshIntersection(planeNormal*height,planeNormal,segmentMeshCollisionPairContainer,this->sID);
				//----------- intersect the mesh with parallel planes to this glide plane -------------
				intersectMeshWithParallelPlanes (planeNormal,height);
			}
		}
		
		/* Desctructor  **********************************************/
		~GlidePlane(){
			std::cout<<"Deleting GlidePlane "<<this->sID<<std::endl;
			assert(this->glidePlaneMap.erase((VectorDimPlusOneD()<< planeNormal, height).finished())==1 && "CANNOT ERASE GLIDE PLANE  FROM STATIC glidePlaneMap.");
			assert(this->empty() && "DELETING NON EMPTY GLIDE PLANE.");
		}
		
		/* addToGLidePlane  ******************************************/
		void addToGLidePlane(SegmentType* const pS){
			assert(this->insert(pS).second && "COULD NOT INSERT SEGMENT POINTER IN SEGMENT CONTAINER.");
		}
		
		/* removeFromGlidePlane  *************************************/
		void removeFromGlidePlane(SegmentType* const pS){
			assert(this->erase(pS)==1 && "COULD NOT ERASE SEGMENT POINTER IN SEGMENT CONTAINER.");
		}
		
		/* begin  **************************************************/
		typename std::set<const SegmentType*>::const_iterator begin() const {
			return std::set<const SegmentType*>::begin();
		}
		
		/* end  **************************************************/
		typename std::set<const SegmentType*>::const_iterator end() const {
			return std::set<const SegmentType*>::end();
		}
		
		/* stress ************************************************/
		MatrixDimD stress(const VectorDimD& Rfield) const {
			MatrixDimD temp(MatrixDimD::Zero());
			for (typename std::set<const SegmentType*>::const_iterator sIter = this->begin(); sIter != this->end() ;++sIter){
		    temp+=(*sIter)->stress_source(Rfield);
		  	}
			return temp;
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
		
		//===========================================================================
		// function to intersect the mesh with planes parallel to the current slip plane
		//============================================================================
		void intersectMeshWithParallelPlanes (const VectorDimD planeNormal,const double height) {
		  
		  unsigned int nPlanes = 10;            // number of parallel planes on each side of the glise plane
		  double separation = 6.0e00;         // the separation distance between each plane
		  
		  //segmentMeshCollisionPairContainerType temp;
		  
		  for (unsigned int i=1; i<= nPlanes; i++ ) {
		    
		    addParallelPlane(planeNormal*(height+(i*separation)) , planeNormal , separation );
		    addParallelPlane(planeNormal*(height-(i*separation)) , planeNormal , separation );
	
		  }
		  
		}
		
		//=================================================================================
		// function to intersect a plane parallel to this glide plane with the mesh surface
		//=================================================================================
		
		void addParallelPlane(VectorDimD x0, VectorDimD planeNormal, double separation) {
		  
		  segmentMeshCollisionPairContainerType temp;
		  		    
		  shared.domain.get_planeMeshIntersection(x0,planeNormal,temp,this->sID);
		    
		  planeMeshIntersectionType contourPointsVector = sortIntersectionLines(temp);
		  
		  generateIntegrationPoints (contourPointsVector, separation);
		}
		
		//===========================================================================
		// function to sort the plane-mesh intersection line to be in a sequence
		//===========================================================================
		
		planeMeshIntersectionType sortIntersectionLines(segmentMeshCollisionPairContainerType temp_in) {
		  
		  double tol = 1.0e-8;
		  
		  segmentMeshCollisionPairContainerType temp = temp_in;
		  
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
		    
		    if (!found) {
		      
		      typename segmentMeshCollisionPairContainerType::iterator itt;

		      for (itt = temp.begin(); itt != temp.end() ;++itt){
			
			if ((P-(*itt).second.first).norm()<tol) {
			  P = (*itt).second.second;
			  linesVector.push_back(std::make_pair((*itt).first,std::make_pair((*itt).second.first , (*itt).second.second)));
			  iTriNext = (*itt).first;
			  found = true;
			  break;
			}
			
			else if ((P-(*itt).second.second).norm()<tol) {
			  P = (*itt).second.first;
			  linesVector.push_back(std::make_pair((*itt).first,std::make_pair((*itt).second.second , (*itt).second.first)));
			  iTriNext = (*itt).first;
			  found = true;
			  break;
			}
			
		      }
		      
		      if (found) temp.erase(iTriNext);
		    
		    }
		  
		    assert(found && "unable to find next segment on which the glide plane intersects the mesh boundary");
		    
		    remaining = temp.size();
		}
		
		assert (linesVector.size() == temp_in.size() && "missing line segments during the glide plane intersection lines sorting process " );
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
		      if(shared.domain.triContainer[contourPointsVector[i].first]->localQuadPnts.find(this->sID) != shared.domain.triContainer[contourPointsVector[i].first]->localQuadPnts.end()) {
			shared.domain.triContainer[contourPointsVector[i].first]->localQuadPnts.find(this->sID)->second.push_back(P);
		      }
		      else {
			//std::vector<VectorDimD> QuadPointsSet;
			//QuadPointsSet.push_back(P);
			std::vector<VectorDimD> QuadPointsSet (1,P);
			shared.domain.triContainer[contourPointsVector[i].first]->localQuadPnts.insert( std::make_pair( this->sID , QuadPointsSet ) );
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

