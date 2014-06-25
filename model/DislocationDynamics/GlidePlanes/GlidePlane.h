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
#include <deque>
#include <set>
#include <assert.h>

#include <Eigen/Core>
#include <Eigen/StdVector>

#include <model/Utilities/StaticID.h>
#include <model/DislocationDynamics/DislocationNetworkTraits.h>
#include <model/DislocationDynamics/GlidePlanes/GlidePlaneObserver.h>
#include <model/DislocationDynamics/DislocationSharedObjects.h>
#include <model/DislocationDynamics/VirtualBoundarySlipSurface.h>
#include <model/MPI/MPIcout.h>
#include <model/Mesh/SimplexObserver.h>


namespace model {
	
	/*************************************************************/
	/*************************************************************/
	template <typename SegmentType>
	class GlidePlane :
    /* base class   */ public  StaticID<GlidePlane<SegmentType> >,
	/* base class   */ private std::set<const SegmentType*>,
	/* base class   */ private std::set<const VirtualBoundarySlipSurface<SegmentType>*>,
	/* base class   */ private GlidePlaneObserver<SegmentType>
    {
		
    public:
        typedef GlidePlane<SegmentType> GlidePlaneType;
        typedef typename GlidePlaneObserver<SegmentType>::GlidePlaneSharedPtrType GlidePlaneSharedPtrType;
        
        enum{dim=TypeTraits<SegmentType>::dim};
        
		typedef std::set<const SegmentType*> SegmentContainerType;
		typedef std::set<const VirtualBoundarySlipSurface<SegmentType>*> BoundarySegmentContainerType;
        
		typedef Eigen::Matrix<double,dim,dim> MatrixDimD;
		typedef Eigen::Matrix<double,dim,1> VectorDimD;
		typedef Eigen::Matrix<double,dim+1,1> VectorDimPlusOneD;
		typedef std::pair<VectorDimD,VectorDimD> segmentMeshCollisionPair;
		typedef std::deque<segmentMeshCollisionPair> SegmentMeshCollisionPairContainerType;
		typedef std::pair<unsigned int, segmentMeshCollisionPair> planeTraingleIntersection;
		typedef std::vector<planeTraingleIntersection, Eigen::aligned_allocator<planeTraingleIntersection> > planeMeshIntersectionType;
        
    private:
		static DislocationSharedObjects<SegmentType> shared;
		
        
        /**********************************************************************/
        SegmentMeshCollisionPairContainerType getPlaneMeshIntersection(const VectorDimD& x0, const VectorDimD& n) const
        {/*!@param[in] x0 a point on *this GlidePlane
          *!@param[n]  n  the unit plane normal to *this GlidePlane
          *\returns a container of segments representing the intersection of *this GlidePLane with the SimplexMesh
          */
            SegmentMeshCollisionPairContainerType temp;
            if (shared.boundary_type)
            {
                for(typename SimplexObserver<dim,dim-1>::const_iterator fIter =SimplexObserver<dim,dim-1>::simplexBegin();
                    /*                                               */ fIter!=SimplexObserver<dim,dim-1>::simplexEnd();++fIter)
                {
                    if(fIter->second->isBoundarySimplex())
                    {
                        std::deque<VectorDimD> intersectionPoints;
                        for(unsigned int e=0;e<dim;++e)
                        {// this part of the code is specific to dim=3;
                            const VectorDimD& v0(fIter->second->child(e).child(0).P0);
                            const VectorDimD& v1(fIter->second->child(e).child(1).P0);
                            const double u ((x0-v0).dot(n) / (v1-v0).dot(n));
                            if(u>=0 && u<=1.0)
                            {
                                const VectorDimD P(v0 + u*(v1-v0));
                                bool isDifferent=true;
                                for (unsigned int k=0;k<intersectionPoints.size();++k)
                                {
                                    isDifferent*= (P-intersectionPoints[k]).squaredNorm()>FLT_EPSILON;
                                }
                                
                                if (isDifferent)
                                {
                                    intersectionPoints.emplace_back(P);
                                }
                            }
                        }
                        assert(intersectionPoints.size()<3);
                        if (intersectionPoints.size()==2)
                        {
                            temp.emplace_back(intersectionPoints[0],intersectionPoints[1]);
                        }
                    }
                }
            }
            return temp;
        }
        
        
        
	public:
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        
        
        //! The unit vector normal to *this GlidePlane
		const VectorDimD planeNormal;
		
        //! The height from the origin along the planeNormal
        const double height;
        
		
		//! Make this const and initialize at constructor time
		const SegmentMeshCollisionPairContainerType segmentMeshCollisionPairContainer;

		/**********************************************************************/
		GlidePlane(const VectorDimD& planeNormal_in, const double& height_in) :
        /* init list */ planeNormal(planeNormal_in),
		/* init list */ height(height_in),
        /* init list */ segmentMeshCollisionPairContainer(getPlaneMeshIntersection(planeNormal*height,planeNormal))
        {
			assert(std::fabs(planeNormal.norm()-1.0)<=DBL_EPSILON && "GLIDE PLANE NORMAL IS NOT UNIT");
			assert(this->glidePlaneMap.insert(std::make_pair((VectorDimPlusOneD()<< planeNormal, height).finished(),this)).second && "CANNOT INSERT GLIDE PLANE  IN STATIC glidePlaneMap.");
		}
		
        /**********************************************************************/
		~GlidePlane()
        {
			model::cout<<"Deleting GlidePlane "<<this->sID<<std::endl;
			assert(this->glidePlaneMap.erase((VectorDimPlusOneD()<< planeNormal, height).finished())==1 && "CANNOT ERASE GLIDE PLANE  FROM STATIC glidePlaneMap.");
			assert(        SegmentContainerType::empty() && "DELETING NON-EMPTY GLIDE PLANE.");
			assert(BoundarySegmentContainerType::empty() && "DELETING NON-EMPTY GLIDE PLANE.");
		}
        
        /**********************************************************************/
		GlidePlaneSharedPtrType getSharedPointer() const
        {
			GlidePlaneSharedPtrType temp;
			if (!SegmentContainerType::empty())
            {
				temp=(*SegmentContainerType::begin())->pGlidePlane;
			}
			else
            { // no segments in the container
				if (!BoundarySegmentContainerType::empty())
                {
					temp=(*BoundarySegmentContainerType::begin())->pGlidePlane;
				}
				else
                {
					assert(0 && "GLIDE PLANE CONTAINS NO SEGMENTS AND NO BOUNDARY SEGMENTS");
				}
			}
            return temp;
		}
		
        /**********************************************************************/
		void addToGLidePlane(SegmentType* const pS)
        {
			assert(SegmentContainerType::insert(pS).second && "COULD NOT INSERT SEGMENT POINTER IN SEGMENT CONTAINER.");
		}
		
        /**********************************************************************/
		void removeFromGlidePlane(SegmentType* const pS)
        {
			assert(SegmentContainerType::erase(pS)==1 && "COULD NOT ERASE SEGMENT POINTER FROM SEGMENT CONTAINER.");
		}
        
        /**********************************************************************/
		void addToGLidePlane(VirtualBoundarySlipSurface<SegmentType>* const pS)
        {
			assert(BoundarySegmentContainerType::insert(pS).second && "COULD NOT INSERT BOUNDARY SEGMENT POINTER IN SEGMENT CONTAINER.");
		}
		
        /**********************************************************************/
		void removeFromGlidePlane(VirtualBoundarySlipSurface<SegmentType>* const pS)
        {
			assert(BoundarySegmentContainerType::erase(pS)==1 && "COULD NOT ERASE BOUNDARY SEGMENT POINTER FROM SEGMENT CONTAINER.");
		}
		
        /**********************************************************************/
		typename std::set<const SegmentType*>::const_iterator begin() const
        {
			return SegmentContainerType::begin();
		}
		
        /**********************************************************************/
		typename std::set<const SegmentType*>::const_iterator end() const
        {
			return SegmentContainerType::end();
		}
		
        /**********************************************************************/
        MatrixDimD stress(const VectorDimD& Rfield) const
        {/*  stress from all segments on this glide plane
          */
            MatrixDimD temp(MatrixDimD::Zero());
            for (typename std::set<const SegmentType*>::const_iterator sIter = SegmentContainerType::begin(); sIter != SegmentContainerType::end() ;++sIter){
                temp+=(*sIter)->stress_source(Rfield);
            }
            return temp;
        }
        
        /* friend T& operator << **********************************************/
		template <class T>
		friend T& operator << (T& os, const GlidePlaneType& gp)
        {
            unsigned int ii = 0 ;
            typename SegmentMeshCollisionPairContainerType::const_iterator itt;
            for (itt = gp.segmentMeshCollisionPairContainer.begin(); itt != gp.segmentMeshCollisionPairContainer.end() ;++itt){
                os << gp.sID<< " "<< ii <<" "
                << itt->first.transpose()<<" "
                << itt->second.transpose()<<"\n";
                ii++;
            }
            return os;
        }
        
    };
    
    /*************************************************************/
    /*************************************************************/
} // namespace model
#endif

