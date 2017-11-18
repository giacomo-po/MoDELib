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

#include <deque>
#include <chrono>
#include <memory>
#include <map>
#include <set>
#include <assert.h>
#include <Eigen/Core>
#include <model/Utilities/StaticID.h>
#include <model/DislocationDynamics/GlidePlanes/GlidePlaneObserver.h>
#include <model/LatticeMath/LatticePlane.h>
#include <model/MPI/MPIcout.h>
#include <model/Mesh/PlaneMeshIntersection.h>

#ifndef NDEBUG
#define VerboseGlidePlane(N,x) if(verboseGlidePlane>=N){model::cout<<x;}
#else
#define VerboseGlidePlane(N,x)
#endif

namespace model
{
    
    /**************************************************************************/
    /**************************************************************************/
    template <int dim>
    struct GlidePlane : public StaticID<GlidePlane<dim> >,
    /* base class    */ public LatticePlane,
    /* base class    */ private std::set<const std::shared_ptr<GlidePlane<dim>>*>
    {
        
        typedef GlidePlane<dim> GlidePlaneType;
        typedef GlidePlaneObserver<dim> GlidePlaneObserverType;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef typename GlidePlaneObserverType::GlidePlaneKeyType GlidePlaneKeyType;
        typedef typename PlaneMeshIntersection<dim>::PlaneMeshIntersectionContainerType PlaneMeshIntersectionContainerType;

        static int verboseGlidePlane;
        GlidePlaneObserverType* const glidePlaneObserver;
//        const Grain<dim>& grain; // REMOVED TO ALLOW CONSTRUCTION FROM A GRAINBOUNDARY
        const std::pair<size_t,size_t> grainIDs;
        const GlidePlaneKeyType glidePlaneKey;
        const PlaneMeshIntersectionContainerType meshIntersections;
        
        
        /**********************************************************************/
        GlidePlane(GlidePlaneObserverType* const gpo,
                   const SimplicialMesh<dim>& mesh,
                   const Grain<dim>& grain_in,
                   const VectorDim& P,
                   const VectorDim& N) :
        /* init */ LatticePlane(P,grain_in.reciprocalLatticeDirection(N)), // BETTER TO CONSTRUCT N WITH PRIMITIVE VECTORS ON THE PLANE
        /* init */ glidePlaneObserver(gpo),
//        /* init */ grain(grain_in),
        /* init */ grainIDs(grain_in.grainID,grain_in.grainID),
        /* init */ glidePlaneKey(GlidePlaneObserverType::getGlidePlaneKey(grain_in,P,N)),
        /* init */ meshIntersections(PlaneMeshIntersection<dim>(mesh,this->P,this->n.cartesian().normalized(),grain_in.grainID))
        {
            VerboseGlidePlane(1,"Creating GlidePlane "<<this->sID<<" ("<<glidePlaneKey.transpose()<<")"<<std::endl;);
            glidePlaneObserver->addGlidePlane(this);
            assert(fabs(this->unitNormal.norm()-1.0)<DBL_EPSILON && "GlidePlane has non-unit normal.");
        }
        
//        /**********************************************************************/
//        GlidePlane(GlidePlaneObserverType* const gpo,
//                   const SimplicialMesh<dim>& mesh,
//                   const GrainBoundary<dim>& gb_in,
//                   const VectorDim& P,
//                   const VectorDim& N) :
//        /* init */ LatticePlane(P,grain_in.reciprocalLatticeDirection(N)), // BETTER TO CONSTRUCT N WITH PRIMITIVE VECTORS ON THE PLANE
//        /* init */ glidePlaneObserver(gpo),
//        /* init */ grain(grain_in),
//        /* init */ glidePlaneKey(GlidePlaneObserverType::getGlidePlaneKey(grain_in,P,N)),
//        /* init */ meshIntersections(PlaneMeshIntersection<dim>(mesh,this->P,this->n.cartesian().normalized(),grain_in.grainID))
//        {
//            static_assert(false, "ALLOW CONSTRUCTOR FROM A GRAIN BOUNDARY AS OPPOSED TO A GRAIN. THIS WILL ALLOW GB PLANES FROM CSL AND DSCL");
//            VerboseGlidePlane(1,"Creating GlidePlane "<<this->sID<<" ("<<glidePlaneKey.transpose()<<")"<<std::endl;);
//            glidePlaneObserver->addGlidePlane(this);
//            assert(fabs(this->unitNormal.norm()-1.0)<DBL_EPSILON && "GlidePlane has non-unit normal.");
//        }
        
        /**********************************************************************/
        GlidePlane(const GlidePlane<dim>& other) = delete;
        
        /**********************************************************************/
        ~GlidePlane()
        {
            VerboseGlidePlane(1,"Destroying GlidePlane "<<this->sID<<" ("<<glidePlaneKey.transpose()<<")"<<std::endl;);
            glidePlaneObserver->removeGlidePlane(this);
        }
        
        /**********************************************************************/
        std::shared_ptr<GlidePlane<dim>> sharedPlane() const
        {/*!\returns a shared pointer to this
          */
            assert(this->size());
            assert((*this->begin())->get()==this);
            return **this->begin();
        }
        
//        /**********************************************************************/
//        template<typename LoopType>
//        void addLoop(const LoopType* const pL)
//        {/*!@\param[in] pS a row pointer to a DislocationSegment
//          * Adds pS to *this GLidePlane
//          */
//            const bool success=this->insert(&pL->_glidePlane).second;
//            assert( success && "COULD NOT INSERT LOOP POINTER IN GLIDE PLANE.");
//        }
//        
//        /**********************************************************************/
//        template<typename LoopType>
//        void removeLoop(const LoopType* const pL)
//        {/*!@\param[in] pS a row pointer to a DislocationSegment
//          * Removes pS from *this GLidePlane
//          */
//            const int success=this->erase(&pL->_glidePlane);
//            assert(success==1 && "COULD NOT ERASE LOOP POINTER FROM GLIDE PLANE.");
//        }
        
        /**********************************************************************/
        void addParentSharedPtr(const std::shared_ptr<GlidePlane<dim>>* const pL)
        {/*!@\param[in] pS a row pointer to a DislocationSegment
          * Adds pS to *this GLidePlane
          */
            const bool success=this->insert(pL).second;
            assert( success && "COULD NOT INSERT LOOP POINTER IN GLIDE PLANE.");
        }
        
        /**********************************************************************/
//        template<typename LoopType>
        void removeParentSharedPtr(const std::shared_ptr<GlidePlane<dim>>* const pL)
        {/*!@\param[in] pS a row pointer to a DislocationSegment
          * Removes pS from *this GLidePlane
          */
            const int success=this->erase(pL);
            assert(success==1 && "COULD NOT ERASE LOOP POINTER FROM GLIDE PLANE.");
        }
        
        /**********************************************************************/
        template <class T>
        friend T& operator << (T& os, const GlidePlaneType& gp)
        {
            size_t kk=0;
            for (const auto& x : gp.meshIntersections)
            {
                os<<gp.sID<< " "<<kk<<" "<< x.first->child(0).xID<< " "<< x.first->child(1).xID <<" "<<x.second.transpose()<<"\n";
                kk++;
            }
            
            return os;
        }
        
    };
    
    template <int dim>
    int GlidePlane<dim>::verboseGlidePlane=0;
    
}
#endif

