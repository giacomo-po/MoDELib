/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpF@ucla.edu>.
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
#include <StaticID.h>
#include <SimplexTraits.h>
#include <SimplicialMesh.h>
#include <GlidePlaneFactory.h>
#include <LatticePlane.h>
#include <MPIcout.h>
//#include <PlaneMeshIntersection.h>
#include <MeshPlane.h>

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
    struct GlidePlane :
    /* base class    */ public LatticePlane,
    /* base class    */ public MeshPlane<dim>
//    /* base class    */ private std::set<const std::shared_ptr<GlidePlane<dim>>*>
    {

        typedef GlidePlane<dim> GlidePlaneType;
        typedef GlidePlaneFactory<dim> GlidePlaneFactoryType;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        //typedef std::array<long int,dim+3> GlidePlaneKeyType;
        typedef typename TypeTraits<GlidePlaneType>::KeyType KeyType;
        typedef KeyType GlidePlaneKeyType;


//        /**********************************************************************/
//        static GlidePlaneKeyType getGlidePlaneKey(const int& grainID,
//                                                  const LatticePlane& lp
//        //                                                const VectorDimD& N
//        )
//        {/*!\param[in] grain the grain on which the GlidePlane is defined
//          * \param[in] P a point on the plane
//          * \param[in] N the normal to the plane
//          * \returns the key which uniquely identifies the plane.
//          * The type of the key is a tuple with entries (grainID,r,h), where r
//          * is the ReciprocalLatticeDirection corresponding to N, and h=P.dot(r)
//          * is an integer indicating the "heigth" of the plane from the origin,
//          * in integer multiples of the interplanar distance d=1/|r|.
//          */
//            GlidePlaneKeyType temp;
//            temp[0]=grainID;
//            temp[1]=grainID;
//            for(int d=0;d<dim;++d)
//            {
//                temp[2+d]=lp.n(d); // reciprocal lattice direction components
//            }
//            temp[2+dim]=lp.h; // integer height (units of reciprocal lattice direction)
//            return temp;
//        }

        static int verboseGlidePlane;
        const  GlidePlaneFactoryType& glidePlaneFactory;
        const Grain<dim>& grain;
        const GlidePlaneKeyType key;

//        /**********************************************************************/
//        GlidePlane(const GlidePlaneFactoryType* const gpF,
//                   const VectorDim& P,
//                   const VectorDim& N) :
//        /* init */ LatticePlane(P,_grain.reciprocalLatticeDirection(N)) // BETTER TO CONSTRUCT N WITH PRIMITIVE VECTORS ON THE PLANE
//        /* init */,MeshPlane<dim>(mesh,_grain.grainID,this->planeOrigin(),this->n.cartesian())
//        /* init */,GlidePlaneFactory(gpF)
//        /* init */,grain(_grain)
//        /* init */,glidePlaneKey(getGlidePlaneKey(grain.grainID,*this))
//        {
//            VerboseGlidePlane(1,"Creating GlidePlane "<<this->sID<<std::endl;);
//            GlidePlaneFactory->addGlidePlane(this);
//        }

        /**********************************************************************/
        GlidePlane(const GlidePlaneFactoryType& gpF,
                   const GlidePlaneKeyType& key_in) :
        /* init */ LatticePlane(key_in.h,ReciprocalLatticeDirection<dim>(key_in.r,gpF.poly.grain(key_in.grainID))) // BETTER TO CONSTRUCT N WITH PRIMITIVE VECTORS ON THE PLANE
        /* init */,MeshPlane<dim>(gpF.poly.mesh,key_in.grainID,this->planeOrigin(),this->n.cartesian())
        /* init */,glidePlaneFactory(gpF)
        /* init */,grain(gpF.poly.grain(key_in.grainID))
        /* init */,key(key_in)
        {
            VerboseGlidePlane(1,"Creating GlidePlane "<<this->sID<<std::endl;);
//            GlidePlaneFactory->addGlidePlane(this);
        }
        
        /**********************************************************************/
        GlidePlane(const GlidePlane<dim>& other) = delete;
        
        /**********************************************************************/
        ~GlidePlane()
        {
//            VerboseGlidePlane(1,"Destroying GlidePlane "<<this->sID<<" ("<<glidePlaneKey.transpose()<<")"<<std::endl;);
            VerboseGlidePlane(1,"Destroying GlidePlane "<<this->sID<<std::endl;);
//            GlidePlaneFactory->removeGlidePlane(this);
        }

//        /**********************************************************************/
//        std::shared_ptr<GlidePlane<dim>> sharedPlane() const
//        {/*!\returns a shared pointer to this
//          */
//            assert(this->size());
//            assert((*this->begin())->get()==this);
//            return **this->begin();
//        }

//        /**********************************************************************/
//        void addParentSharedPtr(const std::shared_ptr<GlidePlane<dim>>* const pL)
////                                const bool& isGlissile,
////                                const size_t& loopID)
//        {/*!@\param[in] pS a row pointer to a DislocationSegment
//          * Adds pS to *this GLidePlane
//          */
//            const bool success=this->insert(pL).second;
//            assert( success && "COULD NOT INSERT LOOP POINTER IN GLIDE PLANE.");
//
////            if(isGlissile)
////            {
////                glissileLoopIDs.insert(loopID);
////            }
////            else
////            {
////                sessileLoopIDs.insert(loopID);
////            }
//        }

//        /**********************************************************************/
////        template<typename LoopType>
//        void removeParentSharedPtr(const std::shared_ptr<GlidePlane<dim>>* const pL)
////                                   const bool& isGlissile,
////                                   const size_t& loopID)
//        {/*!@\param[in] pS a row pointer to a DislocationSegment
//          * Removes pS from *this GLidePlane
//          */
//            const int success=this->erase(pL);
//            assert(success==1 && "COULD NOT ERASE LOOP POINTER FROM GLIDE PLANE.");
//
////            if(isGlissile)
////            {
////                glissileLoopIDs.erase(loopID);
////            }
////            else
////            {
////                sessileLoopIDs.erase(loopID);
////            }
//        }



    };

    template <int dim>
    int GlidePlane<dim>::verboseGlidePlane=0;

}
#endif



//        /**********************************************************************/
//        template <class T>
//        friend T& operator << (T& os, const GlidePlaneType& gp)
//        {
//            size_t kk=0;
//            for (const auto& x : gp.meshIntersections)
//            {
////                os<<gp.sID<< " "<<kk<<" "<< x.first->child(0).xID<< " "<< x.first->child(1).xID <<" "<<x.second.transpose()<<"\n";
//                os<<gp.sID<< " "<<x.face->sID<<" "<< x.P0.transpose()<< " "<< x.P1.transpose()<<"\n";
//                kk++;
//            }
//
//            return os;
//        }


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
