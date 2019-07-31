/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2011 by Mamdouh Mohamed <mamdouh.s.mohamed@gmail.com>,
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_GlidePlaneFactory_H_
#define model_GlidePlaneFactory_H_

#include <algorithm>
#include <map>
#include <tuple>
#include <memory> // std::shared_ptr (c++11)

#include <Eigen/Dense>
#include <Eigen/Cholesky>

#include <SimplicialMesh.h>
#include <CompareVectorsByComponent.h>
#include <LatticeMath.h>
#include <ReciprocalLatticeDirection.h>
#include <PlanePlaneIntersection.h>
#include <MeshPlane.h>
#include <CRTP.h>
#include <TypeTraits.h>
#include <Polycrystal.h>

namespace model
{
    
    
    template<typename Derived>
    struct KeyConstructableSharedPtrFactory : public CRTP<Derived>
    /*                                     */,public std::map<typename TypeTraits<Derived>::KeyType,
    /*                                                      */ const std::shared_ptr<typename TypeTraits<Derived>::ValueType>,
    /*                                                      */ typename TypeTraits<Derived>::CompareType>
    {
        
        typedef typename TypeTraits<Derived>::KeyType   KeyType;
        typedef typename TypeTraits<Derived>::ValueType ValueType;
        typedef std::shared_ptr<ValueType> SharedPtrType;
        typedef std::map<KeyType,const SharedPtrType> MapType;
        typedef typename MapType::size_type SizeType;
        
        
//        KeyConstructableSharedPtrFactory& factory()
//        {
//            return *this;
//        }
//
//        const KeyConstructableSharedPtrFactory& factory() const
//        {
//            return *this;
//        }
        
//        /**************************************************************************/
//        template <typename ...ArgTypes>
//        SharedPtrType operator[](const ArgTypes&... args)
//        {
//            return operator[](KeyType(args...));
//        }
        
        /**************************************************************************/
        SharedPtrType get(const KeyType& key)
        {
            const auto iter(this->find(key));
            if(iter!=this->end())
            {// key found
                return iter->second;
            }
            else
            {
                return this->emplace(key,new ValueType(this->derived(),key)).first->second;
            }

        }

        
//        /**************************************************************************/
//        SharedPtrType get(const KeyType& key)
//        {
//            return MapType::emplace(std::piecewise_construct,
//                                   std::forward_as_tuple(key),
//                                   std::forward_as_tuple(this->derived(),key)
//                                   ).first->second;
//        }
        
//        /**************************************************************************/
//        SizeType erase(const SharedPtrType& shared) FINISH THIS PART
//        {
//            const auto iter(MapType::find(shared->key()));
//            if(iter!=MapType::end())
//            {
//                if(iter->second->use_count()<=2)
//                {// current map
//
//                }
//            }
//            else
//            {
//                return 0;
//            }
//        }
        
    };
    
    
    template <int dim>
    struct GlidePlane; // class predeclaration
    
    template<int dim>
    struct GlidePlaneFactory;
    
    template<int dim>
    struct GlidePlaneKey : public std::array<long int,dim+3>
    {
        
        typedef Eigen::Matrix<double,dim,1> VectorDimD;
        typedef Eigen::Matrix<long int,dim,1> VectorDimI;
        
        const long int& grainID;
        const Eigen::Map<const VectorDimI> r;
        const long int& h;
//
        GlidePlaneKey(const std::array<long int,dim+3>& array) :
        /* init */ std::array<long int,dim+3>(array)
        /* init */,grainID(this->operator[](0))
        /* init */,r(&this->data()[2])
        /* init */,h(this->operator[](5))
        {}
        
        GlidePlaneKey(const GlidePlaneKey& other) :
        /* init */ std::array<long int,dim+3>(other)
        /* init */,grainID(this->operator[](0))
        /* init */,r(&this->data()[2])
        /* init */,h(this->operator[](5))
        {}

        /**********************************************************************/
        GlidePlaneKey(const int& grainID_in,
                      const LatticePlane& lp) :
        /* init */ grainID(this->operator[](0))
        /* init */,r(&this->data()[2])
        /* init */,h(this->operator[](5))
        {/*!\param[in] grain the grain on which the GlidePlane is defined
          * \param[in] P a point on the plane
          * \param[in] N the normal to the plane
          * \returns the key which uniquely identifies the plane.
          * The type of the key is a tuple with entries (grainID,r,h), where r
          * is the ReciprocalLatticeDirection corresponding to N, and h=P.dot(r)
          * is an integer indicating the "heigth" of the plane from the origin,
          * in integer multiples of the interplanar distance d=1/|r|.
          */
            this->operator[](0)=grainID_in;
            this->operator[](1)=grainID_in;
            for(int d=0;d<dim;++d)
            {
                this->operator[](2+d)=lp.n(d); // reciprocal lattice direction components
            }
            this->operator[](2+dim)=lp.h; // integer height (units of reciprocal lattice direction)
        }
        
        /**********************************************************************/
        GlidePlaneKey(const int& grainID,
                      const VectorDimD& P,
                      const ReciprocalLatticeDirection<dim>& r) :
        /* delegate */ GlidePlaneKey(grainID,LatticePlane(P,r))
        {
        }
                
        /**********************************************************************/
        GlidePlaneKey(const int& grainID,
                      const LatticeVector<dim>& L,
                      const ReciprocalLatticeDirection<dim>& r) :
        /* delegate */ GlidePlaneKey(grainID,LatticePlane(L,r))
        {
            
        }
        
    };
    
    template<int dim>
    struct GlidePlaneFactoryTraits
    {
        typedef GlidePlane<dim> ValueType;
        typedef GlidePlaneKey<dim> KeyType;
        typedef std::less<std::array<long int,dim+3>> CompareType;

    };
    
    template<int dim>
    struct TypeTraits<GlidePlaneFactory<dim>> : public GlidePlaneFactoryTraits<dim>
    {

    };
    
    template<int dim>
    struct TypeTraits<GlidePlane<dim>> : public GlidePlaneFactoryTraits<dim>
    {
        
    };
    
    /**************************************************************************/
    
//    struct GlidePlaneFactory : private std::map<std::array<long int,dim+3>,const GlidePlane<dim>* const>
//    /*                       */,private std::map<std::pair<size_t,size_t>,PlanePlaneIntersection<dim>>
    template<int dim>
    struct GlidePlaneFactory : public KeyConstructableSharedPtrFactory<GlidePlaneFactory<dim>>
    /*                      */,private std::map<std::pair<size_t,size_t>,PlanePlaneIntersection<dim>>
    {

        typedef GlidePlaneFactory<dim> GlidePlaneFactoryType;
        typedef GlidePlane<dim> GlidePlaneType;
        typedef MeshPlane<dim> MeshPlaneType;
        typedef Eigen::Matrix<double,dim,1> VectorDimD;
        typedef typename GlidePlaneType::GlidePlaneKeyType GlidePlaneKeyType;
        typedef KeyConstructableSharedPtrFactory<GlidePlaneFactory<dim>> GlidePlaneMapType;
        typedef std::shared_ptr<GlidePlaneType> GlidePlaneSharedPtrType;
        typedef PlanePlaneIntersection<dim> PlanePlaneIntersectionType;
        typedef std::map<std::pair<size_t,size_t>,PlanePlaneIntersectionType> MeshPlaneIntersectionContainerType;
        
    public:
        
        const Polycrystal<dim>& poly;
        
        GlidePlaneFactory(const Polycrystal<dim>& poly_in) :
        /* init */ poly(poly_in)
        {
            
        }
        
        /**********************************************************************/
        MeshPlaneIntersectionContainerType& glidePlaneIntersections()
        {
            return *this;
        }
        
        /**********************************************************************/
        const MeshPlaneIntersectionContainerType& glidePlaneIntersections() const
        {
            return *this;
        }
        
        /**********************************************************************/
        const PlanePlaneIntersectionType& glidePlaneIntersection(const MeshPlaneType* const p1,
                                                                 const MeshPlaneType* const p2)
        {/*@param[in] p1 first plane
          *@param[in] p2 second plane
          *\returns the infinite line of interseciton between the two planes,
          * if already computed. Otherwise it computes the interseciton, stores it,
          * and returns it.
          */
            const auto key=std::make_pair(std::max(p1->sID,p2->sID),std::min(p1->sID,p2->sID));
            const auto iter=glidePlaneIntersections().find(key);
            if(iter==glidePlaneIntersections().end())
            {
                const bool success=glidePlaneIntersections().emplace(std::piecewise_construct,
                                                                     std::make_tuple(key),
                                                                     std::make_tuple(p1->P,p1->unitNormal,p2->P,p2->unitNormal)
                                                                     ).second;
                assert(success);
            }
            return glidePlaneIntersections().at(key);
        }
        
        /**********************************************************************/
        GlidePlaneMapType& glidePlanes()
        {
            return *this;
        }
        
        /**********************************************************************/
        const GlidePlaneMapType& glidePlanes() const
        {
            return *this;
        }
        
//        /**********************************************************************/
//        std::shared_ptr<GlidePlaneType> sharedGlidePlane(const SimplicialMesh<dim>& mesh,
//                                                         const Grain<dim>& grain,
//                                                         const VectorDimD& P,
//                                                         const VectorDimD& N)
//        {
//            const LatticePlane lp(P,grain.reciprocalLatticeDirection(N));
//            const GlidePlaneKeyType key=GlidePlaneType::getGlidePlaneKey(grain.grainID,lp);
//            const auto planeIter=glidePlanes().find(key);
//            return (planeIter!=glidePlanes().end())? planeIter->second->sharedPlane() : std::shared_ptr<GlidePlaneType>(new GlidePlaneType(this,mesh,grain,P,N));
//            //            return (planeIter!=glidePlanes().end())? planeIter->second->loops().begin()->second->_glidePlane :
//            //            /*                            */ std::shared_ptr<GlidePlaneType>(new GlidePlaneType(this,mesh,grain,P,N));
//        }
        
//        /**********************************************************************/
//        void addGlidePlane(const GlidePlaneType* const pL)
//        {/*!@\param[in] pS a row pointer to a DislocationSegment
//          * Adds pS to *this GLidePlane
//          */
//            const bool success=glidePlanes().emplace(pL->glidePlaneKey,pL).second;
//            assert( success && "COULD NOT INSERT GLIDE PLANE POINTER IN GLIDE PLANE OBSERVER.");
//        }
//
//        /**********************************************************************/
//        void removeGlidePlane(const GlidePlaneType* const pL)
//        {/*!@\param[in] pS a row pointer to a DislocationSegment
//          * Removes pS from *this GLidePlane
//          */
//            const int success=glidePlanes().erase(pL->glidePlaneKey);
//            assert(success==1 && "COULD NOT ERASE GLIDE PLANE POINTER FROM GLIDE PLANE OBSERVER.");
//        }

    };
}
#endif






//        /**********************************************************************/
//        template <class T>
//        friend T& operator << (T& os, const GlidePlaneFactoryType& gpo)
//        {
//            for (const auto& glidePlane : gpo.glidePlanes())
//            {
//                if(glidePlane.second->glissileLoopIDs.size())
//                {
//                    os << (*glidePlane.second);
//                }
//            }
//            return os;
//        }

//        /**********************************************************************/
//        static GlidePlaneKeyType getGlidePlaneKey(const Grain<dim>& grain,
//                                                  const VectorDimD& P,
//                                                  const VectorDimD& N)
//        {/*!\param[in] grain the grain on which the GlidePlane is defined
//          * \param[in] P a point on the plane
//          * \param[in] N the normal to the plane
//          * \returns the key which uniquely identifies the plane.
//          * The type of the key is a tuple with entries (grainID,r,h), where r
//          * is the ReciprocalLatticeDirection corresponding to N, and h=P.dot(r)
//          * is an integer indicating the "heigth" of the plane from the origin,
//          * in integer multiples of the interplanar distance d=1/|r|.
//          */
//            const ReciprocalLatticeDirection<dim> r(grain.reciprocalLatticeDirection(N));
//            return (GlidePlaneKeyType()<<grain.grainID,r,LatticePlane::height(LatticePlane::computeHeight(r,P))).finished();
//        }

//        /**********************************************************************/
//        static int sign(const long int& i)
//        {
//            if(i>0)
//            {
//                return 1;
//            }
//            else if(i<0)
//            {
//                return -1;
//            }
//            else
//            {
//                return 0;
//            }
//
//        }

//        /**********************************************************************/
//        static GlidePlaneKeyType getGlidePlaneKey(const Lattice<dim>& lattice,
//                                                  const int& grainID1,
////                                                  const int& grainID2,
//                                                  const VectorDimD& P,
//                                                  const VectorDimD& N)
//        {/*!\param[in] grain the grain on which the GlidePlane is defined
//          * \param[in] P a point on the plane
//          * \param[in] N the normal to the plane
//          * \returns the key which uniquely identifies the plane.
//          * The type of the key is a tuple with entries (grainID,r,h), where r
//          * is the ReciprocalLatticeDirection corresponding to N, and h=P.dot(r)
//          * is an integer indicating the "heigth" of the plane from the origin,
//          * in integer multiples of the interplanar distance d=1/|r|.
//          */
//            const ReciprocalLatticeDirection<dim> r(lattice.reciprocalLatticeDirection(N));
//            //            return (GlidePlaneKeyType()<<grainID1,grainID2,r,LatticePlane::height(LatticePlane::computeHeight(r,P))).finished();
//            //            return (GlidePlaneKeyType()<<grainID1,grainID2,r,LatticePlane::height(r,P)).finished();
//            const long int h(LatticePlane::height(r,P));
//            GlidePlaneKeyType temp;
//            temp[0]=grainID1;
//            temp[1]=grainID1;
//            const int signh(sign(h));
//            for(d=0;d<dim;++d)
//            {
//                temp[2+d]=r(d)*signh;
//            }
//            temp[2+dim]=h*signh;
//            return temp;
////            return (GlidePlaneKeyType()<<grainID1,grainID2,r*sign(h),h*sign(h)).finished(); // make sure that key heigh is always positive for uniqueness
//
//        }
//
//        /**********************************************************************/
//        static GlidePlaneKeyType getGlidePlaneKey(const int& grainID1,
//                                                  const int& grainID2)
//        {/*!\param[in] grain the grain on which the GlidePlane is defined
//          * \param[in] P a point on the plane
//          * \param[in] N the normal to the plane
//          * \returns the key which uniquely identifies the plane.
//          * The type of the key is a tuple with entries (grainID,r,h), where r
//          * is the ReciprocalLatticeDirection corresponding to N, and h=P.dot(r)
//          * is an integer indicating the "heigth" of the plane from the origin,
//          * in integer multiples of the interplanar distance d=1/|r|.
//          */
////            const ReciprocalLatticeDirection<dim> r(lattice.reciprocalLatticeDirection(N));
//            //            return (GlidePlaneKeyType()<<grainID1,grainID2,r,LatticePlane::height(LatticePlane::computeHeight(r,P))).finished();
//            //            return (GlidePlaneKeyType()<<grainID1,grainID2,r,LatticePlane::height(r,P)).finished();
////            const long int h(LatticePlane::height(r,P));
//            GlidePlaneKeyType temp;
//            temp[0]=grainID1;
//            temp[1]=grainID1;
////            const int signh(sign(h));
//            for(d=0;d<dim;++d)
//            {
//                temp[2+d]=0;
//            }
//            temp[2+dim]=0;
//            return temp;
//            //            return (GlidePlaneKeyType()<<grainID1,grainID2,r*sign(h),h*sign(h)).finished(); // make sure that key heigh is always positive for uniqueness
//
//        }

//        /**********************************************************************/
//        static GlidePlaneKeyType getGlidePlaneKey(const Lattice<dim>& lattice,
//                                                  const int& grainID1,
//                                                  const int& grainID2,
//                                                  const VectorDimD& P,
//                                                  const VectorDimD& N)
//        {/*!\param[in] grain the grain on which the GlidePlane is defined
//          * \param[in] P a point on the plane
//          * \param[in] N the normal to the plane
//          * \returns the key which uniquely identifies the plane.
//          * The type of the key is a tuple with entries (grainID,r,h), where r
//          * is the ReciprocalLatticeDirection corresponding to N, and h=P.dot(r)
//          * is an integer indicating the "heigth" of the plane from the origin,
//          * in integer multiples of the interplanar distance d=1/|r|.
//          */
//            const ReciprocalLatticeDirection<dim> r(lattice.reciprocalLatticeDirection(N));
//            //            return (GlidePlaneKeyType()<<grainID1,grainID2,r,LatticePlane::height(LatticePlane::computeHeight(r,P))).finished();
//            //            return (GlidePlaneKeyType()<<grainID1,grainID2,r,LatticePlane::height(r,P)).finished();
//            const long int h(LatticePlane::height(r,P));
//            return (GlidePlaneKeyType()<<grainID1,grainID2,r*sign(h),h*sign(h)).finished(); // make sure that key heigh is always positive for uniqueness
//
//        }




//
//        /**********************************************************************/
//        static std::deque<std::pair<VectorDimD,VectorDimD>,Eigen::aligned_allocator<std::pair<VectorDimD,VectorDimD>>> planeSegmentIntersection(const VectorDimD& P0,
//                                                                                                                                                const VectorDimD& N,
//                                                                                                                                                const VectorDimD& v0,
//                                                                                                                                                const VectorDimD& v1)
//        {
//            //            std::deque<VectorDim> temp;
//
//            const double nNorm(N.norm());
//            assert(nNorm>FLT_EPSILON);
//            const VectorDimD n(N/nNorm);
//
//            std::deque<std::pair<VectorDimD,VectorDimD>,Eigen::aligned_allocator<std::pair<VectorDimD,VectorDimD>>> temp;
//
//
//            // check intersection of v0->v1 with plane
//            // x=v0+u(v1-v0)
//            // (x-P0).n=0
//            // (v0+u(v1-v0)-P0).n=0
//            // u=(P0-v0).n/(v1-v0).n;
//            const double edgeNorm=(v1-v0).norm();
//            if(edgeNorm<FLT_EPSILON)
//            {
//                VectorDimD
//
//            }
//            else
//            {
//                const double den=(v1-v0).dot(n);
//                const double num=(P0-v0).dot(n);
//                const double P0v0norm=(P0-v0).norm();
//
//                const double numCheck= (P0v0norm<FLT_EPSILON)? 0.0 : num/P0v0norm;
//
//                if (fabs(den/edgeNorm)>FLT_EPSILON)
//                {
//                    // edge intersects plane
//                    const double u=num/den;
//
//                    if(fabs(u)<FLT_EPSILON)
//                    {
//                        temp.emplace_back(v0,v0);
//                    }
//                    else if (u>=FLT_EPSILON && u<=1.0-FLT_EPSILON)
//                    {
//                        const VectorDimD x((1.0-u)*v0 + u*v1);
//                        temp.emplace_back(x,x);
//                    }
//                    else if (fabs(1.0-u)<FLT_EPSILON)
//                    {
//                        temp.emplace_back(v1,v1);
//                    }
//                    else
//                    {// no roots
//
//                    }
//
//                }
//                else
//                {
//                    if (fabs(numCheck)>FLT_EPSILON)
//                    {// edge is parallel to plane, no intersection
//
//                    }
//                    else
//                    {// edge is coplanar
//                        temp.emplace_back(v0,v1);
//                        //                        temp.emplace_back(v1,&edge.child(1));
//                    }
//                }
//            }
//
//
//
//
//            return temp;
//        }






//        /**********************************************************************/
//        static GlidePlaneKeyType getGlidePlaneKey(const int& grainID,
//                                                const LatticePlane& lp
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
//            //            const ReciprocalLatticeDirection<dim> r(lattice.reciprocalLatticeDirection(N));
//            //            //            return (GlidePlaneKeyType()<<grainID1,grainID2,r,LatticePlane::height(LatticePlane::computeHeight(r,P))).finished();
//            //            //            return (GlidePlaneKeyType()<<grainID1,grainID2,r,LatticePlane::height(r,P)).finished();
//            //            const long int h(LatticePlane::height(r,P));
//            GlidePlaneKeyType temp;
//            temp[0]=grainID;
//            temp[1]=grainID;
//            //            const int signh(sign(h));
//            for(int d=0;d<dim;++d)
//            {
//                temp[2+d]=lp.n(d);
//            }
//            temp[2+dim]=lp.h;
//            return temp;
//            //            return (GlidePlaneKeyType()<<grainID1,grainID2,r*sign(h),h*sign(h)).finished(); // make sure that key heigh is always positive for uniqueness
//
//        }

//        /**********************************************************************/
//        static GlidePlaneKeyType getGlidePlaneKey(const int& grainID1,
//                                                const int& grainID2)
//        {/*!\param[in] grain the grain on which the GlidePlane is defined
//          * \param[in] P a point on the plane
//          * \param[in] N the normal to the plane
//          * \returns the key which uniquely identifies the plane.
//          * The type of the key is a tuple with entries (grainID,r,h), where r
//          * is the ReciprocalLatticeDirection corresponding to N, and h=P.dot(r)
//          * is an integer indicating the "heigth" of the plane from the origin,
//          * in integer multiples of the interplanar distance d=1/|r|.
//          */
//            //            const ReciprocalLatticeDirection<dim> r(lattice.reciprocalLatticeDirection(N));
//            //            return (GlidePlaneKeyType()<<grainID1,grainID2,r,LatticePlane::height(LatticePlane::computeHeight(r,P))).finished();
//            //            return (GlidePlaneKeyType()<<grainID1,grainID2,r,LatticePlane::height(r,P)).finished();
//            //            const long int h(LatticePlane::height(r,P));
//            GlidePlaneKeyType temp;
//            temp[0]=grainID1;
//            temp[1]=grainID2;
//            //            const int signh(sign(h));
//            for(int d=0;d<dim;++d)
//            {
//                temp[2+d]=0;
//            }
//            temp[2+dim]=0;
//            return temp;
//            //            return (GlidePlaneKeyType()<<grainID1,grainID2,r*sign(h),h*sign(h)).finished(); // make sure that key heigh is always positive for uniqueness
//
//        }
