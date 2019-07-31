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
#include <KeyConstructableSharedPtrFactory.h>
#include <GlidePlaneKey.h>

namespace model
{
    
    
    
    
    template <int dim>
    struct GlidePlane; // class predeclaration
    
    template<int dim>
    struct GlidePlaneFactory;
    
 
    
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
    template<int dim>
    struct GlidePlaneFactory : public KeyConstructableWeakPtrFactory<GlidePlaneFactory<dim>>
    /*                      */,private std::map<std::pair<size_t,size_t>,PlanePlaneIntersection<dim>>
    {

        typedef GlidePlaneFactory<dim> GlidePlaneFactoryType;
        typedef GlidePlane<dim> GlidePlaneType;
        typedef MeshPlane<dim> MeshPlaneType;
        typedef Eigen::Matrix<double,dim,1> VectorDimD;
        typedef typename GlidePlaneType::GlidePlaneKeyType GlidePlaneKeyType;
        typedef KeyConstructableWeakPtrFactory<GlidePlaneFactory<dim>> GlidePlaneMapType;
        typedef std::shared_ptr<GlidePlaneType> GlidePlaneSharedPtrType;
        typedef PlanePlaneIntersection<dim> PlanePlaneIntersectionType;
        typedef std::map<std::pair<size_t,size_t>,PlanePlaneIntersectionType> MeshPlaneIntersectionContainerType;
        
    public:
        
        const Polycrystal<dim>& poly;
        
        GlidePlaneFactory(const Polycrystal<dim>& poly_in) :
        /* init */ poly(poly_in)
        {
            
        }
        
//        GlidePlaneFactory<dim>& glidePlaneFactory()
//        {
//            return *this;
//        }
//        
//        const GlidePlaneFactory<dim>& glidePlaneFactory() const
//        {
//            return *this;
//        }
        
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
        
    };
}
#endif



