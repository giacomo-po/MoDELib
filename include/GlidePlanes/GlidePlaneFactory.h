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
#include <LatticeModule.h>
#include <ReciprocalLatticeDirection.h>
#include <PlanePlaneIntersection.h>
#include <MeshPlane.h>
#include <CRTP.h>
#include <TypeTraits.h>
#include <Polycrystal.h>
#include <WeakPtrFactories.h>
#include <GlidePlaneKey.h>

namespace model
{
    template <int dim>
    struct GlidePlane; // class predeclaration
   
    template<int dim>
    struct GlidePlaneFactory : public KeyConstructableWeakPtrFactory<GlidePlaneFactory<dim>,GlidePlane<dim>>
    /*                      */,private std::map<std::pair<size_t,size_t>,PlanePlaneIntersection<dim>>
    {

        typedef GlidePlaneFactory<dim> GlidePlaneFactoryType;
        typedef GlidePlane<dim> GlidePlaneType;
        typedef MeshPlane<dim> MeshPlaneType;
        typedef Eigen::Matrix<double,dim,1> VectorDimD;
        typedef typename GlidePlaneType::GlidePlaneKeyType GlidePlaneKeyType;
        typedef KeyConstructableWeakPtrFactory<GlidePlaneFactory<dim>,GlidePlaneType> GlidePlaneMapType;
        typedef std::shared_ptr<GlidePlaneType> GlidePlaneSharedPtrType;
        typedef PlanePlaneIntersection<dim> PlanePlaneIntersectionType;
        typedef std::map<std::pair<size_t,size_t>,PlanePlaneIntersectionType> MeshPlaneIntersectionContainerType;
        
    public:
        
        const Polycrystal<dim>& poly;
        
        GlidePlaneFactory(const Polycrystal<dim>& poly_in);

        MeshPlaneIntersectionContainerType& glidePlaneIntersections();
        const MeshPlaneIntersectionContainerType& glidePlaneIntersections() const;
        const PlanePlaneIntersectionType& glidePlaneIntersection(const MeshPlaneType* const p1,
                                                                 const MeshPlaneType* const p2);
        GlidePlaneMapType& glidePlanes();
        const GlidePlaneMapType& glidePlanes() const;
        
    };
}
#endif
