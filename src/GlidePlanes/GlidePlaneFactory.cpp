/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2011 by Mamdouh Mohamed <mamdouh.s.mohamed@gmail.com>,
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_GlidePlaneFactory_cpp_
#define model_GlidePlaneFactory_cpp_

#include <GlidePlaneFactory.h>
#include <GlidePlane.h>

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

    template<int dim>
    GlidePlaneFactory<dim>::GlidePlaneFactory(const Polycrystal<dim>& poly_in) :
    /* init */ poly(poly_in)
    {
        
    }

    template<int dim>
    typename GlidePlaneFactory<dim>::MeshPlaneIntersectionContainerType& GlidePlaneFactory<dim>::glidePlaneIntersections()
    {
        return *this;
    }

    template<int dim>
    const typename GlidePlaneFactory<dim>::MeshPlaneIntersectionContainerType& GlidePlaneFactory<dim>::glidePlaneIntersections() const
    {
        return *this;
    }

    template<int dim>
    const typename GlidePlaneFactory<dim>::PlanePlaneIntersectionType& GlidePlaneFactory<dim>::glidePlaneIntersection(const MeshPlaneType* const p1,
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
            const bool success(glidePlaneIntersections().emplace(std::piecewise_construct,
                                                                 std::make_tuple(key),
                                                                 std::make_tuple(p1->P,p1->unitNormal,p2->P,p2->unitNormal)
                                                                 ).second);
            if(!success)
            {
                throw std::runtime_error("GlidePlaneFactory::glidePlaneIntersection not found.");
            }
        }
        return glidePlaneIntersections().at(key);
    }

    template<int dim>
    typename GlidePlaneFactory<dim>::GlidePlaneMapType& GlidePlaneFactory<dim>::glidePlanes()
    {
        return *this;
    }

    template<int dim>
    const typename GlidePlaneFactory<dim>::GlidePlaneMapType& GlidePlaneFactory<dim>::glidePlanes() const
    {
        return *this;
    }

template struct GlidePlaneFactory<3>;

}
#endif



