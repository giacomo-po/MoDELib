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
#include <GlidePlaneFactory.h>

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
    {

        typedef GlidePlane<dim> GlidePlaneType;
        typedef GlidePlaneFactory<dim> GlidePlaneFactoryType;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef typename TypeTraits<GlidePlaneType>::KeyType KeyType;
        typedef KeyType GlidePlaneKeyType;

        static int verboseGlidePlane;

        const  GlidePlaneFactoryType& glidePlaneFactory;
        const Grain<dim>& grain;
        const GlidePlaneKeyType key;

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
        }
        
        /**********************************************************************/
        GlidePlane(const GlidePlane<dim>& other) = delete;
        
        /**********************************************************************/
        ~GlidePlane()
        {
            VerboseGlidePlane(1,"Destroying GlidePlane "<<this->sID<<std::endl;);
        }

    };

    template <int dim>
    int GlidePlane<dim>::verboseGlidePlane=0;

}
#endif

