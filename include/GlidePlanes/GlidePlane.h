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

//#include <PlaneMeshIntersection.h>
#include <MeshPlane.h>
#include <GlidePlaneFactory.h>


#ifndef NDEBUG
#define VerboseGlidePlane(N,x) if(verboseGlidePlane>=N){std::cout<<x;}
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
        typedef GlidePlaneKey<dim> KeyType;
        typedef KeyType GlidePlaneKeyType;

        static int verboseGlidePlane;

        const  GlidePlaneFactoryType& glidePlaneFactory;
        const Grain<dim>& grain;
        const GlidePlaneKeyType key;

        /**********************************************************************/
        GlidePlane(const GlidePlaneFactoryType* const gpF,
                   const GlidePlaneKeyType& key_in) ;
        GlidePlane(const GlidePlane<dim>& other) = delete;
        ~GlidePlane();
        std::vector<std::shared_ptr<SlipSystem>> slipSystems() const;

    };

    template <int dim>
    int GlidePlane<dim>::verboseGlidePlane=0;

}
#endif

