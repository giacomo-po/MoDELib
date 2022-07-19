/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpF@ucla.edu>.
 * Copyright (C) 2011 by Mamdouh Mohamed <mamdouh.s.mohamed@gmail.com>,
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_GLIDEPLANE_cpp
#define model_GLIDEPLANE_cpp

#include <GlidePlane.h>

namespace model
{

        /**********************************************************************/
        template <int dim>
        GlidePlane<dim>::GlidePlane(const GlidePlaneFactoryType* const gpF,
                   const GlidePlaneKeyType& key_in) :
        /* init */ LatticePlane(key_in.planeIndex(),ReciprocalLatticeDirection<dim>(key_in.reciprocalDirectionComponents(),*gpF->poly.grain(key_in.latticeID()).singleCrystal)) // BETTER TO CONSTRUCT N WITH PRIMITIVE VECTORS ON THE PLANE
        /* init */,MeshPlane<dim>(gpF->poly.mesh,key_in.latticeID(),this->planeOrigin(),this->n.cartesian())
        /* init */,glidePlaneFactory(*gpF)
        /* init */,grain(gpF->poly.grain(key_in.latticeID()))
        /* init */,key(key_in)
        {
            VerboseGlidePlane(1,"Creating GlidePlane "<<this->sID<<std::endl;);
        }
        
        /**********************************************************************/
        
        /**********************************************************************/
        template <int dim>
        GlidePlane<dim>::~GlidePlane()
        {
            VerboseGlidePlane(1,"Destroying GlidePlane "<<this->sID<<std::endl;);
        }

template <int dim>
std::vector<std::shared_ptr<SlipSystem>> GlidePlane<dim>::slipSystems() const
{
    std::vector<std::shared_ptr<SlipSystem>> temp;
    for(const auto& ss : grain.singleCrystal->slipSystems())
    {
        if(this->n.cross(ss->n).squaredNorm()==0)
        {
            temp.push_back(ss);
        }
    }
    
    return temp;
    
}



    template struct GlidePlane<3>;
}
#endif

