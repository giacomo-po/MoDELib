/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_GrainBoundary_H_
#define model_GrainBoundary_H_

#include <cfloat> // FLT_EPSILON
#include <Eigen/Core>
#include <model/Math/RoundEigen.h>
#include <model/Mesh/SimplicialMesh.h>
#include <model/Mesh/MeshRegionObserver.h>
#include <model/DislocationDynamics/Polycrystals/Grain.h>

namespace model
{
    
    
    
    template <int dim>
    class GrainBoundary
    {
        
        typedef MeshRegionBoundary<Simplex<dim,dim-1> > MeshRegionBoundaryType;


//        typedef Eigen::Matrix<  double,dim,1> VectorDimD;
//        typedef Eigen::Matrix<double,dim,dim> MatrixDimD;
        
        

        
    public:
        
        
        const MeshRegionBoundaryType& regionBoundary;
        
        const Grain<dim>& grainFirst;
        const Grain<dim>& grainSecond;
        
        /**********************************************************************/
        GrainBoundary(const MeshRegionBoundaryType& regionbnd_in,
                      const Grain<dim>& grainFirst_in,
                      const Grain<dim>& grainSecond_in) :
        /* init */ regionBoundary(regionbnd_in),
        /* init */ grainFirst(grainFirst_in),
        /* init */ grainSecond(grainSecond_in)
        {
            model::cout<<"Creating GrainBoundary ("<<regionBoundary.regionBndID.first<<" "<<regionBoundary.regionBndID.second<<")"<<std::endl;
            
            
        }
        
    };
    
    
} // end namespace
#endif

