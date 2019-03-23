/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_Plane_H_
#define model_Plane_H_

#include <cfloat>
#include <tuple>
//#include <map>
#include <Eigen/Dense>


namespace model
{
    
    template <int dim>
    struct Plane
    {
        
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        
        const VectorDim P;
        const VectorDim unitNormal;
        
        /**********************************************************************/
        Plane(const VectorDim& p,const VectorDim& n) :
        /* init */ P(p),
        /* init */ unitNormal(n.normalized())
        {
            assert(n.norm()>FLT_EPSILON && "Plane must have non-zero normal");
        }
        
        /**********************************************************************/
        bool contains(const VectorDim& P0) const
        {
            const double PP0((P-P0).norm());
            return PP0<FLT_EPSILON? true : (fabs((P0-P).dot(unitNormal))<FLT_EPSILON*PP0);
        }
        
        /**********************************************************************/
        VectorDim snapToPlane(const VectorDim& P0) const
        {
            return P0-(P0-P).dot(unitNormal)*unitNormal;
        }
        
        /**********************************************************************/
        double distanceTo(const VectorDim& P0) const
        {
            return fabs((P0-P).dot(unitNormal));
        }

        
    };
    
}
#endif
