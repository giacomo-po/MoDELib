/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PlanePlaneIntersection_H_
#define model_PlanePlaneIntersection_H_

#include <tuple>
#include <map>
#include <Eigen/Dense>
#include <Eigen/Cholesky>
#include <Plane.h>

namespace model
{
    
    template <int dim>
    struct PlanePlaneIntersection
    {
        enum IntersectionType {PARALLEL,COINCIDENT,INCIDENT};
        typedef Eigen::Matrix<double,dim,1> VectorDimD;
        typedef std::tuple<IntersectionType,VectorDimD,VectorDimD> SolutionType;
        
        
        /**********************************************************************/
        static bool checkIntersection(const VectorDimD& C,
                                      const VectorDimD& p1,
                                      const VectorDimD& n1,
                                      const VectorDimD& p2,
                                      const VectorDimD& n2,
                                      const double& tol);
        
        /**********************************************************************/
        static SolutionType findIntersection(const VectorDimD& p1,
                                             const VectorDimD& N1,
                                             const VectorDimD& p2,
                                             const VectorDimD& N2,
                                             const double& tol);
        
        
    public:
        
        
        const SolutionType sol;
        const IntersectionType& type;
        const VectorDimD& P;
        const VectorDimD& d;
        
        /**********************************************************************/
        PlanePlaneIntersection(const VectorDimD& p1,
                              const VectorDimD& n1,
                              const VectorDimD& p2,
                              const VectorDimD& n2,
                               const double& tolerance=FLT_EPSILON);

        /**********************************************************************/
        PlanePlaneIntersection(const Plane<dim>& plane0,
                              const Plane<dim>& plane1,
                               const double& tolerance=FLT_EPSILON);
        
    };
    
}
#endif
