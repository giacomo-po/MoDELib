/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PlaneLineIntersection_H_
#define model_PlaneLineIntersection_H_

#include <tuple>
#include <map>
#include <Eigen/Dense>
#include <Eigen/Cholesky>

namespace model
{
    
    template <int dim>
    struct PlaneLineIntersection
    {
        
        enum IntersectionType {PARALLEL,COINCIDENT,INCIDENT};
        
        typedef Eigen::Matrix<double,dim,1> VectorDimD;
        
        typedef std::tuple<IntersectionType,VectorDimD,VectorDimD> SolutionType;
        const SolutionType sol;
        
        /**********************************************************************/
        static SolutionType findIntersection(const VectorDimD& pP,
                                             const VectorDimD& N,
                                             const VectorDimD& pL,
                                             const VectorDimD& D);
        
        
    public:
        
        const IntersectionType& type;
        const VectorDimD& P;
        const VectorDimD& d;
        
        /**********************************************************************/
        PlaneLineIntersection(const VectorDimD& pP,
                              const VectorDimD& N,
                              const VectorDimD& pL,
                              const VectorDimD& D);
        
    };
    
    /******************************************************************/
    /******************************************************************/
} /* namespace model */
#endif
