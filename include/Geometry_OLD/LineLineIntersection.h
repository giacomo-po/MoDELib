/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LineLineIntersection_H_
#define model_LineLineIntersection_H_

#include <tuple>
#include <map>
#include <stdexcept>
#include <Eigen/Dense>

namespace model
{
    
    template <int dim>
    class LineLineIntersection
    {
        
    public:
        enum IntersectionType {PARALLEL,COINCIDENT,INCIDENT,SKEW};
        
        
        typedef Eigen::Matrix<double,dim,1> VectorDimD;
        typedef std::tuple<IntersectionType,VectorDimD,VectorDimD> SolutionType;
        
    private:
        const SolutionType sol;
        
        /**********************************************************************/
        static SolutionType findIntersections(const VectorDimD& A0,
                                              const VectorDimD& D0,
                                              const VectorDimD& A1,
                                              const VectorDimD& D1);
        
    public:
        
        const IntersectionType& type;
        const VectorDimD& x0;
        const VectorDimD& x1;
        
        /**********************************************************************/
        LineLineIntersection(const VectorDimD& A0,
                             const VectorDimD& D0,
                             const VectorDimD& A1,
                             const VectorDimD& D1);
        
    };
    
    /******************************************************************/
    /******************************************************************/
} /* namespace model */
#endif
