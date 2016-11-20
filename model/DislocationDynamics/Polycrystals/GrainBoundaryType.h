/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_GrainBoundaryType_H_
#define model_GrainBoundaryType_H_

#include <string>
#include <math.h>       /* fabs */
#include <cfloat>       /* FLT_EPSILON */

namespace model
{
    
    
    
    template <int dim>
    struct GrainBoundaryType 
    {
        
        typedef Eigen::Matrix<  double,dim,1> VectorDimD;

        const std::string name;
        const VectorDimD t;
        const VectorDimD n1;
        const VectorDimD n2;
        const double energyDensity;
        const double dislocationSpacing;
        const VectorDimD Burgers;

        GrainBoundaryType(const std::string& name_in,
                          const VectorDimD& t_in,
                          const VectorDimD& n1_in,
                          const VectorDimD& n2_in,
                          const double& energy_in,
                          const double& spacing_in,
                          const VectorDimD& b_in) :
        /* init list */ name(name_in),
        /* init list */ t(t_in),
        /* init list */ n1(n1_in),
        /* init list */ n2(n2_in),
        /* init list */ energyDensity(energy_in),
        /* init list */ dislocationSpacing(spacing_in),
        /* init list */ Burgers(b_in)
        {
            assert(fabs(t.dot(n1))<FLT_EPSILON);
            assert(fabs(t.dot(n2))<FLT_EPSILON);
            assert(energyDensity>=0.0);
            assert(dislocationSpacing>0.0);
        }
        
        
    };
    
} // end namespace
#endif

