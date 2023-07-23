/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DDtimeIntegrator_cpp_
#define model_DDtimeIntegrator_cpp_

#include <cfloat>
#include <DDtimeIntegrator.h>

namespace model
{

    /******************************************************************/
    DDtimeIntegrator::DDtimeIntegrator(const std::string& fileName):
    /*                  init                    */  dxMax(TextFileParser(fileName).readScalar<double>("dxMax", true))
    /*                  init                    */, shearWaveSpeedFraction(1.0e-7)
    /*                  init                    */, timeIntegrationMethod(TextFileParser(fileName).readScalar<int>("timeIntegrationMethod", true))
    /*                  init                    */, dtMax(TextFileParser(fileName).readScalar<double>("timeStep", true))
    /*                  init                    */, dpdMax(TextFileParser(fileName).readScalar<double>("dpdMax", true))
{

        if (dxMax < FLT_EPSILON)
        {
            throw std::runtime_error("dxMax must be > FLT_EPSILON.");
        }
        
        if (dtMax < FLT_EPSILON)
        {
            throw std::runtime_error("dtMax must be > FLT_EPSILON.");
        }

    }
    
//    template struct DDtimeIntegrator<0>;

} // end namespace
#endif

