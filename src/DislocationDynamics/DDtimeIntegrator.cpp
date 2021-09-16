/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DDtimeIntegrator_cpp_
#define model_DDtimeIntegrator_cpp_

#include <DDtimeIntegrator.h>

namespace model
{

    /******************************************************************/
    DDtimeIntegrator<0>::DDtimeIntegrator(const std::string& fileName): dxMax(TextFileParser(fileName).readScalar<double>("dxMax", true))
    /*                  init                    */, timeIntegrationMethod(TextFileParser(fileName).readScalar<int>("timeIntegrationMethod", true))
    /*                  init                    */, dt(TextFileParser(fileName).readScalar<double>("timeStep", true))
    {
        assert(dxMax > 0.0);
        if (timeIntegrationMethod == 0)
        {
            assert(dt > 0.0 && "Time step should be greater than zero for constant time stepping.");
        }
    }
    
    template struct DDtimeIntegrator<0>;

} // end namespace
#endif

