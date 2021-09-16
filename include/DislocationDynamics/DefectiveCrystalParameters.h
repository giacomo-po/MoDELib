/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2017 by Yinan Cui <cuiyinan@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DefectiveCrystalParameters_H_
#define model_DefectiveCrystalParameters_H_

#include <string>
#include <TextFileParser.h>
#include <IDreader.h>



namespace model
{
    
    struct DefectiveCrystalParameters
    {
        
        enum SimulationType{FINITE_NO_FEM=0,FINITE_FEM=1,PERIODIC_IMAGES=2,PERIODIC_FEM=3};

        
        const int simulationType;
        const bool useDislocations;
        const bool useCracks;

        const int periodicImages_x;
        const int periodicImages_y;
        const int periodicImages_z;
        const long int Nsteps;
        const int timeIntegrationMethod;

        const std::string externalLoadControllerName;
        const double virtualSegmentDistance;

        
        long int runID;
        double totalTime;
        double dt;
        
        /**********************************************************************/
        void manageRestart();
        
    public:
        
        
        /**********************************************************************/
        DefectiveCrystalParameters(int& , char* []) ;
        
        bool isPeriodicSimulation() const;
    };
}
#endif
