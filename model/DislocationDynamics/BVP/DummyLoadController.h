/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2012 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _model_DummyLoadController_h
#define _model_DummyLoadController_h

#include <model/Utilities/TerminalColors.h>
#include <model/MPI/MPIcout.h>


namespace model
{
    
    /**************************************************************************/
    template<typename TrialFunctionType>
    struct LoadController
    {
        
        LoadController(const TrialFunctionType&)
        {
        
        }
        
        /**********************************************************************/
        template <typename DislocationNetworkType>
        void init(const DislocationNetworkType&)
        {
            std::cout<<greenColor<<"Initializing DummyLoadController"<<defaultColor<<std::endl;
        }
        
    };
    
} // end namespace
#endif
