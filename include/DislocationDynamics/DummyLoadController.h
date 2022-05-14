/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2012 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _model_DummyLoadController_h
#define _model_DummyLoadController_h

#include <Eigen/Dense>
#include <TerminalColors.h>



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
        void addDirichletConditions(const DislocationNetworkType& ) const
        {
        }
        
        /**********************************************************************/
        template <typename DislocationNetworkType>
        void init(const DislocationNetworkType&) const
        {
        }
        
        /**********************************************************************/
        template <typename DislocationNetworkType>
        void update(const DislocationNetworkType& ) const
        {
        }
        

        /**********************************************************************/
        template <typename DislocationNetworkType>
        Eigen::VectorXd globalVector(const DislocationNetworkType&) const
        {
            return Eigen::VectorXd();
        }
        
        /**************************************************************************/
        template <typename DislocationNetworkType>
        void output(const DislocationNetworkType& ,
                    const long int& ,
                    std::ofstream& ,
                    std::ofstream& ) const
        {
            
        }
        
    };
    
} // end namespace
#endif
