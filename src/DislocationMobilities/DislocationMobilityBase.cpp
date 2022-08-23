/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by David Rivera <drivera2@ucla.edu>.
 * Copyright (C) 2013 by Giacomo Po   <gpo@ucla.edu>.
 * Copyright (C) 2013 by Tamer Crosby <tcrosby@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _model_DislocationMobility_cpp_
#define _model_DislocationMobility_cpp_

#include <DislocationMobilityBase.h>

namespace model
{
    
 
        /**********************************************************************/
        StochasticForceGenerator::StochasticForceGenerator(const DDtraitsIO& traitsIO):
        			 stochasticForceSeed(TextFileParser(traitsIO.ddFile).readScalar<int>("stochasticForceSeed",true))
        			, seed(stochasticForceSeed<0?std::chrono::system_clock::now().time_since_epoch().count():stochasticForceSeed)
        {
         
#ifdef _OPENMP
            const size_t nThreads = omp_get_max_threads();
#else
            const size_t nThreads = 1;
#endif
            generators.resize(nThreads,T(seed));
            distributions.resize(nThreads,std::normal_distribution<double>(0.0,1.0));
        
        }
        
        /**********************************************************************/
        double StochasticForceGenerator::stochasticVelocity(const double& kB,
                               const double& T,
                               const double& B,
                               const double& L,
                               const double& dt) 
        {
#ifdef _OPENMP
            return distributions[omp_get_thread_num()](generators[omp_get_thread_num()])*sqrt(2.0*kB*T/B/L/dt);
#else
            return distributions[0](generators[0])*sqrt(2.0*kB*T/B/L/dt);
#endif
        }
 
    /**************************************************************************/


        
        /**********************************************************************/
        DislocationMobilityBase::DislocationMobilityBase(const std::string& name_in) :
        /* init */ name(name_in)
        {
            std::cout<<greenBoldColor<<"Creating DislocationMobility "<<this->sID<<" ("<<name<<")"<<defaultColor<<std::endl;
        }
        
        
    
    
}
#endif
