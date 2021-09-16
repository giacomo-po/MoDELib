/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by David Rivera <drivera2@ucla.edu>.
 * Copyright (C) 2013 by Giacomo Po   <gpo@ucla.edu>.
 * Copyright (C) 2013 by Tamer Crosby <tcrosby@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _model_DislocationMobility_h_
#define _model_DislocationMobility_h_

#include <iostream>
#include <random>
#include <cmath>
#include <vector>
#include <chrono>
#include <assert.h>
#include <Eigen/Dense>
#include <Eigen/Geometry>
 // defines std::cout
#include <TerminalColors.h> // defines mode::cout
#include <PolycrystallineMaterialBase.h>
#include <StaticID.h>

namespace model
{
    
    /**************************************************************************/
    /**************************************************************************/
    struct StochasticForceGenerator
    {
        typedef std::default_random_engine T;
        std::vector<T> generators;
        std::vector<std::normal_distribution<double>> distributions;
        const int stochasticForceSeed;
        const int seed;
        /**********************************************************************/
        StochasticForceGenerator():
        			 stochasticForceSeed(TextFileParser("inputFiles/DD.txt").readScalar<int>("stochasticForceSeed",true))
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
/*
        static void init(const int& seed)
        {
#ifdef _OPENMP
            const size_t nThreads = omp_get_max_threads();
#else
            const size_t nThreads = 1;
#endif
            generators.resize(nThreads,T(seed));
            distributions.resize(nThreads,std::normal_distribution<double>(0.0,1.0));
        }
*/
        
        /**********************************************************************/
        double stochasticVelocity(const double& kB,
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
        
    };
    
/*
#ifndef _MODEL_GREATWHITE_
    std::vector<std::default_random_engine> StochasticForceGenerator::generators;
    std::vector<std::normal_distribution<double>> StochasticForceGenerator::distributions;
#else
#ifdef _MODEL_GREATWHITE_standalone
    std::vector<std::default_random_engine> StochasticForceGenerator::generators;
    std::vector<std::normal_distribution<double>> StochasticForceGenerator::distributions;
#endif

#endif
  */  
    /**************************************************************************/
    /**************************************************************************/
    struct DislocationMobilityBase : public StaticID<DislocationMobilityBase>
    				, public StochasticForceGenerator
    {
        typedef Eigen::Matrix<double,3,3> MatrixDim;
        typedef Eigen::Matrix<double,3,1> VectorDim;
        const std::string name;

        
        /**********************************************************************/
        DislocationMobilityBase(const std::string& name_in) :
        /* init */ name(name_in)
        {
            std::cout<<greenBoldColor<<"Creating DislocationMobility "<<this->sID<<" ("<<name<<")"<<defaultColor<<std::endl;
        }
        
        /**********************************************************************/
        virtual ~DislocationMobilityBase(){};
        
        /**********************************************************************/
        virtual double velocity(const MatrixDim& S,
                                const VectorDim& b,
                                const VectorDim& , // xi
                                const VectorDim& n,
                                const double& T,
                                const double& dL,
                                const double& dt,
                                const bool& use_stochasticForce) =0 ;
        
    };
    
}
#endif
