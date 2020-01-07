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
#include <assert.h>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <MPIcout.h> // defines model::cout
#include <TerminalColors.h> // defines mode::cout
#include <DislocatedMaterialBase.h>

namespace model
{
    
    /**************************************************************************/
    /**************************************************************************/
    struct StochasticForceGenerator
    {
        typedef std::default_random_engine T;
        static std::vector<T> generators;
        static std::vector<std::normal_distribution<double>> distributions;
        
        /**********************************************************************/
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
        
        /**********************************************************************/
        static double velocity(const double& kB,
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
    
#ifndef _MODEL_GREATWHITE_
    std::vector<std::default_random_engine> StochasticForceGenerator::generators;
    std::vector<std::normal_distribution<double>> StochasticForceGenerator::distributions;
#endif
    
    /**************************************************************************/
    /**************************************************************************/
    struct DislocationMobilityBase : public StaticID<DislocationMobilityBase>
    {
        typedef Eigen::Matrix<double,3,3> MatrixDim;
        typedef Eigen::Matrix<double,3,1> VectorDim;
        const std::string name;
        
        /**********************************************************************/
        DislocationMobilityBase(const std::string& name_in) :
        /* init */ name(name_in)
        {
            model::cout<<greenBoldColor<<"Creating DislocationMobility "<<this->sID<<" ("<<name<<")"<<defaultColor<<std::endl;
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
                                const bool& use_stochasticForce) const=0 ;
        
    };
    
}
#endif
