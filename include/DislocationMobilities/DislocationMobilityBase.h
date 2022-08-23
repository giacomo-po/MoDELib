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
#include <DDtraitsIO.h>

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
        StochasticForceGenerator(const DDtraitsIO& traitsIO);
        
        /**********************************************************************/
        double stochasticVelocity(const double& kB,
                               const double& T,
                               const double& B,
                               const double& L,
                               const double& dt);
        
    };
    
    /**************************************************************************/
    /**************************************************************************/
    struct DislocationMobilityBase : public StaticID<DislocationMobilityBase>
//    				, public StochasticForceGenerator
    {
        typedef Eigen::Matrix<double,3,3> MatrixDim;
        typedef Eigen::Matrix<double,3,1> VectorDim;
        const std::string name;

        
        /**********************************************************************/
        DislocationMobilityBase(const std::string& name_in) ;
        
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
                                const std::shared_ptr<StochasticForceGenerator>& sfg) =0 ;
        
    };
    
}
#endif
