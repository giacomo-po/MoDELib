/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by David Rivera <drivera2@ucla.edu>.
 * Copyright (C) 2013 by Giacomo Po   <gpo@ucla.edu>.
 * Copyright (C) 2013 by Tamer Crosby <tcrosby@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _model_GrainBoundaryMobility_h_
#define _model_GrainBoundaryMobility_h_

#include <iostream>
#include <cmath>
#include <assert.h>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <BCClattice.h>
#include <FCClattice.h>
 // defines mode::cout
#include <TerminalColors.h> // defines mode::cout

namespace model
{
    


    
    /**************************************************************************/
    /**************************************************************************/
    struct GrainBoundaryMobility
    {
        
        //! Boltzmann constant in [eV]
        static constexpr double kB=8.617e-5;

        
        typedef Eigen::Matrix<double,3,3> MatrixDim;
        typedef Eigen::Matrix<double,3,1> VectorDim;
        
        const double B0;
        const double B1;
        const double dH0;
        const double tauP;
        
        /**********************************************************************/
        template<typename DislocationMobilityType>
        GrainBoundaryMobility(const DislocationMobilityType& dm,
                              const double& spacing,
                              const double& gbE) :
        
        /* init */ B0(dm.B0*100.0),
        /* init */ B1(dm.B1*100.0),
        /* init */ dH0(0.0),
        /* init */ tauP(1.0)
        {/*! Empty constructor is required by constexpr
          */
        }
        
        /**********************************************************************/
        double velocity(const MatrixDim& S,
                        const VectorDim& b,
                        const VectorDim& , // xi
                        const VectorDim& n,
                        const double& T) const
        {
            const double tauB=std::fabs(b.transpose()*S*n);
            const double tau=tauB/b.norm();
            const double num(1.0-tau/tauP);
            double exponential=1.0;
            if(num>0.0)
            {
                exponential=exp(-dH0*num/(kB*T));
            }
            
            
            return tauB/(B0+B1*T)*exponential;
        }
        
        
        
    };
    
    
}

#endif


