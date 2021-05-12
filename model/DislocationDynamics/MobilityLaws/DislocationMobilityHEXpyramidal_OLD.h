/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by David Rivera <drivera2@ucla.edu>.
 * Copyright (C) 2013 by Giacomo Po   <gpo@ucla.edu>.
 * Copyright (C) 2013 by Tamer Crosby <tcrosby@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _model_DislocationMobilityHEXpyramidal_h_
#define _model_DislocationMobilityHEXpyramidal_h_

#include <DislocationMobilityBase.h>

namespace model
{
    
    struct DislocationMobilityHEXpyramidal : public DislocationMobilityBase
    {
        
        typedef Eigen::Matrix<double,3,3> MatrixDim;
        typedef Eigen::Matrix<double,3,1> VectorDim;
        
        static constexpr double kB_SI=1.38064852e-23; // [J/K]
        
        const double B0e;
        const double B1e;
        const double B0s;
        const double B1s;
        const double kB;
        
        /**********************************************************************/
        DislocationMobilityHEXpyramidal(const DislocatedMaterialBase& material) :
        /* init */ DislocationMobilityBase("HEX-pyramidal mobility for "+material.materialName),
        /* init */ B0e(TextFileParser(material.materialFile).readScalar<double>("B0e_SI",true)*material.cs_SI/(material.mu_SI*material.b_SI)),
        /* init */ B1e(TextFileParser(material.materialFile).readScalar<double>("B1e_SI",true)*material.cs_SI/(material.mu_SI*material.b_SI)),
        /* init */ B0s(TextFileParser(material.materialFile).readScalar<double>("B0s_SI",true)*material.cs_SI/(material.mu_SI*material.b_SI)),
        /* init */ B1s(TextFileParser(material.materialFile).readScalar<double>("B1s_SI",true)*material.cs_SI/(material.mu_SI*material.b_SI)),
        //        /* init */ kB(kB_SI/mu_SI/std::pow(b_SI,3))
        /* init */ kB(kB_SI/material.mu_SI/std::pow(material.b_SI,3))
        {
            
        }
        /**********************************************************************/
        double velocity(const MatrixDim& S,
                        const VectorDim& b,
                        const VectorDim& xi,
                        const VectorDim& n,
                        const double& T,
                        const double& dL,
                        const double& dt,
                        const bool& use_stochasticForce) const override
        {
            double ve=std::fabs(b.transpose()*S*n)/(B0e+B1e*T);
            double vs=std::fabs(b.transpose()*S*n)/(B0s+B1s*T);
            
            if(use_stochasticForce)
            {
                ve+=StochasticForceGenerator::velocity(kB,T,B0e+B1e*T,dL,dt);
                vs+=StochasticForceGenerator::velocity(kB,T,B0s+B1s*T,dL,dt);
            }
            
            // Interpolate ve and vs
            const double cos2=std::pow(b.normalized().dot(xi),2);
            return vs*cos2+ve*(1.0-cos2);
        }
        
    };
    
}
#endif
