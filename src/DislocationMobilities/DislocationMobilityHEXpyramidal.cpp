/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by David Rivera <drivera2@ucla.edu>.
 * Copyright (C) 2013 by Giacomo Po   <gpo@ucla.edu>.
 * Copyright (C) 2013 by Tamer Crosby <tcrosby@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _model_DislocationMobilityHEXpyramidal_cpp_
#define _model_DislocationMobilityHEXpyramidal_cpp_

#include <DislocationMobilityHEXpyramidal.h>

namespace model
{
    
    
    
        /**********************************************************************/
        DislocationMobilityHEXpyramidal::DislocationMobilityHEXpyramidal(const PolycrystallineMaterialBase& material) :
        /* init */ DislocationMobilityBase("HEX-pyramidal mobility for "+material.materialName),
        /* init */ B0e(TextFileParser(material.materialFile).readScalar<double>("B0e_SI_py",true)*material.cs_SI/(material.mu_SI*material.b_SI)),
        /* init */ B1e(TextFileParser(material.materialFile).readScalar<double>("B1e_SI_py",true)*material.cs_SI/(material.mu_SI*material.b_SI)),
        /* init */ B0s(TextFileParser(material.materialFile).readScalar<double>("B0s_SI_py",true)*material.cs_SI/(material.mu_SI*material.b_SI)),
        /* init */ B1s(TextFileParser(material.materialFile).readScalar<double>("B1s_SI_py",true)*material.cs_SI/(material.mu_SI*material.b_SI)),
        //        /* init */ kB(kB_SI/mu_SI/std::pow(b_SI,3))
        /* init */ kB(kB_SI/material.mu_SI/std::pow(material.b_SI,3))
        {
            
        }
        /**********************************************************************/
        double DislocationMobilityHEXpyramidal::velocity(const MatrixDim& S,
                        const VectorDim& b,
                        const VectorDim& xi,
                        const VectorDim& n,
                        const double& T,
                        const double& dL,
                        const double& dt,
                        const std::shared_ptr<StochasticForceGenerator>& sfg)
        {
            double ve=std::fabs(b.transpose()*S*n)/(B0e+B1e*T);
            double vs=std::fabs(b.transpose()*S*n)/(B0s+B1s*T);
            
            if(sfg)
            {
                ve+=sfg->stochasticVelocity(kB,T,B0e+B1e*T,dL,dt);
                vs+=sfg->stochasticVelocity(kB,T,B0s+B1s*T,dL,dt);
            }
            
            // Interpolate ve and vs
            // const double cos2=std::pow(b.normalized().dot(xi),2);
            // return vs*cos2+ve*(1.0-cos2);
            const double cos1=std::fabs(b.normalized().dot(xi));
            const double sgm1=1.0/(1.0+exp(-2.0*(cos1-0.7)/0.05));
            return vs*sgm1+ve*(1.0-sgm1);
        }
   
    
}
#endif
