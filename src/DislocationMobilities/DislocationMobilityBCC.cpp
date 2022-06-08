/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by David Rivera <drivera2@ucla.edu>.
 * Copyright (C) 2013 by Giacomo Po   <gpo@ucla.edu>.
 * Copyright (C) 2013 by Tamer Crosby <tcrosby@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _model_DislocationMobilityBCC_cpp_
#define _model_DislocationMobilityBCC_cpp_

#include <numbers>

#include <DislocationMobilityBCC.h>

namespace model
{
    
        /**********************************************************************/
        DislocationMobilityBCC::DislocationMobilityBCC(const double& b_SI,
                            const double& mu_SI,
                            const double& cs_SI,
                            const double& B0e_SI, const double& B1e_SI,
                            const double& B0s_SI, const double& B1s_SI,
                            const double& Bk_SI,
                            const double& dH0_SI,
                            const double& p_in,
                            const double& q_in,
                            const double& T0_in,
                            const double& tauC_in,
                            const double& a0_in,
                            const double& a1_in,
                            const double& a2_in,
                            const double& a3_in) :
        /* init */ DislocationMobilityBase("BCC DislocationMobility"),
        /* init */ h(2.0*sqrt(2.0)/3.0), // units of b
        /* init */ w(25.0), // units of b
        /* init */ B0e(B0e_SI*cs_SI/(mu_SI*b_SI)),
        /* init */ B1e(B1e_SI*cs_SI/(mu_SI*b_SI)),
        /* init */ B0s(B0s_SI*cs_SI/(mu_SI*b_SI)),
        /* init */ B1s(B1s_SI*cs_SI/(mu_SI*b_SI)),
        /* init */ Bk(  Bk_SI*cs_SI/(mu_SI*b_SI)),
        /* init */ dH0(dH0_SI),
        /* init */ p(p_in),
        /* init */ q(q_in),
        /* init */ T0(T0_in),
        /* init */ tauC(tauC_in/mu_SI),
        /* init */ a0(a0_in),
        /* init */ a1(a1_in),
        /* init */ a2(a2_in),
        /* init */ a3(a3_in),
        /* init */ kB(kB_SI/mu_SI/std::pow(b_SI,3))
        {/*! Empty constructor is required by constexpr
          */
            
        }

        
        /**********************************************************************/
        DislocationMobilityBCC::DislocationMobilityBCC(const PolycrystallineMaterialBase& material) :
                /* init */ DislocationMobilityBase("BCC mobility for "+material.materialName),
        /* init */ h(2.0*sqrt(2.0)/3.0), // units of b
        /* init */ w(25.0), // units of b
        /* init */ B0e(TextFileParser(material.materialFile).readScalar<double>("B0e_SI",true)*material.cs_SI/(material.mu_SI*material.b_SI)),
        /* init */ B1e(TextFileParser(material.materialFile).readScalar<double>("B1e_SI",true)*material.cs_SI/(material.mu_SI*material.b_SI)),
        /* init */ B0s(TextFileParser(material.materialFile).readScalar<double>("B0s_SI",true)*material.cs_SI/(material.mu_SI*material.b_SI)),
        /* init */ B1s(TextFileParser(material.materialFile).readScalar<double>("B1s_SI",true)*material.cs_SI/(material.mu_SI*material.b_SI)),
        /* init */ Bk (TextFileParser(material.materialFile).readScalar<double>("Bk_SI", true)*material.cs_SI/(material.mu_SI*material.b_SI)),
        /* init */ dH0(TextFileParser(material.materialFile).readScalar<double>("dH0_eV",true)),
        /* init */ p(TextFileParser(material.materialFile).readScalar<double>("p",true)),
        /* init */ q(TextFileParser(material.materialFile).readScalar<double>("q",true)),
        /* init */ T0(TextFileParser(material.materialFile).readScalar<double>("Tf",true)*material.Tm),
        /* init */ tauC(TextFileParser(material.materialFile).readScalar<double>("tauC_SI",true)/material.mu_SI),
        /* init */ a0(TextFileParser(material.materialFile).readScalar<double>("a0",true)),
        /* init */ a1(TextFileParser(material.materialFile).readScalar<double>("a1",true)),
        /* init */ a2(TextFileParser(material.materialFile).readScalar<double>("a2",true)),
        /* init */ a3(TextFileParser(material.materialFile).readScalar<double>("a3",true)),
        /* init */ kB(kB_SI/material.mu_SI/std::pow(material.b_SI,3))
        {/*! Empty constructor is required by constexpr
          */

        }
        

        /**********************************************************************/
        double DislocationMobilityBCC::sigmoid(const double & x)
        {
            return 2.0/(1.0+exp(2.0*x));
        }
        
        /**********************************************************************/
        double DislocationMobilityBCC::velocity(const MatrixDim& S,
                        const VectorDim& b,
                        const VectorDim& xi,
                        const VectorDim& n,
                        const double& T,
                        const double& dL,
                        const double& dt,
                        const std::shared_ptr<StochasticForceGenerator>& sfg)
        {
            
            const double bNorm=b.norm();
            const VectorDim s = b/bNorm;
            const VectorDim n1 = Eigen::AngleAxisd(std::numbers::pi/3.0,s)*n;
            
            // Compute components of non-Schmid model
            const double tau=s.transpose()*S*n; // magnitude of resolved shear stress
            const double tauOrt=n.cross(s).transpose()*S*n;
            const double tau1=s.transpose()*S*n1; // resolved schear stress on
            const double tauOrt1=n1.cross(s).transpose()*S*n1;
            
            const double num=std::fabs(tau+a1*tau1);
            const double den=a0*tauC*sigmoid((a2*tauOrt+a3*tauOrt1)/a0/tauC);

            assert(den>0.0 && "den must be > 0.");
            
            const double Theta=num/den;
            const double dg = (Theta<1.0)? (std::pow(1.0-std::pow(Theta,p),q)-T/T0) : 0.0;
            const double dg1 = (dg>0.0)? dg : 0.0;
            const double expCoeff = exp(-dH0*dg1/(2.0*kB_eV*T));
            
            // Compute screw drag coeff
            const double sgm=0.5*sigmoid(-0.5*(0.05-dg1)/0.05);
            const double Bs=Bk*w/(2.0*h)*(1.0-sgm)+(B0s+B1s*T)*sgm; //kink-dominated to drag-dominated interpolation
            
            // Compute screw velocity
            double vs=std::fabs(tau)*bNorm/Bs*expCoeff;
            
            // Compute edge velocity
            double ve=std::fabs(tau)*bNorm/(B0e+B1e*T);
            
            if(sfg)
            {
                vs+=sfg->stochasticVelocity(kB,T,Bs,dL,dt);
                ve+=sfg->stochasticVelocity(kB,T,B0e+B1e*T,dL,dt);
            }
            
            // Interpolate ve and vs
            const double cos2=std::pow(s.dot(xi),2);
            return vs*cos2+ve*(1.0-cos2);
        }

    
}
#endif
