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
#include <BCClattice.h>
#include <FCClattice.h>
#include <MPIcout.h> // defines mode::cout
#include <TerminalColors.h> // defines mode::cout
#include <Material.h>

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
    
    std::vector<std::default_random_engine> StochasticForceGenerator::generators;
    std::vector<std::normal_distribution<double>> StochasticForceGenerator::distributions;
    
    /**************************************************************************/
    /**************************************************************************/
    struct DislocationMobilityBase
    {
        typedef Eigen::Matrix<double,3,3> MatrixDim;
        typedef Eigen::Matrix<double,3,1> VectorDim;
        
        DislocationMobilityBase(const std::string& name)
        {
            model::cout<<greenBoldColor<<"Creating "<<name<<defaultColor<<std::endl;
        }
        
        virtual double velocity(const MatrixDim& S,
                                const VectorDim& b,
                                const VectorDim& , // xi
                                const VectorDim& n,
                                const double& T,
                                const double& dL,
                                const double& dt,
                                const bool& use_stochasticForce) const=0 ;
        
        //        virtual ~DislocationMobilityBase(){};
    };
    
    /**************************************************************************/
    /**************************************************************************/
    template <typename CrystalStructure>
    struct DislocationMobility : DislocationMobilityBase
    {
        //        static_assert(false, "CrystalStructure must be FCC or BCC");
    };
    
    template <>
    struct DislocationMobility<FCClattice<3>> : DislocationMobilityBase
    {
        
        typedef Eigen::Matrix<double,3,3> MatrixDim;
        typedef Eigen::Matrix<double,3,1> VectorDim;
        
        static constexpr double kB_SI=1.38064852e-23; // [J/K]
        
        const double B0e;
        const double B1e;
        const double B0s;
        const double B1s;
        const double kB;
        
        
        
        
//        /**********************************************************************/
//        DislocationMobility(const double& b_SI,
//                            const double& mu_SI,
//                            const double& cs_SI,
//                            const double& B1e_SI,
//                            const double& B1s_SI) :
//                /* init */ DislocationMobilityBase("FCC DislocationMobility"),
//        /* init */ B0(0.0),
//        /* init */ B1(0.5*(B1e_SI+B1s_SI)*cs_SI/(mu_SI*b_SI)),
//        /* init */ kB(kB_SI/mu_SI/std::pow(b_SI,3))
//        
//        {/*! Empty constructor is required by constexpr
//          */
//            
//            
//            
//        }
        
        /**********************************************************************/
        DislocationMobility(const Material<3,Isotropic>& material) :
        /* init */ DislocationMobilityBase("FCC DislocationMobility for "+material.materialName),
        /* init */ B0e(TextFileParser(material.materialFile).readScalar<double>("B0e_SI",true)*material.cs_SI/(material.mu_SI*material.b_SI)),
        /* init */ B1e(TextFileParser(material.materialFile).readScalar<double>("B1e_SI",true)*material.cs_SI/(material.mu_SI*material.b_SI)),
        /* init */ B0s(TextFileParser(material.materialFile).readScalar<double>("B0s_SI",true)*material.cs_SI/(material.mu_SI*material.b_SI)),
        /* init */ B1s(TextFileParser(material.materialFile).readScalar<double>("B1s_SI",true)*material.cs_SI/(material.mu_SI*material.b_SI)),
//        /* init */ kB(kB_SI/mu_SI/std::pow(b_SI,3))
        /* init */ kB(kB_SI/material.mu_SI/std::pow(material.b_SI,3))
        {/*! Empty constructor is required by constexpr
          */
            
            
            
        }
        
//        /**********************************************************************/
//        double velocity(const MatrixDim& S,
//                        const VectorDim& b,
//                        const VectorDim& , // xi
//                        const VectorDim& n,
//                        const double& T,
//                        const double& dL,
//                        const double& dt,
//                        const bool& use_stochasticForce) const override
//        {
//            double v=std::fabs(b.transpose()*S*n)/(B0+B1*T);
//            if(use_stochasticForce)
//            {
//                v+=StochasticForceGenerator::velocity(kB,T,B0+B1*T,dL,dt);
//            }
//            
//            // Interpolate ve and vs
//            const double cos2=std::pow(b.normalized().dot(xi),2);
//            const double sin2=1.0-cos2;
//            return vs*cos2+ve*sin2;
//
//        }

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
    
    /**************************************************************************/
    /**************************************************************************/
    template <>
    struct DislocationMobility<BCClattice<3>> : DislocationMobilityBase
    {
        
        typedef Eigen::Matrix<double,3,3> MatrixDim;
        typedef Eigen::Matrix<double,3,1> VectorDim;
        
        
        static constexpr double kB_eV=8.617e-5;       // Boltzmann constant in [eV]
        static constexpr double kB_SI=1.38064852e-23; // Boltzmann constant in SI units [J/K]
        
        const double h;
        const double w;
        const double B0e;
        const double B1e;
        const double B0s;
        const double B1s;
        const double Bk;
        const double dH0;
        const double p;
        const double q;
        const double T0;
        const double tauC;
        const double a0;
        const double a1;
        const double a2;
        const double a3;
        const double kB;
        
        /**********************************************************************/
        DislocationMobility(const double& b_SI,
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
                            const double& a3_in
        //                            const double& a4_in
        ) :
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
        //        /* init */ a4(a4_in),
        /* init */ kB(kB_SI/mu_SI/std::pow(b_SI,3))
        {/*! Empty constructor is required by constexpr
          */
            
        }

        
        /**********************************************************************/
        DislocationMobility(const Material<3,Isotropic>& material) :
                /* init */ DislocationMobilityBase("BCC DislocationMobility for "+material.materialName),
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
        static double sigmoid(const double & x)
        {
            return 2.0/(1.0+exp(2.0*x));
        }
        
        /**********************************************************************/
        double velocity(const MatrixDim& S,
                        const VectorDim& b,
                        const VectorDim& xi,
                        const VectorDim& n,
                        const double& T,
                        const double& dL,
                        const double& dt,
                        const bool& use_stochasticForce) const
        {
            
            const double bNorm=b.norm();
            const VectorDim s = b/bNorm;
            const VectorDim n1 = Eigen::AngleAxisd(M_PI/3.0,s)*n;
            
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
            
            if(use_stochasticForce)
            {
                vs+=StochasticForceGenerator::velocity(kB,T,Bs,dL,dt);
                
                ve+=StochasticForceGenerator::velocity(kB,T,B0e+B1e*T,dL,dt);
            }
            
            // Interpolate ve and vs
            const double cos2=std::pow(s.dot(xi),2);
//            const double sin2=1.0-cos2;
            return vs*cos2+ve*(1.0-cos2);
        }
        
    };
    
}
#endif



//        static double sigmoid(const double & x)
//        {
//            return 1.0/(1.0+exp(-x));
//        }
//
//        /**********************************************************************/
//        double velocity(const MatrixDim& S,
//                        const VectorDim& b,
//                        const VectorDim& xi,
//                        const VectorDim& n,
//                        const double& T,
//                        const double& dL,
//                        const double& dt,
//                        const bool& use_stochasticForce) const
//        {
//
//            const double bNorm=b.norm();
//            const VectorDim s = b/bNorm;
//            const VectorDim n1 = Eigen::AngleAxisd(M_PI/3.0,s)*n;
//
//            // Compute components of non-Schmid model
//            const double tau=std::fabs(s.transpose()*S*n); // magnitude of resolved shear stress
//            const double tauOrt=n.cross(s).transpose()*S*n;
//            const double tau1=std::fabs(s.transpose()*S*n1); // resolved schear stress on
//            const double tauOrt1=n1.cross(s).transpose()*S*n1;
//
//            const double num=tau+a1*tau1;
//            //
//            assert(num>=0.0 && "num must be >= 0.");
//            const double den=(1.0+0.5*a1)*tauC*(a4+a0*sigmoid(-(a2*tauOrt+a3*tauOrt1)/tauC));
//            assert(den>0.0 && "den must be > 0.");
//
//            const double Theta=num/den;
//            const double dg = (Theta<1.0)? (std::pow(1.0-std::pow(Theta,p),q)-T/T0) : 0.0;
//            const double dg1 = (dg>0.0)? dg : 0.0;
//            const double expCoeff = exp(-dH0*dg1/(2.0*kB_eV*T));
//
//            // Compute screw drag coeff
//            const double sgm=sigmoid((0.05-dg1)/0.05);
//            const double Bs=Bk*w/(2.0*h)*(1.0-sgm)+(B0s+B1s*T)*sgm; //kink-dominated to drag-dominated interpolation
//
//            // Compute screw velocity
//            double vs=tau*bNorm/Bs*expCoeff;
//
//            // Compute edge velocity
//            double ve=tau*bNorm/(B0e+B1e*T);
//
//            if(use_stochasticForce)
//            {
//                vs+=StochasticForceGenerator::velocity(kB,T,Bs,dL,dt);
//
//                ve+=StochasticForceGenerator::velocity(kB,T,B0e+B1e*T,dL,dt);
//            }
//
//            // Interpolate ve and vs
//            const double cos2=std::pow(b.normalized().dot(xi),2);
//            const double sin2=1.0-cos2;
//            return vs*cos2+ve*sin2;
//        }
