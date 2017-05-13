/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_MATERIAL_H_
#define model_MATERIAL_H_

#include <cmath>
#include <model/DislocationDynamics/Materials/PeriodicElement.h>
#include <model/DislocationDynamics/Materials/CrystalOrientation.h>
#include <model/MPI/MPIcout.h> // defines mode::cout
#include <model/Utilities/TerminalColors.h> // defines mode::cout



namespace model
{
    
    /**************************************************************************/
    /**************************************************************************/
    template<typename SymmetryType>
    class Material { };
    
    /**************************************************************************/
    /**************************************************************************/
    template<>
    class Material<Isotropic>
    {
        
        
        static int selectedMaterial;
        
        static constexpr PeriodicElement<12,Isotropic> Mg=PeriodicElement<12,Isotropic>();
        static constexpr PeriodicElement<13,Isotropic> Al=PeriodicElement<13,Isotropic>();
        static constexpr PeriodicElement<28,Isotropic> Ni=PeriodicElement<28,Isotropic>();
        static constexpr PeriodicElement<29,Isotropic> Cu=PeriodicElement<29,Isotropic>();
        static constexpr PeriodicElement<26,Isotropic>  Fe=PeriodicElement<26,Isotropic>();
        static constexpr PeriodicElement<74,Isotropic>  W=PeriodicElement<74,Isotropic>();
        
        
        /**********************************************************************/
        template<int Z>
        static void select()
        {
            typedef PeriodicElement<Z,Isotropic> IM;
            
            mu=1.0;     // mu is used to normalized stress
            b =1.0;     // b is used to normalized length
            cs=1.0;     // sc*=sqrt(mu*/rho*)
            
            b_real=PeriodicElement<Z,Isotropic>::b;
            //            B =1.0;     // B=A*T [Pa*sec] is effectively used to normalize time
            //            Binv=1.0/B;
            //            rho=IM::mu*std::pow(IM::b/(IM::Ae*T),2)*IM::rho;  //! rho* = mu*(b/B)^2 * rho, B=A*T
            //            //        cs=IM::B/IM::b*std::pow(IM::rho*IM::mu,-0.5);
            //            cs=sqrt(mu/rho); // sc*=sqrt(mu*/rho*)
            ////            cs =1.0;    // shear wave speed is used to normalized velocity
            nu=IM::nu;
            lambda=2.0*mu*nu/(1.0-2.0*nu);
            C1=1.0-nu;
            C2=1.0/(4.0*M_PI*C1);
            C3=1.0-2.0*nu;
            C4=0.5*C2;
            
            /* still unused, to be used in BCC mobility */
            kB=1.38e-23/IM::mu/std::pow(IM::b,3); // [-]
            //            dH0=PeriodicElement<Z,Isotropic>::dH0/IM::mu/std::pow(IM::b,3); // [-]
            //            p=PeriodicElement<Z,Isotropic>::p; // [-]
            //            q=PeriodicElement<Z,Isotropic>::q; // [-]
            ////            Ae=PeriodicElement<Z,Isotropic>::Ae / IM::mu / (IM::b/sqrt(IM::mu/IM::rho));   // [1/K]
            //            Be=1.0;   // [1/K]
            //            Bs=PeriodicElement<Z,Isotropic>::As / PeriodicElement<Z,Isotropic>::Ae;   // [1/K]
            //            tauP=PeriodicElement<Z,Isotropic>::tauP / IM::mu ;  // [-]
            //            Ta=PeriodicElement<Z,Isotropic>::Ta;  // [K]
            
            model::cout<<magentaColor<<"Material is now: "<<IM::name<<defaultColor<<std::endl;
            model::cout<<greenColor<<"  units of stress (shear modulus): mu="<<IM::mu<<" [Pa]"<<std::endl;
            model::cout<<greenColor<<"  units of length (Burgers vector): b="<<IM::b<<" [m]"<<std::endl;
            model::cout<<greenColor<<"  units of speed (shear-wave speed): cs="<<IM::cs<<" [m/s]"<<std::endl;
            model::cout<<greenColor<<"  units of time: b/cs="<<IM::b/IM::cs<<" [sec]"<<defaultColor<<std::endl;
        }
        
        
    public:
        
        //        enum{
        //            _Al=13,
        //            _Fe=26,
        //            _Ni=28,
        //            _Cu=29,
        //            _W=74
        //        };
        
        static double b;
        static double b_real;
        static double mu;
        static double cs;
        
        static double T;
        //        static double rho;
        static double nu;
        static double lambda;
        static double C1;
        static double C2;
        static double C3;
        static double C4;
        static double kB;
        
        static Eigen::Matrix<double,Eigen::Dynamic,2> dH0;
        
        /**********************************************************************/
        static void select(const unsigned int& Z)
        {/*!\param[in] Z the atomic number of the element to be selected
          *
          * Selects the element Z as the current material. All material properties
          * are updated accordingly.
          */
            switch (Z)
            {
                case Mg.Z:
                    selectedMaterial=Mg.Z;
                    select<Mg.Z>();
                    break;
                case Al.Z:
                    selectedMaterial=Al.Z;
                    select<Al.Z>();
                    break;
                case Ni.Z:
                    selectedMaterial=Ni.Z;
                    select<Ni.Z>();
                    break;
                case Cu.Z:
                    selectedMaterial=Cu.Z;
                    select<Cu.Z>();
                    break;
                case Fe.Z:
                    selectedMaterial=Fe.Z;
                    select<Fe.Z>();
                    break;
                case W.Z:
                    selectedMaterial=W.Z;
                    select<W.Z>();
                    break;
                default:
                    assert(0 && "Material not implemented.");
                    break;
            }
            
        }
        
        /**********************************************************************/
        template <int dim>
        static void rotateCrystal(const Eigen::Matrix<double,dim,dim>& C2G)
        {/*!\param[in] C2G the crystal-to-global rotation matrix
          *
          * Rotates the default plane normals in
          * PeriodicElement<Z,Isotropic>::CrystalStructure, where Z is the current
          * selected material.
          */
            
            switch (selectedMaterial)
            {
                case Mg.Z:
                    CrystalOrientation<dim>::template rotate<PeriodicElement<Mg.Z,Isotropic>>(C2G);
                    break;
                case Al.Z:
                    CrystalOrientation<dim>::template rotate<PeriodicElement<Al.Z,Isotropic>>(C2G);
                    break;
                case Ni.Z:
                    CrystalOrientation<dim>::template rotate<PeriodicElement<Ni.Z,Isotropic>>(C2G);
                    break;
                case Cu.Z:
                    CrystalOrientation<dim>::template rotate<PeriodicElement<Cu.Z,Isotropic>>(C2G);
                    break;
                case Fe.Z:
                    CrystalOrientation<dim>::template rotate<PeriodicElement<Fe.Z,Isotropic>>(C2G);
                    break;
                case W.Z:
                    CrystalOrientation<dim>::template rotate<PeriodicElement<W.Z,Isotropic>>(C2G);
                    break;
                    
                default:
                    assert(0 && "Material not implemented.");
                    break;
            }
        }
        
        /**********************************************************************/
        static double velocity(const Eigen::Matrix<double,3,3>& S,
                               const Eigen::Matrix<double,3,1>& b,
                               const Eigen::Matrix<double,3,1>& xi,
                               const Eigen::Matrix<double,3,1>& n)
        {
            switch (selectedMaterial)
            {
                case Mg.Z:
                    return Mg.dm.velocity(S,b,xi,n,T);
                case Al.Z:
                    return Al.dm.velocity(S,b,xi,n,T);
                case Ni.Z:
                    return Ni.dm.velocity(S,b,xi,n,T);
                case Cu.Z:
                    return Cu.dm.velocity(S,b,xi,n,T);
                case Fe.Z:
                    return Fe.dm.velocity(S,b,xi,n,T);
                case W.Z:
                    return W.dm.velocity(S,b,xi,n,T);
                    
                default:
                    assert(0 && "velocity function not implemented.");
                    return 0.0;
            }
        }
        
    };
    
    // Static data
    int Material<Isotropic>::selectedMaterial=29;
    double Material<Isotropic>::b=1.0;    // dimensionless Burgers vector
    double Material<Isotropic>::b_real=0.2556e-9;    // dimensionless Burgers vector
    double Material<Isotropic>::mu=1.0;
    double Material<Isotropic>::cs=1.0;
    double Material<Isotropic>::T=300.0;  // Temperature [K]
    double Material<Isotropic>::nu=0.34;
    double Material<Isotropic>::lambda=2.0*1.0*0.34/(1.0-2.0*0.34);
    double Material<Isotropic>::C1=1.0-0.34;    // 1-nu
    double Material<Isotropic>::C2=1.0/(4.0*M_PI*(1.0-0.34));  // 1/(4*pi*(1-nu))
    double Material<Isotropic>::C3=1.0-2.0*0.34; // 1-2*nu
    double Material<Isotropic>::C4=1.0/(8.0*M_PI*(1.0-0.34));  // 1/(4*pi*(1-nu))
    double Material<Isotropic>::kB=1.38e-23/48.0e9/pow(0.2556e-9,3);  // Boltzmann constant
    
    /**************************************************************************/
    /**************************************************************************/
} // namespace model
#endif


//        /**********************************************************************/
//        static double velocity(const Eigen::Matrix<double,3,3>& S,
//                               const Eigen::Matrix<double,3,1>& b,
//                               const Eigen::Matrix<double,3,1>& xi,
//                               const Eigen::Matrix<double,3,1>& n)
//        {
//            switch (selectedMaterial)
//            {
//                case Al:
//                    return PeriodicElement<Al,Isotropic>::dm.velocity(S,b,xi,n,T);
//                case Ni:
//                    return PeriodicElement<Ni,Isotropic>::dm.velocity(S,b,xi,n,T);
//                case Cu:
//                    return PeriodicElement<Cu,Isotropic>::dm.velocity(S,b,xi,n,T);
//                    //              case Fe:
//                    //                  return PeriodicElement<Cu,Isotropic>::dm.velocity(S,b,xi,n,T);
//                case W:
//                    return PeriodicElement<W,Isotropic>::dm.velocity(S,b,xi,n,T);
//
//                default:
//                    assert(0 && "velocity function not implemented.");
//                    return 0.0;
//            }
//        }
