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
//#include <model/DislocationDynamics/Materials/Copper.h>
//#include <model/DislocationDynamics/Materials/Aluminum.h>
//#include <model/DislocationDynamics/Materials/Nickel.h>
#include <model/DislocationDynamics/Materials/PeriodicElement.h>
#include <model/DislocationDynamics/Materials/CrystalOrientation.h>



namespace model {

    
    
    
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
        
    public:
        
        enum{Al=13, Ni=28, Cu=29};
        
        static double mu;
        static double b;
        static double B;
        static double Binv;
        static double rho;
        static double cs;
        static double nu;
        static double lambda;
        static double C1;
        static double C2;
        static double C3;
        static double C4;
        static double T; // temperature in [K]
        static double kT; // temperature in [K]
        static double tauIII; // stage III tau
        static double vAct; // activation volume for cross slip
        
        /* select **************************************/
        //template<template<typename T> class MaterialType>
        template<int Z>
        static void select()
        {
        //    typedef MaterialType<Isotropic> IM;
                typedef PeriodicElement<Z,Isotropic> IM;

            mu=1.0;
            b =1.0;
            B =1.0;
            Binv=1.0/B;
            rho=IM::mu*std::pow(IM::b/IM::B,2)*IM::rho;  //! rho* = mu*(b/B)^2 * rho
            //        cs=IM::B/IM::b*std::pow(IM::rho*IM::mu,-0.5);
            cs=std::pow(mu/rho,0.5); // sc*=sqrt(mu*/rho*)
            nu=IM::nu;
            lambda=2.0*mu*nu/(1.0-2.0*nu);
            C1=1.0-nu;
            C2=1.0/(4.0*M_PI*C1);
            C3=1.0-2.0*nu;
            C4=0.5*C2;
            kT=1.38e-23/IM::mu/std::pow(IM::b,3)*T;
            
            std::string magentaColor    = "\033[0;35m";   // a magenta color
            std::string defaultColor    = "\033[0m";	   // the default color for the console
            std::cout<<magentaColor<<"Material is now: "<<IM::name<<defaultColor<<std::endl;

        }
        
        /* select **************************************/
        static void select(const unsigned int& Z)
        {
        
            switch (Z) {
                case Al:
                    selectedMaterial=Al;
                    select<Al>();
                    break;
                case Ni:
                    selectedMaterial=Ni;
                    select<Ni>();
                    break;
                case Cu:
                    selectedMaterial=Cu;
                    select<Cu>();
                    break;
                    
                default:
                    assert(0 && "Material not implemented.");
                    break;
            }
                        
        }
        
        /* select **************************************/
        template <int dim>
        static void rotateCrystal(const Eigen::Matrix<double,dim,dim>& C2G)
        {
            
            switch (selectedMaterial) {
                case Al:
//                    CrystalOrientation<dim>::template rotate<Aluminum<Isotropic>::CrystalStructure>(C2G);
                    CrystalOrientation<dim>::template rotate<PeriodicElement<Al,Isotropic>::CrystalStructure>(C2G);
                    break;
                case Ni:
//                    CrystalOrientation<dim>::template rotate<Nickel<Isotropic>::CrystalStructure>(C2G);
                    CrystalOrientation<dim>::template rotate<PeriodicElement<Ni,Isotropic>::CrystalStructure>(C2G);
                    break;
                case Cu:
//                    CrystalOrientation<dim>::template rotate<Copper<Isotropic>::CrystalStructure>(C2G);
                    CrystalOrientation<dim>::template rotate<PeriodicElement<Cu,Isotropic>::CrystalStructure>(C2G);
                    break;
                default:
                    assert(0 && "Material not implemented.");
                    break;
            }
        }
        
    };
    
    
    //double Material<Isotropic>::a2=1.0;  // square of core size a
    int Material<Isotropic>::selectedMaterial=29;
    double Material<Isotropic>::mu=1.0;
    double Material<Isotropic>::b=1.0;
    double Material<Isotropic>::B=1.0;
    double Material<Isotropic>::Binv=1.0;
    double Material<Isotropic>::rho=1.0;
    double Material<Isotropic>::cs=1.0;
    double Material<Isotropic>::nu=0.34;
    double Material<Isotropic>::lambda=2.0*1.0*0.34/(1.0-2.0*0.34);
    double Material<Isotropic>::C1=1.0-0.34;    // 1-nu
    double Material<Isotropic>::C2=1.0/(4.0*M_PI*(1.0-0.34));  // 1/(4*pi*(1-nu))
    double Material<Isotropic>::C3=1.0-2.0*0.34; // 1-2*nu
    double Material<Isotropic>::C4=1.0/(8.0*M_PI*(1.0-0.34));  // 1/(4*pi*(1-nu))
    double Material<Isotropic>::T=300.0;  // temperature in [K]
    double Material<Isotropic>::kT=1.38e-23/48.0e9/std::pow(0.2556e-9,3)*300.0;  // 
	double Material<Isotropic>::tauIII=0.667e-3;  // critical resolved shear stress in stage III normalized by mu
    double Material<Isotropic>::vAct=300.0;  // activation volume [b^3]
    
    /**************************************************************************/
    /**************************************************************************/
} // namespace model
#endif
