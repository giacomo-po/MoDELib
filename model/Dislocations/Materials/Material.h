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
//#include <model/Dislocations/Materials/Copper.h>
#include <model/Dislocations/Materials/Aluminum.h>
#include <model/Dislocations/Materials/Nickel.h>
#include <model/Dislocations/Materials/CrystalOrientation.h>


namespace model {

    
    
    
    /**************************************************************************/
    /**************************************************************************/
    template<typename SymmetryType>
    class Material { };
    
    /**************************************************************************/
    /**************************************************************************/
    template<>
    class Material<Isotropic> {
        
        
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
        static double C1;
        static double C2;
        static double C3;
        static double C4;
        
        /* select **************************************/
        template<template<typename T> class MaterialType>
        static void select(){
            typedef MaterialType<Isotropic> IM;
            mu=1.0;
            b =1.0;
            B =1.0;
            Binv=1.0/B;
            rho=IM::mu*std::pow(IM::b/IM::B,2)*IM::rho;  //! rho* = mu*(b/B)^2 * rho
            //        cs=IM::B/IM::b*std::pow(IM::rho*IM::mu,-0.5);
            cs=std::pow(mu/rho,0.5); // sc*=mu*/rho*
            nu=IM::nu;
            C1=1.0-nu;
            C2=1.0/(4.0*M_PI*C1);
            C3=1.0-2.0*nu;
            C4=0.5*C2;
                        
            std::string magentaColor    = "\033[0;35m";   // a magenta color
            std::string defaultColor    = "\033[0m";	   // the default color for the console
            std::cout<<magentaColor<<"Material is now: "<<IM::name<<defaultColor<<std::endl;

        }
        
        /* select **************************************/
        static void select(const unsigned int& Z){
        
            switch (Z) {
                case Al:
                    selectedMaterial=Al;
                    select<Aluminum>();
                    break;
                case Ni:
                    selectedMaterial=Ni;
                    select<Nickel>();
                    break;
                case Cu:
                    selectedMaterial=Cu;
                    select<Copper>();
                    break;
                    
                default:
                    assert(0 && "Material not implemented.");
                    break;
            }
                        
        }
        
        /* select **************************************/
        template <int dim>
        static void rotateCrystal(const Eigen::Matrix<double,dim,dim>& C2G){
            
            switch (selectedMaterial) {
                case Al:
                    CrystalOrientation<dim>::template rotate<Aluminum<Isotropic>::CrystalStructure>(C2G);
                    break;
                case Ni:
                    CrystalOrientation<dim>::template rotate<Nickel<Isotropic>::CrystalStructure>(C2G);
                    break;
                case Cu:
                    CrystalOrientation<dim>::template rotate<Copper<Isotropic>::CrystalStructure>(C2G);
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
    double Material<Isotropic>::C1=1.0-0.34;    // 1-nu
    double Material<Isotropic>::C2=1.0/(4.0*M_PI*(1.0-0.34));  // 1/(4*pi*(1-nu))
    double Material<Isotropic>::C3=1.0-2.0*0.34; // 1-2*nu
    double Material<Isotropic>::C4=1.0/(8.0*M_PI*(1.0-0.34));  // 1/(4*pi*(1-nu))
	
    /**************************************************************************/
    /**************************************************************************/
} // namespace model
#endif










//	class Material {
//
//
//	public:
//
//		const double nu;
//		const double mu;	// =1
//		const double b;		// =1
//		const double B;		// =1
//		const double Binv;		// =1
//
//		//! rho_star = mu*b^2/B^2 * rho
//		const double rho;
//
//		//! cs_star = B/b * sqrt(1/(rho*mu))
//		const double cs;
//
//		const double C1;
//		const double C2;
//		const double C3;
//		const double C4;
//
//
//		/////////////////////////////////////////////////////////////
//		Material(const double & nu_in, const double & mu_in,
//						  const double & b_in, const double & B_in,
//						  const double & rho_in) :	nu(nu_in),
//		/*										*/	mu(1.0),
//		/*										*/	b(1.0),
//		/*										*/	B(1.0),
//		/*										*/	Binv(1.0/B),
//		/*										*/	rho(mu_in*std::pow(b_in/B_in,2)*rho_in),
//		/*										*/	cs(B_in/b_in*std::pow(rho_in*mu_in,-0.5)),
//		/*										*/	C1(1.0-nu),
//		/*										*/	C2(mu/(4.0*M_PI*C1)),
//		/*										*/	C3(1.0-2.0*nu),
//		/*										*/	C4(1.0/(8.0*M_PI*C1)){
//
//			std::cout<<"DIMENSIONLESS MATERIAL PROPERTIES:"<<std::endl;
//			std::cout<<"	Poisson Ratio nu= "<<nu<<std::endl;
//			std::cout<<"	Shear Modulus mu= "<<mu<<std::endl;
//			std::cout<<"	Burgers Vector b=  "<<b<<std::endl;
//			std::cout<<"	Mobility B=  "<<B<<std::endl;
//			std::cout<<"	Density rho="<<rho<<std::endl;
//			std::cout<<"	Shear Wave Velocity cs="<<cs<<std::endl;
//
//		}
//
//	};

