/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_ISOTROPICMATERIAL_H_
#define model_ISOTROPICMATERIAL_H_

#include <Eigen/Core>


namespace model {
	
	
	class IsotropicMaterial {
		
		
	public:
		
		const double nu;	
		const double mu;	// =1
		const double b;		// =1
		const double B;		// =1
		const double Binv;		// =1
		
		//! rho_star = mu*b^2/B^2 * rho 
		const double rho;
		
		//! cs_star = B/b * sqrt(1/(rho*mu))
		const double cs;	

		const double C1;
		const double C2;
		const double C3;
		const double C4;

		
		/////////////////////////////////////////////////////////////
		IsotropicMaterial(const double & nu_in, const double & mu_in, 
						  const double & b_in, const double & B_in,
						  const double & rho_in) :	nu(nu_in), 
		/*										*/	mu(1.0), 
		/*										*/	b(1.0), 
		/*										*/	B(1.0),
		/*										*/	Binv(1.0/B),
		/*										*/	rho(mu_in*std::pow(b_in/B_in,2)*rho_in),
		/*										*/	cs(B_in/b_in*std::pow(rho_in*mu_in,-0.5)),
		/*										*/	C1(1.0-nu),
		/*										*/	C2(mu/(4.0*M_PI*C1)),
		/*										*/	C3(1.0-2.0*nu),
		/*										*/	C4(1.0/(8.0*M_PI*C1)){
			
			std::cout<<"DIMENSIONLESS MATERIAL PROPERTIES:"<<std::endl;
			std::cout<<"	Poisson Ratio nu= "<<nu<<std::endl;
			std::cout<<"	Shear Modulus mu= "<<mu<<std::endl;
			std::cout<<"	Burgers Vector b=  "<<b<<std::endl;
			std::cout<<"	Mobility B=  "<<B<<std::endl;
			std::cout<<"	Density rho="<<rho<<std::endl;
			std::cout<<"	Shear Wave Velocity cs="<<cs<<std::endl;
			
		}
		
	};
	
	//////////////////////////////////////////////////////////////
} // namespace model 
#endif

