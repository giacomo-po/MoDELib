/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_COPPER_H_
#define model_COPPER_H_

#include <Eigen/Core>

#include <model/Dislocations/Materials/Crystal.h>
#include <model/Dislocations/Materials/IsotropicMaterial.h>

namespace model {
	
	class Copper :	public model::Crystal<3,FCC>,
	/*			*/	public model::IsotropicMaterial{
		
	public:	
		/////////////////////////////////////////////////////////////
		Copper() : IsotropicMaterial::IsotropicMaterial(0.34,48e9,0.2556e-9,1.0e-4,8940.0){	
		//			Copper() : IsotropicMaterial::IsotropicMaterial(0.375,48e9,0.2556e-9,1.0e-4,8940.0){	

			/*! Values for copper in SI units:
			* nu = 0.34;
			* mu = 48e9;
			* a = 0.36e-9;
			* b = a/sqrt(2)=0.2556e-9
			* B = 1.0e-4;
			* rho = 8940;
			*/
			
		}
		
	};
	
	//////////////////////////////////////////////////////////////
} // namespace model 
#endif

