/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Tamer Crsoby     <tamercrosby@gmail.com>, 
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */



#ifndef model_STRAINRATECONTROLLER_H_
#define model_STRAINRATECONTROLLER_H_


#include <model/Dislocations/DislocationNetwork.h>


namespace model {
	
	struct StrainRateController {
		
		
		double zSurf;
		double targetStraiRate;
		
		
		//double zTractionInt; 
		//double kI;
		double kP;
		
		double zTraction;
		
		
	DislocationNetwork<3,1,CatmullRom,centripetal,16,GaussLegendre,Copper> DN;

		
		
		StrainRateController(){
		
			double disp(0.0);
		// Compute the actual dispacement at zSurface
			for (typename bvpfe::Domain::NodeContainerType::const_iterator iter=DN.shared.domain.nodeContainer.begin();iter!=DN.shared.domain.nodeContainer.end();++iter){
				if (!iter->isBoundaryNode) continue;
				if ( iter->P(2)!=zSurf) continue;
				//    std::cout<<"u = "<<iter->u.transpose()+iter->uInf.transpose()<<std::endl;
				Eigen::Matrix<double,dim,1> UT( iter->u.transpose() + iter->uInf.transpose() );
				counter++;
				disp +=  UT(2);
			}
			double actualStraiRate(disp/dt/zSurf);
			
			double straiRateError(targetStraiRate-actualStraiRate);
			
			//zTractionInt += kI*straiRateError*dt;

			
			double zTractionRate(kP*straiRateError/*+zTractionInt*/);
					
			zTraction += zTractionRate*dt;
			
			
			// call update b.c.
			
		
		}
		
		
		
	};
	
	//////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////
} // namespace model
#endif


