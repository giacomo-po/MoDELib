/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_CRYSTALSTRUCTURES_H_
#define model_CRYSTALSTRUCTURES_H_

#include <vector>
#include <Eigen/Dense>

#include <model/DislocationDynamics/Materials/SlipSystem.h>

namespace model {
    
    /**************************************************************************/
    /**************************************************************************/
    struct FCC {
        
        template <int dim>
        static std::vector<Eigen::Matrix<double,dim,1> > getPlaneNormals(){
            typedef Eigen::Matrix<double,dim,1> VectorDim;
            typedef std::vector<VectorDim> PlaneNormalContainerType;
            
            VectorDim alpha(-1.0, 1.0,-1.0);
			VectorDim  beta( 1.0,-1.0,-1.0);
			VectorDim gamma(-1.0,-1.0, 1.0);
			VectorDim delta( 1.0, 1.0, 1.0);
            
            PlaneNormalContainerType temp;
            temp.push_back(alpha);
            temp.push_back(beta);
            temp.push_back(gamma);
            temp.push_back(delta);
            
            return temp;
        }
        
        /**********************************************************************/
        template <int dim>
        static std::vector<SlipSystem<dim> > getSlipSystems(const Eigen::Matrix<double,dim,dim>& C2G=Eigen::Matrix<double,dim,dim>::Identity()){
            typedef Eigen::Matrix<double,dim,1> VectorDim;
            typedef SlipSystem<dim> SlipSystemDim;
            typedef std::vector<SlipSystemDim> SlipSystemContainerType;
            
            
            SlipSystemContainerType temp;
            
            // The alpha-plane
            const VectorDim alpha(-1.0, 1.0,-1.0);
            const VectorDim alpha_BC( 1.0, 0.0,-1.0);
            const VectorDim alpha_CD(-1.0,-1.0, 0.0);
            const VectorDim alpha_DB( 0.0, 1.0, 1.0);
            temp.push_back(SlipSystemDim(C2G*alpha.normalized(),C2G*alpha_BC.normalized()));
            temp.push_back(SlipSystemDim(C2G*alpha.normalized(),C2G*alpha_CD.normalized()));
            temp.push_back(SlipSystemDim(C2G*alpha.normalized(),C2G*alpha_DB.normalized()));

            // The beta-plane
			const VectorDim beta( 1.0,-1.0,-1.0);
            const VectorDim beta_CA( 0.0,-1.0, 1.0);
            const VectorDim beta_AD(-1.0, 0.0,-1.0);
            const VectorDim beta_DC( 1.0, 1.0, 0.0);
            temp.push_back(SlipSystemDim(C2G*beta.normalized(),C2G*beta_CA.normalized()));
            temp.push_back(SlipSystemDim(C2G*beta.normalized(),C2G*beta_AD.normalized()));
            temp.push_back(SlipSystemDim(C2G*beta.normalized(),C2G*beta_DC.normalized()));

            
            // The gamma-plane
			const VectorDim gamma(-1.0,-1.0, 1.0);
            const VectorDim gamma_AB(-1.0, 1.0, 0.0);
            const VectorDim gamma_BD( 0.0,-1.0,-1.0);
            const VectorDim gamma_DA( 1.0, 0.0, 1.0);
            temp.push_back(SlipSystemDim(C2G*gamma.normalized(),C2G*gamma_AB.normalized()));
            temp.push_back(SlipSystemDim(C2G*gamma.normalized(),C2G*gamma_BD.normalized()));
            temp.push_back(SlipSystemDim(C2G*gamma.normalized(),C2G*gamma_DA.normalized()));
            
            // The delta-plane
			const VectorDim delta( 1.0, 1.0, 1.0);
            const VectorDim delta_AC( 0.0, 1.0,-1.0);
            const VectorDim delta_CB(-1.0, 0.0, 1.0);
            const VectorDim delta_BA( 1.0,-1.0, 0.0);
            temp.push_back(SlipSystemDim(C2G*delta.normalized(),C2G*delta_AC.normalized()));
            temp.push_back(SlipSystemDim(C2G*delta.normalized(),C2G*delta_CB.normalized()));
            temp.push_back(SlipSystemDim(C2G*delta.normalized(),C2G*delta_BA.normalized()));

            
            return temp;
        }

        
    };
    /**************************************************************************/
} // namespace model
#endif
