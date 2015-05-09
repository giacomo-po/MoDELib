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

#include <model/LatticeMath/LatticeMath.h>

namespace model
{
    
    /**************************************************************************/
    /**************************************************************************/
    struct FCC
    {
        
        template <int dim>
        static Eigen::Matrix<double,dim,dim> getLatticeBasis()
        {/*!\returns The matrix of lattice vectors (cartesian cooridinates in columns), 
          * in units of the crystallographic Burgers vector.
          */
            
            Eigen::Matrix<double,dim,dim> temp;
            temp << 0.0, 1.0, 1.0,
            /*   */ 1.0, 0.0, 1.0,
            /*   */ 1.0, 1.0, 0.0;
            
            return temp/sqrt(2.0);
        }
        
        template <int dim>
        static std::vector<LatticePlaneBase> reciprocalPlaneNormals()
        {/*!\returns a std::vector of ReciprocalLatticeDirection(s) corresponding
          * the slip plane normals of the FCC lattice
          */
            typedef Eigen::Matrix<long int,dim,1> VectorDimI;
            
            typedef LatticeVector<dim> LatticeVectorType;
            LatticeVectorType a1(VectorDimI(1,0,0));
            LatticeVectorType a2(VectorDimI(0,1,0));
            LatticeVectorType a3(VectorDimI(0,0,1));
            
//            ReciprocalLatticeDirectionType alpha(VectorDimI( 0, 0,-1));
//            ReciprocalLatticeDirectionType  beta(VectorDimI( 0,-1, 0));
//            ReciprocalLatticeDirectionType gamma(VectorDimI(-1, 0, 0)); // is (-1,-1, 1) in cartesian
//            ReciprocalLatticeDirectionType delta(VectorDimI( 1, 1, 1)); // is ( 1, 1, 1) in cartesian
            
            std::vector<LatticePlaneBase> temp;
            temp.emplace_back(a1,a3);           // is (-1, 1,-1) in cartesian
            temp.emplace_back(a3,a2);           // is ( 1,-1,-1) in cartesian
            temp.emplace_back(a2,a1);           // is (-1,-1, 1) in cartesian
            temp.emplace_back(a1-a3,a2-a3);     // is ( 1, 1, 1) in cartesian
            
            return temp;
        }


        
    };
    
    /**************************************************************************/
    /**************************************************************************/
    struct BCC
    {
        
        template <int dim>
        static Eigen::Matrix<double,dim,dim> getLatticeBasis()
        {/*!\returns The matrix of lattice vectors (in columns), in units of the
          * crystallographic Burgers vector.
          */
            
            Eigen::Matrix<double,dim,dim> temp;
            temp << -1.0,  1.0,  1.0,
            /*   */  1.0, -1.0,  1.0,
            /*   */  1.0,  1.0, -1.0;
            
            return temp/sqrt(3.0);
        }
        
        template <int dim>
        static std::vector<LatticePlaneBase> reciprocalPlaneNormals()
        {/*!\returns a std::vector of ReciprocalLatticeDirection(s) corresponding
          * the slip plane normals of the FCC lattice
          */
            typedef Eigen::Matrix<long int,dim,1> VectorDimI;
            
            typedef LatticeVector<dim> LatticeVectorType;
            LatticeVectorType a1(VectorDimI(1,0,0));
            LatticeVectorType a2(VectorDimI(0,1,0));
            LatticeVectorType a3(VectorDimI(0,0,1));
            LatticeVectorType  y(VectorDimI(1,1,1));
            
            std::vector<LatticePlaneBase> temp;
            temp.emplace_back(a3,a1); // is ( 1, 0, 1) in cartesian
            temp.emplace_back( y,a2); // is ( 1, 0,-1) in cartesian
            temp.emplace_back(a2,a3); // is ( 0, 1, 1) in cartesian
            temp.emplace_back( y,a1); // is ( 0,-1, 1) in cartesian
            temp.emplace_back(a1,a2); // is ( 1, 1, 0) in cartesian
            temp.emplace_back( y,a3); // is (-1, 1, 0) in cartesian
            
            return temp;
        }
        

    };    
    /**************************************************************************/

} // namespace model
#endif


//        template <int dim>
//        static std::vector<Eigen::Matrix<double,dim,1> > getPlaneNormals()
//        {
//            typedef Eigen::Matrix<double,dim,1> VectorDim;
//            typedef std::vector<VectorDim> PlaneNormalContainerType;
//            PlaneNormalContainerType temp;
//
//            VectorDim alpha(-1.0, 1.0,-1.0);
//			VectorDim  beta( 1.0,-1.0,-1.0);
//			VectorDim gamma(-1.0,-1.0, 1.0);
//			VectorDim delta( 1.0, 1.0, 1.0);
//
//            temp.push_back(alpha);
//            temp.push_back(beta);
//            temp.push_back(gamma);
//            temp.push_back(delta);
//
//            return temp;
//        }

//        template <int dim>
//        static std::vector<Eigen::Matrix<double,dim,1> > getPlaneNormals()
//        {
//            typedef Eigen::Matrix<double,dim,1> VectorDim;
//            typedef std::vector<VectorDim> PlaneNormalContainerType;
//            PlaneNormalContainerType temp;
//
//             // 6 {110} slip planes with 2 <111> slip directiosn = 12 slip systems
//            VectorDim system1_alpha(0.0,1.0,1.0);
//			VectorDim system1_beta( 1.0,0.0,1.0);
//			VectorDim system1_gamma(1.0,-1.0,0.0);
//			VectorDim system1_delta(0.0, 1.0, -1.0);
//            VectorDim system1_epsilon(1.0,1.0, 0.0);
//			VectorDim system1_zeta( 1.0,0.0, -1.0);
//
//            temp.push_back(system1_alpha);
//            temp.push_back(system1_beta);
//            temp.push_back(system1_gamma);
//            temp.push_back(system1_delta);
//            temp.push_back(system1_epsilon);
//            temp.push_back(system1_zeta);
//
//            // 12 {211} slip planes with 1 <111> slip directiosn = 12 slip systems
//            VectorDim system2_alpha(2.0,-1.0,1.0);
//			VectorDim system2_beta( 1.0,-2.0,-1.0);
//			VectorDim system2_gamma(1.0,1.0,2.0);
//			VectorDim system2_delta(2.0, 1.0, 1.0);
//            VectorDim system2_epsilon(1.0,2.0,-1.0);
//			VectorDim system2_zeta( 1.0, -1.0,2.0);
//            VectorDim system2_eta(2.0,1.0,-1.0);
//			VectorDim system2_iota( 1.0,2.0,1.0);
//			VectorDim system2_kappa(1.0,-1.0,-2.0);
//			VectorDim system2_lambda(2.0,-1.0,-1.0);
//            VectorDim system2_mu(1.0,-2.0,1.0);
//			VectorDim system2_nu( 1.0,1.0,-2.0);
//
//            temp.push_back(system2_alpha);
//            temp.push_back(system2_beta);
//            temp.push_back(system2_gamma);
//            temp.push_back(system2_delta);
//            temp.push_back(system2_epsilon);
//            temp.push_back(system2_zeta);
//            temp.push_back(system2_eta);
//            temp.push_back(system2_iota);
//            temp.push_back(system2_kappa);
//            temp.push_back(system2_lambda);
//            temp.push_back(system2_mu);
//            temp.push_back(system2_nu);
//
//
//            return temp;
//        }

//        /**************************************************************************/
//        template <int dim>
//        static std::vector<SlipSystem<dim> > getSlipSystems(const Eigen::Matrix<double,dim,dim>& C2G=Eigen::Matrix<double,dim,dim>::Identity())
//        {
//            typedef Eigen::Matrix<double,dim,1> VectorDim;
//            typedef SlipSystem<dim> SlipSystemDim;
//            typedef std::vector<SlipSystemDim> SlipSystemContainerType;
//
//
//            SlipSystemContainerType temp;
//                                        /********  {110} <111> ******/
//            // System 1 -> The alpha-plane
//            const VectorDim system1_alpha(0.0, 1.0,1.0);
//            const VectorDim system1_alpha_1(1.0, 1.0,-1.0);
//            const VectorDim system1_alpha_2(1.0,-1.0, 1.0);
//            temp.push_back(SlipSystemDim(C2G*system1_alpha.normalized(),C2G*system1_alpha_1.normalized()));
//            temp.push_back(SlipSystemDim(C2G*system1_alpha.normalized(),C2G*system1_alpha_2.normalized()));
//
//            // System 1 -> The beta-plane
//			const VectorDim system1_beta(1.0,0.0,1.0);
//            const VectorDim system1_beta_1( 1.0, 1.0,-1.0);
//            const VectorDim system1_beta_2(1.0,-1.0,-1.0);
//            temp.push_back(SlipSystemDim(C2G*system1_beta.normalized(),C2G*system1_beta_1.normalized()));
//            temp.push_back(SlipSystemDim(C2G*system1_beta.normalized(),C2G*system1_beta_2.normalized()));
//
//
//            // System 1 -> The gamma-plane
//			const VectorDim system1_gamma(1.0,-1.0,0.0);
//            const VectorDim system1_gamma_1( 1.0, 1.0,-1.0);
//            const VectorDim system1_gamma_2(1.0,1.0, 1.0);
//            temp.push_back(SlipSystemDim(C2G*system1_gamma.normalized(),C2G*system1_gamma_1.normalized()));
//            temp.push_back(SlipSystemDim(C2G*system1_gamma.normalized(),C2G*system1_gamma_2.normalized()));
//
//            // System 1 -> The delta-plane
//			const VectorDim system1_delta(0.0, 1.0,-1.0);
//            const VectorDim system1_delta_1( 1.0, -1.0,-1.0);
//            const VectorDim system1_delta_2(1.0,1.0, 1.0);
//            temp.push_back(SlipSystemDim(C2G*system1_delta.normalized(),C2G*system1_delta_1.normalized()));
//            temp.push_back(SlipSystemDim(C2G*system1_delta.normalized(),C2G*system1_delta_2.normalized()));
//
//            // System 1 -> The epsilon-plane
//			const VectorDim system1_epsilon(1.0,1.0,0.0);
//            const VectorDim system1_epsilon_1( 1.0,-1.0,-1.0);
//            const VectorDim system1_epsilon_2(1.0,-1.0,1.0);
//            temp.push_back(SlipSystemDim(C2G*system1_epsilon.normalized(),C2G*system1_epsilon_1.normalized()));
//            temp.push_back(SlipSystemDim(C2G*system1_epsilon.normalized(),C2G*system1_epsilon_2.normalized()));
//
//            // System 1 -> The zeta-plane
//			const VectorDim system1_zeta(1.0,0.0,-1.0);
//            const VectorDim system1_zeta_1( 1.0,-1.0,1.0);
//            const VectorDim system1_zeta_2(1.0,1.0, 1.0);
//            temp.push_back(SlipSystemDim(C2G*system1_zeta.normalized(),C2G*system1_zeta_1.normalized()));
//            temp.push_back(SlipSystemDim(C2G*system1_zeta.normalized(),C2G*system1_zeta_2.normalized()));
//
//                                        /********  {211} <111> ******/
//
//            // System 2 -> The alpha-plane
//			const VectorDim system2_alpha(2.0,-1.0,1.0);
//            const VectorDim system2_alpha_1(1.0,1.0,-1.0);
//            temp.push_back(SlipSystemDim(C2G*system2_alpha.normalized(),C2G*system2_alpha_1.normalized()));
//
//            // System 2 -> The beta-plane
//			const VectorDim system2_beta(1.0,-2.0,-1.0);
//            const VectorDim system2_beta_1(1.0,1.0,-1.0);
//            temp.push_back(SlipSystemDim(C2G*system2_beta.normalized(),C2G*system2_beta_1.normalized()));
//
//            // System 2 -> The gamma-plane
//			const VectorDim system2_gamma(1.0,1.0,2.0);
//            const VectorDim system2_gamma_1(1.0,1.0,-1.0);
//            temp.push_back(SlipSystemDim(C2G*system2_gamma.normalized(),C2G*system2_gamma_1.normalized()));
//
//            // System 2 -> The delta-plane
//			const VectorDim system2_delta(2.0,1.0,1.0);
//            const VectorDim system2_delta_1(1.0,-1.0,-1.0);
//            temp.push_back(SlipSystemDim(C2G*system2_delta.normalized(),C2G*system2_delta_1.normalized()));
//
//            // System 2 -> The epsilon-plane
//			const VectorDim system2_epsilon(1.0,2.0,-1.0);
//            const VectorDim system2_epsilon_1(1.0,-1.0,-1.0);
//            temp.push_back(SlipSystemDim(C2G*system2_epsilon.normalized(),C2G*system2_epsilon_1.normalized()));
//
//            // System 2 -> The zeta-plane
//			const VectorDim system2_zeta(1.0,-1.0,2.0);
//            const VectorDim system2_zeta_1(1.0,-1.0,-1.0);
//            temp.push_back(SlipSystemDim(C2G*system2_zeta.normalized(),C2G*system2_zeta_1.normalized()));
//
//            // System 2 -> The eta-plane
//			const VectorDim system2_eta(2.0,1.0,-1.0);
//            const VectorDim system2_eta_1(1.0,-1.0,1.0);
//            temp.push_back(SlipSystemDim(C2G*system2_eta.normalized(),C2G*system2_eta_1.normalized()));
//
//            // System 2 -> The iota-plane
//			const VectorDim system2_iota(1.0,2.0,1.0);
//            const VectorDim system2_iota_1(1.0,-1.0,1.0);
//            temp.push_back(SlipSystemDim(C2G*system2_iota.normalized(),C2G*system2_iota_1.normalized()));
//
//            // System 2 -> The kappa-plane
//			const VectorDim system2_kappa(1.0,-1.0,-2.0);
//            const VectorDim system2_kappa_1(1.0,-1.0,1.0);
//            temp.push_back(SlipSystemDim(C2G*system2_kappa.normalized(),C2G*system2_kappa_1.normalized()));
//
//            // System 2 -> The lambda-plane
//			const VectorDim system2_lambda(2.0,-1.0,-1.0);
//            const VectorDim system2_lambda_1(1.0,1.0,1.0);
//            temp.push_back(SlipSystemDim(C2G*system2_lambda.normalized(),C2G*system2_lambda_1.normalized()));
//
//            // System 2 -> The mu-plane
//			const VectorDim system2_mu(1.0,-2.0,1.0);
//            const VectorDim system2_mu_1(1.0,1.0,1.0);
//            temp.push_back(SlipSystemDim(C2G*system2_mu.normalized(),C2G*system2_mu_1.normalized()));
//
//            // System 2 -> The nu-plane
//			const VectorDim system2_nu(1.0,1.0,-2.0);
//            const VectorDim system2_nu_1(1.0,1.0,1.0);
//            temp.push_back(SlipSystemDim(C2G*system2_nu.normalized(),C2G*system2_nu_1.normalized()));
//
//
//            return temp;
//        }


//        /**********************************************************************/
//        template <int dim>
//        static std::vector<SlipSystem<dim> > getSlipSystems(const Eigen::Matrix<double,dim,dim>& C2G=Eigen::Matrix<double,dim,dim>::Identity())
//        {
//            typedef Eigen::Matrix<double,dim,1> VectorDim;
//            typedef SlipSystem<dim> SlipSystemDim;
//            typedef std::vector<SlipSystemDim> SlipSystemContainerType;
//
//
//            SlipSystemContainerType temp;
//
//            // The alpha-plane
//            const VectorDim alpha(-1.0, 1.0,-1.0);
//            const VectorDim alpha_BC( 1.0, 0.0,-1.0);
//            const VectorDim alpha_CD(-1.0,-1.0, 0.0);
//            const VectorDim alpha_DB( 0.0, 1.0, 1.0);
//            temp.push_back(SlipSystemDim(C2G*alpha.normalized(),C2G*alpha_BC.normalized()));
//            temp.push_back(SlipSystemDim(C2G*alpha.normalized(),C2G*alpha_CD.normalized()));
//            temp.push_back(SlipSystemDim(C2G*alpha.normalized(),C2G*alpha_DB.normalized()));
//
//            // The beta-plane
//			const VectorDim beta( 1.0,-1.0,-1.0);
//            const VectorDim beta_CA( 0.0,-1.0, 1.0);
//            const VectorDim beta_AD(-1.0, 0.0,-1.0);
//            const VectorDim beta_DC( 1.0, 1.0, 0.0);
//            temp.push_back(SlipSystemDim(C2G*beta.normalized(),C2G*beta_CA.normalized()));
//            temp.push_back(SlipSystemDim(C2G*beta.normalized(),C2G*beta_AD.normalized()));
//            temp.push_back(SlipSystemDim(C2G*beta.normalized(),C2G*beta_DC.normalized()));
//
//
//            // The gamma-plane
//			const VectorDim gamma(-1.0,-1.0, 1.0);
//            const VectorDim gamma_AB(-1.0, 1.0, 0.0);
//            const VectorDim gamma_BD( 0.0,-1.0,-1.0);
//            const VectorDim gamma_DA( 1.0, 0.0, 1.0);
//            temp.push_back(SlipSystemDim(C2G*gamma.normalized(),C2G*gamma_AB.normalized()));
//            temp.push_back(SlipSystemDim(C2G*gamma.normalized(),C2G*gamma_BD.normalized()));
//            temp.push_back(SlipSystemDim(C2G*gamma.normalized(),C2G*gamma_DA.normalized()));
//
//            // The delta-plane
//			const VectorDim delta( 1.0, 1.0, 1.0);
//            const VectorDim delta_AC( 0.0, 1.0,-1.0);
//            const VectorDim delta_CB(-1.0, 0.0, 1.0);
//            const VectorDim delta_BA( 1.0,-1.0, 0.0);
//            temp.push_back(SlipSystemDim(C2G*delta.normalized(),C2G*delta_AC.normalized()));
//            temp.push_back(SlipSystemDim(C2G*delta.normalized(),C2G*delta_CB.normalized()));
//            temp.push_back(SlipSystemDim(C2G*delta.normalized(),C2G*delta_BA.normalized()));
//
//
//            return temp;
//        }
