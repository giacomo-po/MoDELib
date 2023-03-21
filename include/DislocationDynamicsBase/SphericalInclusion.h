/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2018 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2018 by Yinan Cui <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SphericalInclusion_H_
#define model_SphericalInclusion_H_


#include <Eigen/Dense>
//#include <Material.h>
#include <StaticID.h>

// http://solidmechanics.org/text/Chapter5_4/Chapter5_4.htm
#include <SlipSystem.h>
#include <SecondPhase.h>
#include <EshelbyInclusionBase.h>

namespace model
{
    
    
    template <int dim>
    class SphericalInclusion : public EshelbyInclusionBase<dim>
    {
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        typedef Eigen::Matrix<double,dim,1>   VectorDim;
        
    public:
//        const double nu;
//        const double mu;
//        const double mobilityReduction;
//        const int phaseID;
//        const std::shared_ptr<SecondPhase<dim>> secondPhase;
        const double L;
        const double M;
//        const double K;

        VectorDim C;
        double a;
        double a2;
//        MatrixDim eT;
//        double eTNorm;
//        MatrixDim pT;
        
        
    public:
        
//        static std::map<const GammaSurface*,GammaSurface> gammaSurfaceMap;
//
//        static void addSlipSystems(const std::vector<std::shared_ptr<SlipSystem>>& slipSystems);
        
        
//        double misfitEnergy(const Eigen::Matrix<double,3,1>& b, const GammaSurface* const matrixGammaSurface) const;

        
        /**********************************************************************/
        SphericalInclusion(const VectorDim& _C,
                         const double& _a,
                         const MatrixDim& _eT,
                         const double& _nu,
                         const double& _mu,
                         const double& _mobilityReduction,
                         const int& _phaseID,
                         const std::shared_ptr<SecondPhase<dim>>& sph);
  
        
        bool contains(const VectorDim& x) const override;
        MatrixDim stress(const VectorDim& x) const override;
        
    };
    
    
//    template <int dim>
//    std::map<const GammaSurface*,GammaSurface> SphericalInclusion<dim>::gammaSurfaceMap;
}
#endif
