/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2018 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2018 by Yinan Cui <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_EshelbyInclusion_H_
#define model_EshelbyInclusion_H_

#ifndef _MODEL_NON_SINGULAR_DD_
#define _MODEL_NON_SINGULAR_DD_ 1
#endif

#include <Eigen/Dense>
//#include <Material.h>
#include <StaticID.h>

//#include <DislocationStress.h>

#include <SlipSystem.h>

namespace model
{
    
    
    template <int dim>
    class EshelbyInclusion : public StaticID<EshelbyInclusion<dim>>
    {
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        typedef Eigen::Matrix<double,dim,1>   VectorDim;
        
    public:
        const double nu;
        const double mu;
        const double mobilityReduction;
        const int typeID;
        const double L;
        const double M;
        const double K;

    private:
        VectorDim C;
        double a;
        double a2;
        MatrixDim eT;
        double eTNorm;
        MatrixDim pT;
        
        
    public:
        
        static std::map<const GammaSurface*,GammaSurface> gammaSurfaceMap;
        
        static void addSlipSystems(const std::vector<std::shared_ptr<SlipSystem>>& slipSystems);
        double misfitEnergy(const Eigen::Matrix<double,3,1>& b, const GammaSurface* const matrixGammaSurface);
        /**********************************************************************/
        EshelbyInclusion(const VectorDim& _C,
                         const double& _a,
                         const MatrixDim& _eT,
                         const double& _nu,
                         const double& _mu,
                         const double& _mobilityReduction,
                         const int& _type) ;
        /**********************************************************************/
        bool contains(const VectorDim& x) const;
        /**********************************************************************/
        MatrixDim stress(const VectorDim& x) const;
    };
    
    
    template <int dim>
    std::map<const GammaSurface*,GammaSurface> EshelbyInclusion<dim>::gammaSurfaceMap;
}
#endif
