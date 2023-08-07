/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2018 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2018 by Yinan Cui <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_EshelbyInclusionBase_cpp_
#define model_EshelbyInclusionBase_cpp_


#include <EshelbyInclusionBase.h>

namespace model
{


//    template <int dim>
//    void EshelbyInclusionBase<dim>::addSlipSystems(const std::vector<std::shared_ptr<SlipSystem>>& slipSystems)
//    {
//
//        for(const auto& slipSystem : slipSystems)
//        {
//            if(gammaSurfaceMap.find(slipSystem->gammaSurface.get())==gammaSurfaceMap.end())
//            {// current slipSystem gammaSurface not found
//
//                //                    GammaSurface temp();
//                //
//                //                    gammaSurfaceMap.emaplace(slipSystem->gammaSurface.get(),temp);
//            }
//        }
//
//    }

//    template <int dim>
//    double EshelbyInclusionBase<dim>::misfitEnergy(const Eigen::Matrix<double,3,1>& b, const GammaSurface* const matrixGammaSurface) const
//    {
//        const auto iter(gammaSurfaceMap.find(matrixGammaSurface));
//        return iter==gammaSurfaceMap.end()? 0.0 : iter->second(b);
//    }

    /**********************************************************************/
    template <int dim>
    EshelbyInclusionBase<dim>::EshelbyInclusionBase(
                                            const MatrixDim& _eT,
                                            const double& _nu,
                                            const double& _mu,
                                            const double& _mobilityReduction,
                                            const int& _phaseID,
                                            const std::shared_ptr<SecondPhase<dim>>& sph) :
    /* init */ nu(_nu)
    /* init */,mu(_mu)
    /* init */,lambda(2.0*mu*nu/(1.0-2.0*nu))
    /* init */,mobilityReduction(_mobilityReduction)
    /* init */,phaseID(_phaseID)
    /* init */,secondPhase(sph)
    /* init */,eT(_eT)
    /* init */,eTNorm(eT.norm())
    /* init */,pT(2.0*mu*(eT+nu/(1.0-2.0*nu)*eT.trace()*MatrixDim::Identity()))
    {
        assert((_eT-_eT.transpose()).norm()<FLT_EPSILON && "eT is not symmetric.");
    }

    template <int dim>
    EshelbyInclusionBase<dim>::~EshelbyInclusionBase(){};
  
    template class EshelbyInclusionBase<3>;

}
#endif
