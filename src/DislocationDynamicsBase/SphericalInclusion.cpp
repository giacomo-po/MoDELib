/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2018 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2018 by Yinan Cui <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SphericalInclusion_cpp_
#define model_SphericalInclusion_cpp_


#include <SphericalInclusion.h>

namespace model
{


//    template <int dim>
//    void SphericalInclusion<dim>::addSlipSystems(const std::vector<std::shared_ptr<SlipSystem>>& slipSystems)
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
//    double SphericalInclusion<dim>::misfitEnergy(const Eigen::Matrix<double,3,1>& b, const GammaSurface* const matrixGammaSurface) const
//    {
//        const auto iter(gammaSurfaceMap.find(matrixGammaSurface));
//        return iter==gammaSurfaceMap.end()? 0.0 : iter->second(b);
//    }

    /**********************************************************************/
    template <int dim>
    SphericalInclusion<dim>::SphericalInclusion(const VectorDim& _C,
                                            const double& _a,
                                            const MatrixDim& _eT,
                                            const double& _nu,
                                            const double& _mu,
                                            const double& _mobilityReduction,
                                            const int& _phaseID,
                                            const std::shared_ptr<SecondPhase<dim>>& sph) :
    /* init */ EshelbyInclusionBase<dim>(_eT,_nu,_mu,_mobilityReduction,_phaseID,sph)
    /* init */,L((5.0*this->nu-1.0)/15.0/(1.0-this->nu))
    /* init */,M((4.0-5.0*this->nu)/15.0/(1.0-this->nu))
    //        /* init */,K((3.0*L+2.0*M)/3.0)
    /* init */,C(_C)
    /* init */,a(_a)
    /* init */,a2(a*a)
    {
//        std::cout<<"Creating SphericalInclusion "<<this->sID<<" (type "<<this->phaseID<<"):\n C="<<C.transpose()<<"\n a="<<a<<"\n eT="<<this->eT<<std::endl;
    }

    /**********************************************************************/
    template <int dim>
    bool SphericalInclusion<dim>::contains(const VectorDim& x) const
    {
        return (x-C).squaredNorm()<a2;
    }

    /**********************************************************************/
    template <int dim>
    typename SphericalInclusion<dim>::MatrixDim SphericalInclusion<dim>::stress(const VectorDim& x) const
    {
        
        if(this->eTNorm>FLT_EPSILON)
        {
            const VectorDim r(x-C);
            const double R2(r.squaredNorm());
            const double R(sqrt(R2));
            
            if(R>a)
            {
                const double R3=std::pow(R,3);
                const double R4=std::pow(R,4);
                const double a2R2=std::pow(a,2)/R2;
                const double a3R3=std::pow(a,3)/R3;
                
                const VectorDim pTr=this->pT*r;
                const double pTrr=pTr.dot(r);
                const double pTt=this->pT.trace();
                return a3R3/2.0/(1-this->nu)*( (10.0*(1.0-2.0*this->nu)+6.0*a2R2)/15.0*this->pT
                                        +(2.0*this->nu-2.0*a2R2)/R2*(pTr*r.transpose()+r*pTr.transpose())
                                        +((3.0*a2R2-5.0*(1.0-2.0*this->nu))/15.0*pTt + (1.0-2.0*this->nu-a2R2)/R2*pTrr)*MatrixDim::Identity()
                                        +(-(5.0-7.0*a2R2)/R4*pTrr+(1.0-a2R2)/R2*pTt)*r*r.transpose()
                                        );
            }
            else
            {
                return 2.0*this->mu*((L+this->nu/(1.0-2.0*this->nu)*(2.0*M+3.0*L))*this->eT.trace()*MatrixDim::Identity()+2.0*M*this->eT)-this->pT;
            }
        }
        else
        {
            return MatrixDim::Zero();
        }
    }

    template class SphericalInclusion<3>;

}
#endif
