/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2018 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2018 by Yinan Cui <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_EshelbyInclusion_cpp_
#define model_EshelbyInclusion_cpp_

#ifndef _MODEL_NON_SINGULAR_DD_
#define _MODEL_NON_SINGULAR_DD_ 1
#endif

#include <EshelbyInclusion.h>

namespace model
{
    
    
        template <int dim>
        void EshelbyInclusion<dim>::addSlipSystems(const std::vector<std::shared_ptr<SlipSystem>>& slipSystems)
        {
            
            for(const auto& slipSystem : slipSystems)
            {
                if(gammaSurfaceMap.find(slipSystem->gammaSurface.get())==gammaSurfaceMap.end())
                {// current slipSystem gammaSurface not found
                    
//                    GammaSurface temp();
//
//                    gammaSurfaceMap.emaplace(slipSystem->gammaSurface.get(),temp);
                }
            }
            
        }
        
        template <int dim>
        double EshelbyInclusion<dim>::misfitEnergy(const Eigen::Matrix<double,3,1>& b, const GammaSurface* const matrixGammaSurface)
        {
            const auto iter(gammaSurfaceMap.find(matrixGammaSurface));
            return iter==gammaSurfaceMap.end()? 0.0 : iter->second(b);
        }
        
        /**********************************************************************/
        template <int dim>
        EshelbyInclusion<dim>::EshelbyInclusion(const VectorDim& _C,
                         const double& _a,
                         const MatrixDim& _eT,
                         const double& _nu,
                         const double& _mu,
                         const double& _mobilityReduction,
                         const int& _type) :
        /* init */ nu(_nu)
        /* init */,mu(_mu)
        /* init */,mobilityReduction(_mobilityReduction)
        /* init */,typeID(_type)
        /* init */,L((5.0*nu-1.0)/15.0/(1.0-nu))
        /* init */,M((4.0-5.0*nu)/15.0/(1.0-nu))
        /* init */,K((3.0*L+2.0*M)/3.0)
        /* init */,C(_C)
        /* init */,a(_a)
        /* init */,a2(a*a)
        /* init */,eT(_eT)
        /* init */,eTNorm(eT.norm())
        /* init */,pT(2.0*mu*(eT+nu/(1.0-2.0*nu)*eT.trace()*MatrixDim::Identity()))
        {
            std::cout<<"Creating EshelbyInclusion "<<this->sID<<" (type "<<typeID<<"):\n C="<<C.transpose()<<"\n a="<<a<<"\n eT="<<eT<<std::endl;
            assert((_eT-_eT.transpose()).norm()<FLT_EPSILON && "eT is not symmetric.");
        }
        
        /**********************************************************************/
        template <int dim>
        bool EshelbyInclusion<dim>::contains(const VectorDim& x) const
        {
            return (x-C).squaredNorm()<a2;
        }
        
        /**********************************************************************/
        template <int dim>
        typename EshelbyInclusion<dim>::MatrixDim EshelbyInclusion<dim>::stress(const VectorDim& x) const
        {
            
            const VectorDim r(x-C);
            const double R2(r.squaredNorm());
            const double R(sqrt(R2));
            
            if(eTNorm>FLT_EPSILON)
            {
                if(R>a)
                {
                    const double R3=std::pow(R,3);
                    const double R4=std::pow(R,4);
                    const double a2R2=std::pow(a,2)/R2;
                    const double a3R3=std::pow(a,3)/R3;
                    
                    const VectorDim pTr=pT*r;
                    const double pTrr=pTr.dot(r);
                    const double pTt=pT.trace();
                    return a3R3/2.0/(1-nu)*( (10.0*(1.0-2.0*nu)+6.0*a2R2)/15.0*pT
                                            +(2.0*nu-2.0*a2R2)/R2*(pTr*r.transpose()+r*pTr.transpose())
                                            +((3.0*a2R2-5.0*(1.0-2.0*nu))/15.0*pTt + (1.0-2.0*nu-a2R2)/R2*pTrr)*MatrixDim::Identity()
                                            +(-(5.0-7.0*a2R2)/R4*pTrr+(1.0-a2R2)/R2*pTt)*r*r.transpose()
                                            );
                }
                else
                {
                    return 2.0*mu*((L+nu/(1.0-2.0*nu)*3.0*K)*eT.trace()*MatrixDim::Identity()+2.0*M*eT);
                }
            }
            else
            {
                return MatrixDim::Zero();
            }
        }
        template class EshelbyInclusion<3>;

}
#endif
