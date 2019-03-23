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
#include <model/DislocationDynamics/Materials/Material.h>
#include <model/Utilities/StaticID.h>

//#include <model/DislocationDynamics/NearestNeighbor/DislocationStress.h>

namespace model
{
    
    
    template <int dim>
    class EshelbyInclusion : public StaticID<EshelbyInclusion<dim>>
    {
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        typedef Eigen::Matrix<double,dim,1>   VectorDim;
        
        
        VectorDim C;
        double a;
        MatrixDim eT;
        const double nu;
        const double mu;
        const int typeID;
        const double L;
        const double M;
        const double K;
        MatrixDim pT;
        
        
    public:
        
        /**********************************************************************/
        EshelbyInclusion(const VectorDim& _C,
                         const double& _a,
                         const MatrixDim& _eT,
                         const double& _nu,
                         const double& _mu,
                         const int& _type) :
        /* init */ C(_C)
        /* init */,a(_a)
        /* init */,eT(_eT)
        /* init */,nu(_nu)
        /* init */,mu(_mu)
        /* init */,typeID(_type)
        /* init */,L((5.0*nu-1.0)/15.0/(1.0-nu))
        /* init */,M((4.0-5.0*nu)/15.0/(1.0-nu))
        /* init */,K((3.0*L+2.0*M)/3.0)
        /* init */,pT(2.0*mu*(eT+nu/(1.0-2.0*nu)*eT.trace()*MatrixDim::Identity()))
        {
            model::cout<<"Creating EshelbyInclusion "<<this->sID<<" (type "<<typeID<<"):\n C="<<C.transpose()<<"\n a="<<a<<"\n eT="<<eT<<std::endl;
            assert((_eT-_eT.transpose()).norm()<FLT_EPSILON && "eT is not symmetric.");
        }
        
        /**********************************************************************/
        MatrixDim stress(const VectorDim& x) const
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
        
        
    };
    
}
#endif

