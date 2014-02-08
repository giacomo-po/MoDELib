/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _model_DislocationEnergy_h_
#define _model_DislocationEnergy_h_

#include <model/ParticleInteraction/FieldBase.h>
#include <model/ParticleInteraction/FieldPoint.h>
#include <model/DislocationDynamics/NearestNeighbor/DislocationStress.h>

namespace model
{
	
	/**************************************************************************/
	/**************************************************************************/
	template<short unsigned int _dim>
	class DislocationEnergy
    /* inheritance */ : public FieldBase<double,1,1>
    {
        
        
    public:
        
        
        typedef DislocationEnergy<_dim> DislocationEnergyType;
        typedef FieldBase<double,1,1> FieldBaseType;
        typedef typename FieldBaseType::MatrixType MatrixType;
        
        
        template <typename DislocationParticleType>
        static MatrixType compute(const DislocationParticleType& source,const DislocationParticleType& field)
        {/*!@param[in] source the DislocationParticle that is source of stress
          * @param[in] field  the DislocationParticle on which stress is computed
          *\returns the stress field produced by source on field
          */
            
            MatrixType temp(MatrixType::Zero());
            
            const Eigen::Matrix<double,_dim,1> r(field.P-source.P);
            const double R2(r.squaredNorm());
            if (R2>0.0)
            {
                const double R(sqrt(R2));
                const double RL(R/DislocationStress<_dim>::a);
                const double eRL(exp(-RL));
                const double lapA(2.0/R*(1.0-eRL)); // laplacian of A
                const double dAR((1.0 + 2.0*DislocationStress<_dim>::a2/R2*((1.0+RL)*eRL-1.0) )/R); // this is 1/R*dA/dR
                
                temp(0,0) = lapA*( 0.5*Material<Isotropic>::C1*source.B.dot(source.T)*field.B.dot( field.T)
                                  /*             */ +    Material<Isotropic>::nu*source.B.dot( field.T)*field.B.dot(source.T)
                                  /*             */ -                            source.B.dot( field.B)*field.T.dot(source.T)
                                  /*             */)
                /*      */ +field.T.dot(source.T)*( (lapA-3.0*dAR)*source.B.dot(r)*field.B.dot(r)/R2
                                                   /*                               */ +dAR*source.B.dot(field.B)
                                                   /*                               */ );
            }
            
            
            //            return  (Material<Isotropic>::C1*(1.0+0.5*coreLsquared/RaSquared)*Burgers.dot(rugauss.col(k))*bf.dot(ruf)
            //					 +2.0*Material<Isotropic>::nu*(1.0+0.5*coreLsquared/RaSquared)*(bf.dot(rugauss.col(k))*Burgers.dot(ruf))
            //					 -(Burgers.dot(bf)*(1.0+coreLsquared/RaSquared)+ Burgers.dot(DR)*bf.dot(DR)*1.0/RaSquared )*ruf.dot(rugauss.col(k))
            //					 )/sqrt(RaSquared);
            
            return -0.5*Material<Isotropic>::C2*temp*source.quadWeight*field.quadWeight; // C2=1/(4*pi*(1-nu))
        }
        
        
	};
    
    
    
    /**************************************************************************/
    /**************************************************************************/
}	// close namespace
#endif
