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
#include <model/DislocationDynamics/ElasticFields/DislocationStress.h>

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
        typedef Eigen::Matrix<double,_dim,_dim> MatrixDim;
        typedef Eigen::Matrix<double,_dim,1>   VectorDim;
        
        static bool use_multipole;
        
        /**********************************************************************/
        static const MatrixType& get(const MatrixType& temp)
        {
            return temp;
        }
        
        /**********************************************************************/
        template <typename ParticleType, typename CellContainerType>
        static MatrixType multipole(const ParticleType& ,const CellContainerType&)
        {/*!@param[in] field  the FieldPoint at which stress is computed
          * @param[in] farCells container of SpatialCell(s) that are not neighbors of field
          *\returns the displacement contribution of the farCells on field.
          *
          */
            assert(0 && "Multiple expansion of energy not implemented yet"); // \todo Finish implementation here

            MatrixType temp(MatrixType::Zero());

            return temp;
        }
        
        
#if _MODEL_NON_SINGULAR_DD_ == 0 // Note that if _MODEL_NON_SINGULAR_DD_ is not #defined, the preprocessor treats it as having the value 0.
        
        template <typename DislocationParticleType>
        static MatrixType compute(const DislocationParticleType& source,const DislocationParticleType& field)
        {/*!@param[in] source the DislocationParticle that is source of stress
          * @param[in] field  the DislocationParticle on which stress is computed
          *\returns the stress field produced by source on field
          */
            VectorDim R(field.P-source.P);
			const double Ra=sqrt(R.squaredNorm()+DislocationStress<_dim>::a2);
            R/=Ra; // normalize R
            return (MatrixType()<<-Material<Isotropic>::C2*source.quadWeight*field.quadWeight/Ra*
            (Material<Isotropic>::C1*source.B.dot(source.T)*field.B.dot(field.T)
             +2.0*Material<Isotropic>::nu*(field.B.dot(source.T)*source.B.dot(field.T))
             -(source.B.dot(field.B)+ source.B.dot(R)*field.B.dot(R))*field.T.dot(source.T)
             )).finished();
        }
        
#elif _MODEL_NON_SINGULAR_DD_ == 1 /* Cai's non-singular theory */
        
        template <typename DislocationParticleType>
//        static MatrixType compute(const DislocationParticleType& source,const DislocationParticleType& field)
        static MatrixType compute(const DislocationParticleType& ,const DislocationParticleType& )
        {/*!@param[in] source the DislocationParticle that is source of stress
          * @param[in] field  the DislocationParticle on which stress is computed
          *\returns the stress field produced by source on field
          */
            
            assert(0 && "FINISH IMPLEMENTATION HERE"); // \todo Finish implementation here
            
            //            MatrixType temp(MatrixType::Zero());
            //
            //            const Eigen::Matrix<double,_dim,1> r(field.P-source.P);
            //
            //
            //
            
            //            return -0.5*Material<Isotropic>::C2*source.quadWeight*field.quadWeight*
            //            (Material<Isotropic>::C1*(1.0+0.5*coreLsquared/RaSquared)*source.B.dot(source.T)*field.B.dot(field.T)
            //                     +2.0*Material<Isotropic>::nu*(1.0+0.5*coreLsquared/RaSquared)*(field.B.dot(source.T)*source.B.dot(field.T))
            //                     -(source.B.dot(field.B)*(1.0+coreLsquared/RaSquared)+ source.B.dot(DR)*field.B.dot(DR)*1.0/RaSquared )*field.T.dot(source.T)
            //                     )/sqrt(RaSquared);
            
        }
        
#elif _MODEL_NON_SINGULAR_DD_ == 2 /* Lazar's non-singular theory */
        
        template <typename DislocationParticleType>
        static MatrixType compute(const DislocationParticleType& source,const DislocationParticleType& field)
        {/*!@param[in] source the DislocationParticle that is source of stress
          * @param[in] field  the DislocationParticle on which stress is computed
          *\returns the stress field produced by source on field
          */
            
            MatrixType temp(MatrixType::Zero());
            
            const VectorDim r(field.P-source.P);
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
            
            
            return -0.5*Material<Isotropic>::C2*temp*source.quadWeight*field.quadWeight; // C2=1/(4*pi*(1-nu))
        }
        
#else
#error Unsupported choice of field regularization
#endif
        
        
	};
    
    template<short unsigned int _dim>
	bool DislocationEnergy<_dim>::use_multipole=true;

    
    /**************************************************************************/
    /**************************************************************************/
}	// close namespace
#endif
