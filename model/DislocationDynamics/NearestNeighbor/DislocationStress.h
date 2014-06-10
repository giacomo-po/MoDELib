/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _model_DislocationStress_h_
#define _model_DislocationStress_h_

#include <model/ParticleInteraction/FieldBase.h>
#include <model/ParticleInteraction/FieldPoint.h>


namespace model
{
	
	/**************************************************************************/
	/**************************************************************************/
	template<short unsigned int _dim>
	class DislocationStress
    /* inheritance */ : public FieldBase<double,_dim,_dim>
    {        
        
    public:
        
        
        
        typedef DislocationStress<_dim> DislocationStressType;
        typedef FieldBase<double,_dim,_dim> FieldBaseType;
        typedef typename FieldBaseType::MatrixType MatrixType;
        
        static  double a;
        static  double a2;
		static const Eigen::Matrix<double,_dim,_dim> I;
        
        
        
#if _MODEL_NON_SINGULAR_DD_ == 1 /* Cai's non-singular theory */
        template <typename DislocationParticleType>
        static MatrixType compute(const DislocationParticleType& source,const DislocationParticleType& field)
        {/*!@param[in] source the DislocationParticle that is source of stress
          * @param[in] field  the DislocationParticle on which stress is computed
          *\returns the stress field produced by source on field
          */
            
            Eigen::Matrix<double,_dim,1> R(field.P-source.P);
			double RaSquared (R.squaredNorm() + a2);
			return   (Material<Isotropic>::C1*(1.0+1.5*a2/RaSquared)*source.T*(source.B.cross(R)).transpose()
                      + 	R*(source.T.cross(source.B)).transpose()
                      +   0.5* R.cross(source.B).dot(source.T) * (I*(1.0+3.0*a2/RaSquared) + 3.0/RaSquared*R*R.transpose())
                      )/std::pow(sqrt(RaSquared),3)*source.quadWeight;
        }
        
        template <typename DislocationParticleType, typename OtherParticleType>
        static MatrixType compute(const DislocationParticleType& source, const OtherParticleType& field)
        {/*!@param[in] source the DislocationParticle that is source of stress
          * @param[in] field  the DislocationParticle on which stress is computed
          *\returns the stress field produced by source on field
          */
            
            Eigen::Matrix<double,_dim,1> R(field.P-source.P);
            double RaSquared (R.squaredNorm() + a2);
            return   (Material<Isotropic>::C1*(1.0+1.5*a2/RaSquared)*source.T*(source.B.cross(R)).transpose()
                      + 	R*(source.T.cross(source.B)).transpose()
                      +   0.5* R.cross(source.B).dot(source.T) * (I*(1.0+3.0*a2/RaSquared) + 3.0/RaSquared*R*R.transpose())
                      )/std::pow(sqrt(RaSquared),3)*source.quadWeight;
            
        }
        
#elif _MODEL_NON_SINGULAR_DD_ == 2 /* Lazar's non-singular theory */
        
//        static double f4(const double& x)
//        {
//            return (RL>0.0)? ((1.0+x)*exp(-x)-1.0)/std::pow(x,2) : -0.5;
//        }
        
        template <typename DislocationParticleType>
        static MatrixType compute(const DislocationParticleType& source,const DislocationParticleType& field)
        {/*!@param[in] source the DislocationParticle that is source of stress
          * @param[in] field  the DislocationParticle on which stress is computed
          *\returns the stress field produced by source on field
          */
            
            MatrixType temp(MatrixType::Zero());
            
            Eigen::Matrix<double,_dim,1> r(field.P-source.P);
            const double R2(r.squaredNorm());
            if (R2>0.0)
            {
                const double R(sqrt(R2));
                
                double f4(1.0/R2);
                double F2(f4);
                double F3(f4);
                
                
                const double RL(R/a);       // this is R/L
                const double LR(a/R);       // this is L/R
                const double eRL(exp(-RL)); // this is exp(-R/L)
                f4*=(1.0-(1.0+RL)*eRL);
                F2*=(1.0-6.0*a2/R2*(1.0-eRL)+(2.0+6.0*LR)*eRL);
                F3*=(1.0-10.0*a2/R2*(1.0-eRL)+(4.0+10.0*LR+2.0/3.0*RL)*eRL);
                
                r/=R; // normalize r
                temp = Material<Isotropic>::C1*source.T*(source.B.cross(r)).transpose()*f4
                /*  */ +r*(source.T.cross(source.B)).transpose()*F2
                /*  */ +0.5* r.cross(source.B).dot(source.T) * (I*(2.0*f4-F2) + 3.0*F3*r*r.transpose());
            }
            
            return temp*source.quadWeight;
        }
        
#else // Note that if _MODEL_NON_SINGULAR_DD_ is not #defined, the preprocessor treats it as having the value 0.
#error Unsupported choice of field regularization
#endif
        
	};
    
    template<short unsigned int _dim>
    double DislocationStress<_dim>::a=1.0;  // square of core size a
    
    // Static data members
	template<short unsigned int _dim>
    double DislocationStress<_dim>::a2=1.0;  // square of core size a
    
	template<short unsigned int _dim>
	const Eigen::Matrix<double,_dim,_dim> DislocationStress<_dim>::I=Eigen::Matrix<double,_dim,_dim>::Identity();  // square of core size a
    
    
    /**************************************************************************/
    /**************************************************************************/
}	// close namespace
#endif
