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
	
	/*!\brief Class template that implements the calculation of the stress field
     * generated by a dislocation. Depending on the value of the macro
     * _MODEL_NON_SINGULAR_DD_, different dislocation theories are used for the 
     * calculation.
     */
	template<short unsigned int _dim>
	struct DislocationStress
    /* inheritance */ : public FieldBase<double,_dim,_dim>
    {
        
        constexpr static int dim=_dim;
        typedef DislocationStress<_dim> DislocationStressType;
        typedef FieldBase<double,_dim,_dim> FieldBaseType;
        typedef typename FieldBaseType::MatrixType MatrixType;
        
        //! Dislocation core size
        static  double a;
        
        //! Dislocation core size squared
        static  double a2;
        
        //! dim x dim identity matrix
		static const Eigen::Matrix<double,_dim,_dim> I;
        
        
#if _MODEL_NON_SINGULAR_DD_ == 0 // Note that if _MODEL_NON_SINGULAR_DD_ is not #defined, the preprocessor treats it as having the value 0.
        
        template <typename DislocationParticleType, typename OtherParticleType>
        static MatrixType compute(const DislocationParticleType& source,const OtherParticleType& field)
        {/*!@param[in] source the DislocationParticle that is the source of stress
          * @param[in] field  the FieldPoint at which stress is computed
          *\returns the Cauchy stress contribution of source on field. This
          * function is used when the macro _MODEL_NON_SINGULAR_DD_ is set to 0
          * (classical elasticity is used).
          *
          * Note that the return value is the asymmetric and dimensionless stress
          * quantity:
          *	\f[
          * \mathbf s=\frac{1}{R_a^2}\left\{
          * -C_1\hat{\mathbf\xi}'\otimes(\hat{\mathbf R}\times\mathbf b)
          * -\hat{\mathbf R}\otimes(\mathbf b\times\hat{\mathbf \xi}')
          * +\frac{1}{2}\hat{\mathbf R}\cdot(\mathbf b\times\hat{\mathbf\xi}')\left[\hat{3\mathbf R}\otimes\hat{\mathbf R}+\mathbf I\right]\right\}\ dL'
          *	\f]
          * from which the stress is computed as
          *	\f[
          *		\mathbf{\sigma}(\mathbf x)=C_2\oint_{\mathcal{L}}\left(\mathbf s+\mathbf s^T \right)
          *	\f]
          *  where \f$\mathbf R=\mathbf x'-\mathbf{x}\f$, \f$R_a=\sqrt{R^2+a^2}\f$, 
          * and \f$\hat{\mathbf R}=\mathbf R/R_a\f$. The parameter \f$a\f$
          * is used to regularize the stress field. The exact (singular) classical theory is
          * obtained for \f$a=0\f$.
          *
          * The derivation of this quantity follows from [1]:
          *	\f[
          * \begin{align}
          * \mathbf{\sigma}(\mathbf x)
          * &=\frac{\mu}{2\pi(1-\nu)}\oint_{\mathcal{L}}\frac{1}{R_a^2}\left\{
          * -\frac{1-\nu}{2}\left[\hat{\mathbf\xi}'\otimes(\hat{\mathbf R}\times\mathbf b)+(\hat{\mathbf R}\times\mathbf b)\otimes\hat{\mathbf\xi}'\right]
          * -\frac{1}{2}\left[(\mathbf b\times\hat{\mathbf \xi}')\otimes\hat{\mathbf R}+\hat{\mathbf R}\otimes(\mathbf b\times\hat{\mathbf \xi}')\right]
          * +\frac{1}{2}\hat{\mathbf R}\cdot(\mathbf b\times\hat{\mathbf\xi}')\left[3\hat{\mathbf R}\otimes\hat{\mathbf R}+\mathbf I\right]
          * \right\}\ dL'\\
          * &=\underbrace{\frac{\mu}{4\pi(1-\nu)}}_{C_2}\oint_{\mathcal{L}}\frac{1}{R_a^2}\left\{
          * -\underbrace{(1-\nu)}_{C_1}\left[\hat{\mathbf\xi}'\otimes(\hat{\mathbf R}\times\mathbf b)+(\hat{\mathbf R}\times\mathbf b)\otimes\hat{\mathbf\xi}'\right]
          * -\left[(\mathbf b\times\hat{\mathbf \xi}')\otimes\hat{\mathbf R}+\hat{\mathbf R}\otimes(\mathbf b\times\hat{\mathbf \xi}')\right]
          * +\hat{\mathbf R}\cdot(\mathbf b\times\hat{\mathbf\xi}')\left[3\hat{\mathbf R}\otimes\hat{\mathbf R}+\mathbf I\right]
          * \right\}\ dL'\\
          * &=C_2\oint_{\mathcal{L}}\left(\mathbf s+\mathbf s^T \right)
          * \end{align}
          *	\f]
          *
          * References:
          *
          * [1] de Wit, R. (1960). The Continuum Theory of Stationary 
          * Dislocations. Solid State Physics, 10, 249–292.
          */
            
            MatrixType temp(MatrixType::Zero());
            Eigen::Matrix<double,_dim,1> r(field.P-source.P);
            const double R2(r.squaredNorm()+a2);
                r/=sqrt(R2); // normalize r
                temp = -Material<Isotropic>::C1*source.T*(r.cross(source.B)).transpose()
                /*  */ -r*(source.B.cross(source.T)).transpose() // the factor 2.0 takes care of the C2 term in DislocationParticle::stress
                /*  */ +0.5*r.dot(source.B.cross(source.T)) * (3.0*r*r.transpose() + I);
            return temp/R2*source.quadWeight;
        }
        
#elif _MODEL_NON_SINGULAR_DD_ == 1 /* Cai's non-singular theory */

        template <typename DislocationParticleType, typename OtherParticleType>
        static MatrixType compute(const DislocationParticleType& source, const OtherParticleType& field)
        {/*!@param[in] source the DislocationParticle that is the source of stress
          * @param[in] field  the FieldPoint at which stress is computed
          *\returns the effective stress contribution of source on field. This
          * function is used when the macro _MODEL_NON_SINGULAR_DD_ is set to 1
          * (Cai's non-singular theory is used).
          *
          * Note that the return value is the asymmetric and dimensionless stress
          * quantity:
          *	\f[
          * ????
          *	\f]
          * from which the stress is computed as
          *	\f[
          *		\mathbf{\sigma}(\mathbf x)=C_2\oint_{\mathcal{L}}\left(\mathbf s+\mathbf s^T \right)
          *	\f]
          *  where \f$\mathbf R=\mathbf x'-\mathbf{x}\f$, \f$R_a=\sqrt{R^2+a^2}\f$,
          * and \f$\hat{\mathbf R}=\mathbf R/R_a\f$.
          *
          * The derivation of this quantity follows from [2]:
          *	\f[
          *
          *	\f]
          *
          * References:
          *
          * [2] Cai, W.et al (2006). A non-singular continuum theory of dislocations.
          * Journal of the Mechanics and Physics of Solids, 54(3), 561–587.
          */
            
            Eigen::Matrix<double,_dim,1> R(field.P-source.P);
            double RaSquared (R.squaredNorm() + a2);
            return   (Material<Isotropic>::C1*(1.0+1.5*a2/RaSquared)*source.T*(source.B.cross(R)).transpose()
                      + 	R*(source.T.cross(source.B)).transpose()
                      +   0.5* R.cross(source.B).dot(source.T) * (I*(1.0+3.0*a2/RaSquared) + 3.0/RaSquared*R*R.transpose())
                      )/std::pow(sqrt(RaSquared),3)*source.quadWeight;
        }
        
#elif _MODEL_NON_SINGULAR_DD_ == 2 /* Lazar's non-singular theory */
        
        template <typename DislocationParticleType, typename OtherParticleType>
        static MatrixType compute(const DislocationParticleType& source,const OtherParticleType& field)
        {/*!@param[in] source the DislocationParticle that is the source of stress
          * @param[in] field  the FieldPoint at which stress is computed
          *\returns the Cauchy stress contribution of source on field. This 
          * function is used when the macro _MODEL_NON_SINGULAR_DD_ is set to 2
          * (gradient elasticity of Helmholtz type is used).
          *
          * Note that the return value is the asymmetric and dimensionless stress 
          * quantity:
          *	\f[
          * \mathbf s=\frac{1}{\ell^2}\left\{
          * C_1\hat{\mathbf\xi}'\otimes(\hat{\mathbf R}\times\mathbf b)\ f^*_4(R/\ell) +
          * \hat{\mathbf R}\otimes(\mathbf b\times\hat{\mathbf \xi}')\ 2f^*_5(R/\ell) +
          * \hat{\mathbf R}\cdot(\mathbf b\times\hat{\mathbf\xi}')\left[\hat{\mathbf R}\otimes\hat{\mathbf R}f_6(R/\ell)+\mathbf If_7(R/\ell)\right]\right\}\ dL'
          *	\f]
          * from which the stress is computed as
          *	\f[
          *		\mathbf{\sigma}(\mathbf x)=C_2\oint_{\mathcal{L}}\left(\mathbf s+\mathbf s^T \right)
          *	\f]
          *  where \f$\mathbf R=\mathbf x'-\mathbf{x}\f$ and \f$\hat{\mathbf R}=\mathbf R/R\f$.
          *
          * The derivation of this quantity follows from [3]:
          *	\f[
          * \begin{align}
          * \mathbf{\sigma}(\mathbf x)
          * &=\frac{\mu}{2\pi(1-\nu)\ell^2}\oint_{\mathcal{L}}\left\{
          * \frac{1-\nu}{2}\left[\hat{\mathbf\xi}'\otimes(\hat{\mathbf R}\times\mathbf b)+(\hat{\mathbf R}\times\mathbf b)\otimes\hat{\mathbf\xi}'\right]f^*_4(R/\ell)
          * +\left[(\mathbf b\times\hat{\mathbf \xi}')\otimes\hat{\mathbf R}+\hat{\mathbf R}\otimes(\mathbf b\times\hat{\mathbf \xi}')\right]f_5(R/\ell)
          * +\hat{\mathbf R}\cdot(\mathbf b\times\hat{\mathbf\xi}')\left[\hat{\mathbf R}\otimes\hat{\mathbf R}f_6(R/\ell)+\mathbf If_7(R/\ell)\right]
          * \right\}\ dL'\\
          * &=\underbrace{\frac{\mu}{4\pi(1-\nu)}}_{C_2}\frac{1}{\ell^2}\oint_{\mathcal{L}}\left\{
          * \underbrace{(1-\nu)}_{C_1}\left[\hat{\mathbf\xi}'\otimes(\hat{\mathbf R}\times\mathbf b)+(\hat{\mathbf R}\times\mathbf b)\otimes\hat{\mathbf\xi}'\right]f^*_4(R/\ell)
          * +\left[(\mathbf b\times\hat{\mathbf \xi}')\otimes\hat{\mathbf R}+\hat{\mathbf R}\otimes(\mathbf b\times\hat{\mathbf \xi}')\right]2f_5(R/\ell)
          * +2\hat{\mathbf R}\cdot(\mathbf b\times\hat{\mathbf\xi}')\left[\hat{\mathbf R}\otimes\hat{\mathbf R}f_6(R/\ell)+\mathbf If_7(R/\ell)\right]
          * \right\}\ dL'\\
          * &=C_2\oint_{\mathcal{L}}\left(\mathbf s+\mathbf s^T \right)
          * \end{align}
          *	\f]
          * In the previous equations, the scalar functions \f$f^*_4\ldots f^*_7\f$ are:
          *	\f[
          * \begin{align}
          * f^*_4(x)&=\frac{df^*_1}{dx}=\frac{(1+x)\text{e}^{-x} - 1 }{x^2}\\
          * f^*_5(x)&=\frac{df^*_2}{dx}=\frac{f^*_3}{x}=\frac{3[1-(1+x)\text{e}^{-x}] - x^2\, \left( \frac{1}{2}+\text{e}^{- x}\right) }{x^4}\\
          * f^*_6(x)&=f^*_4(x)-5f^*_5(x)\\
          * f^*_7(x)&=f^*_5(x)-f^*_4(x)
          * \end{align}
          *	\f]
          *
          * References:
          *
          * [3] Po, G et al (2014). Singularity-free
          * Dislocation Dynamics with strain gradient elasticity.
          * Journal Of The Mechanics And Physics Of Solids, 68  161–178.
          */
            
            MatrixType temp(MatrixType::Zero());
            Eigen::Matrix<double,_dim,1> r(field.P-source.P);
            const double R(r.norm());
            if (R>0.0)
            {
                const double x = R/a;       // this is R/L
                const double x2 = std::pow(x ,2);
                const double ex = exp(-x); // this is exp(-R/L)
                const double f4 = ((1.0+x)*ex-1.0)/x2; // as in Po et at 2014
                const double f5 = (3.0*(1.0-(1.0+x)*ex)-x2*(0.5+ex))/std::pow(x2,2);
                r/=R; // normalize r
                temp = Material<Isotropic>::C1*source.T*(r.cross(source.B)).transpose()*f4
                /*  */ +r*(source.B.cross(source.T)).transpose()*2.0*f5
                /*  */ +r.dot(source.B.cross(source.T)) * (r*r.transpose()*(f4-5.0*f5) + I*(f5-f4));
            }
            return temp/a2*source.quadWeight;
        }
        
#else
#error Unsupported choice of field regularization
#endif
        
	};
    
    /**************************************************************************/
    // Static data members
    
    //! Dislocation core size
    template<short unsigned int _dim>
    double DislocationStress<_dim>::a=1.0;
    
    //! Dislocation core size squared
	template<short unsigned int _dim>
    double DislocationStress<_dim>::a2=1.0;
    
    //! Identity matrix
	template<short unsigned int _dim>
	const Eigen::Matrix<double,_dim,_dim> DislocationStress<_dim>::I=Eigen::Matrix<double,_dim,_dim>::Identity();
    
}	// close namespace
#endif
