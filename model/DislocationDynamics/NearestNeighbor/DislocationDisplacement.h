/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _model_DislocationDisplacement_h_
#define _model_DislocationDisplacement_h_

#include <model/ParticleInteraction/FieldBase.h>


namespace model
{
	
	/**************************************************************************/
	/**************************************************************************/
	template<short unsigned int _dim>
	class DislocationDisplacement
    /* inheritance */ : public FieldBase<double,_dim,1>
    {
        
//        typedef typename DislocationParticleType::VectorDimD VectorDimD;
//        typedef typename DislocationParticleType::MatrixDim MatrixDim;
        
//        enum{dim=DislocationParticleType::dim};
//        enum{dataSize=dim*dim};
        
        
    public:
        

        
        typedef DislocationDisplacement<_dim> DislocationDisplacementType;
        typedef FieldBase<double,_dim,1> FieldBaseType;
        typedef typename FieldBaseType::MatrixType MatrixType;
        
        static  double a2;
		static const Eigen::Matrix<double,_dim,_dim> I;

        
        template <typename DislocationParticleType>
        static MatrixType compute(const DislocationParticleType& source,const DislocationParticleType& field)
        {
            
            Eigen::Matrix<double,_dim,1> R(field.P-source.P);
			double RaSquared (R.squaredNorm() + a2);
			return   (Material<Isotropic>::C1*(1.0+1.5*a2/RaSquared)*source.T*(source.B.cross(R)).transpose()
                               + 	R*(source.T.cross(source.B)).transpose()
                               +   0.5* R.cross(source.B).dot(source.T) * (I*(1.0+3.0*a2/RaSquared) + 3.0/RaSquared*R*R.transpose())
                               )/std::pow(sqrt(RaSquared),3)*source.quadWeight;
			
            //			return Material<Isotropic>::C2*(temp+temp.transpose());
//			return temp;
//
//            
//            return temp;
        }
        
        
//#ifdef _MODEL_MPI_
//        /**********************************************************************/
//        DislocationDisplacement(DislocationParticleType& dp1,const DislocationParticleType& dp2)
//        {/*! @param[in] dp1 A const reference to ChargedParticle 1
//          *  @param[in] dp2 A const reference to ChargedParticle 2
//          *
//          *  This constructor works with ParticleSystem<true,ChargedParticle>,
//          *  that is the paralle version of ParticleSystem. The result of the
//          *  interaction between dp1 and dp2 is stored in this->resultVector,
//          *  therefore dp1 and dp2 can be const.
//          */
//            
////            dp1._stress += dp2.stress_at(dp1.P);
//            MatrixDim temp(dp2.stress_at(dp1.P));
//            
//        }
//        
////        /**********************************************************************/
////        static ResultType get(const DislocationParticleType& dp1)
////        {
////            return ResultType((ResultType()<<InteractionBase<double,3>::resultVector[dp1.mpiID*3+0],
////                               /*         */ InteractionBase<double,3>::resultVector[dp1.mpiID*3+1],
////                               /*         */ InteractionBase<double,3>::resultVector[dp1.mpiID*3+2]).finished());
////        }
//        
//        /**********************************************************************/
//        static ResultType get(const DislocationParticleType& dp1)
//        {
//            return Eigen::Map<ResultType>(resultVector[dp1.mpiID*dataSize],dataSize);            
////            ((ResultType()<<InteractionBase<double,3>::resultVector[dp1.mpiID*3+0],
////                               /*         */ InteractionBase<double,3>::resultVector[dp1.mpiID*3+1],
////                               /*         */ InteractionBase<double,3>::resultVector[dp1.mpiID*3+2]).finished());
//        }
//        
//#else
//        /**********************************************************************/
//        DislocationDisplacement(DislocationParticleType& dp1, const DislocationParticleType& dp2)
//        {/*! Returns the asymmetric (and dimensionless) part of the stress integrand generated by the current quadrature point.
//          * @param[in] dp1	DislocationParticle 1 (field point)
//          * @param[in] dp2	DislocationParticle 1 (source point)
//          * \returns The stress field generated by dp2 on dp1
//          * 
//          * The return value is calculated according to:
//          * Cai, W., Arsenlis, A., Weinberger, C., & Bulatov, V. (2006). A non-singular continuum theory of dislocations. Journal Of The Mechanics And Physics Of Solids, 54(3), 561–587.
//          *	\f[
//          *		d\mathbf{s} = (1-\nu) \left(1+\frac{3}{2}\frac{a^2}{R_a^2}\right)\frac{\partial \mathbf{r}}{\partial u}\otimes \left(\mathbf{b}\times \mathbf{R}\right)+
//          *		\mathbf{R}\otimes\left(\frac{\partial \mathbf{r}}{\partial u}\times\mathbf{b}\right)+
//          *		\frac{1}{2} \left[ \left(\mathbf{R}\times\mathbf{b}\right)\cdot \frac{\partial \mathbf{r}}{\partial u} \right]\left[\mathbf{I}\left(1+\frac{3a^2}{R_a^2}\right)+\frac{3}{R_a^2}\mathbf{R}\otimes\mathbf{R}\right]
//          *	\f]
//          *  where \f$R_a^2=|\mathbf{R}|^2+a^2\f$ is the modified squared norm of \f$\mathbf{R}\f$.
//          *
//          *	The return value is asymmetric and dimensionless in the sense that the actual stress integrand is:
//          *	\f[
//          *		d\mathbf{\sigma}=\frac{\mu}{4\pi (1-\nu)}\left(d\mathbf{s}+d\mathbf{s}^T\right)
//          *	\f]
//          */
////			const VectorDimD R(dp1.P-dp2.P);
////			const double RaSquared (R.squaredNorm() + a2);
//            //            const VectorDimD BcrossR(dp2.B.cross(R));
////            dp1._stress += (Material<Isotropic>::C1*(1.0+1.5*a2/RaSquared)*dp2.T*(dp2.B.cross(R)).transpose()
////                            + 	R*(dp2.T.cross(dp2.B)).transpose()
////                            +   0.5* R.cross(dp2.B).dot(dp2.T) * (I*(1.0+3.0*a2/RaSquared) + 3.0/RaSquared*R*R.transpose())
////                            )/std::pow(sqrt(RaSquared),3)*dp2.quadWeight;
//
//                        dp1._stress += dp2.stress_at(dp1.P);
//
//            
////            MatrixDim temp((Material<Isotropic>::C1*(1.0+1.5*a2/RaSquared)*dp2.T*(dp2.B.cross(R)).transpose()
////                            + 	R*(dp2.T.cross(dp2.B)).transpose()
////                            +   0.5* R.cross(dp2.B).dot(dp2.T) * (I*(1.0+3.0*a2/RaSquared) + 3.0/RaSquared*R*R.transpose())
////                            )/std::pow(sqrt(RaSquared),3)*dp2.quadWeight);
////            
////            dp1._stress += temp;
//            
//        }
//        
////        /**********************************************************************/
////        static void compute(DislocationParticleType& dp1, const DislocationParticleType& dp2)
////        {/*! Returns the asymmetric (and dimensionless) part of the stress integrand generated by the current quadrature point.
////          * @param[in] dp1	DislocationParticle 1 (field point)
////          * @param[in] dp2	DislocationParticle 1 (source point)
////          * \returns The stress field generated by dp2 on dp1
////          *
////          * The return value is calculated according to:
////          * Cai, W., Arsenlis, A., Weinberger, C., & Bulatov, V. (2006). A non-singular continuum theory of dislocations. Journal Of The Mechanics And Physics Of Solids, 54(3), 561–587.
////          *	\f[
////          *		d\mathbf{s} = (1-\nu) \left(1+\frac{3}{2}\frac{a^2}{R_a^2}\right)\frac{\partial \mathbf{r}}{\partial u}\otimes \left(\mathbf{b}\times \mathbf{R}\right)+
////          *		\mathbf{R}\otimes\left(\frac{\partial \mathbf{r}}{\partial u}\times\mathbf{b}\right)+
////          *		\frac{1}{2} \left[ \left(\mathbf{R}\times\mathbf{b}\right)\cdot \frac{\partial \mathbf{r}}{\partial u} \right]\left[\mathbf{I}\left(1+\frac{3a^2}{R_a^2}\right)+\frac{3}{R_a^2}\mathbf{R}\otimes\mathbf{R}\right]
////          *	\f]
////          *  where \f$R_a^2=|\mathbf{R}|^2+a^2\f$ is the modified squared norm of \f$\mathbf{R}\f$.
////          *
////          *	The return value is asymmetric and dimensionless in the sense that the actual stress integrand is:
////          *	\f[
////          *		d\mathbf{\sigma}=\frac{\mu}{4\pi (1-\nu)}\left(d\mathbf{s}+d\mathbf{s}^T\right)
////          *	\f]
////          */
////			const VectorDimD R(dp1.P-dp2.P);
////			double RaSquared (R.squaredNorm() + a2);
////            
////            dp1._stress += (Material<Isotropic>::C1*(1.0+1.5*a2/RaSquared)*dp2.T*(dp2.B.cross(R)).transpose()
////                            + 	R*(dp2.T.cross(dp2.B)).transpose()
////                            +   0.5* R.cross(dp2.B).dot(dp2.T) * (I*(1.0+3.0*a2/RaSquared) + 3.0/RaSquared*R*R.transpose())
////                            )/std::pow(RaSquared,1.5)*dp2.quadWeight;
////            
////            
////            //			MatrixDim stress(   (Material<Isotropic>::C1*(1.0+1.5*a2/RaSquared)*T*(B.cross(R)).transpose()
////            //                               + 	R*(T.cross(B)).transpose()
////            //                               +   0.5* R.cross(B).dot(T) * (I*(1.0+3.0*a2/RaSquared) + 3.0/RaSquared*R*R.transpose())
////            //                               )/std::pow(RaSquared,1.5)*quadWeight);
////            
////            //            const VectorDimD BcrossR(dp2.B.cross(R));
////            //			MatrixDim stress(   (Material<Isotropic>::C1*(1.0+1.5*a2/RaSquared)*dp2.T*(BcrossR).transpose()
////            //                               + 	R*(dp2.T.cross(dp2.B)).transpose()
////            //                               -   0.5* BcrossR.dot(dp2.T) * (I*(1.0+3.0*a2/RaSquared) + 3.0/RaSquared*R*R.transpose())
////            //                               )/std::pow(RaSquared,1.5)*dp2.quadWeight);
////			
////            //            dp1._stress += stress;
////        }
//        
//        /**********************************************************************/
//        static const MatrixDim& get(const DislocationParticleType& dp1)
//        {
////            return dp1.force;
//        }
//#endif
				
	};
    
    // Static data members
	template<short unsigned int _dim>
    double DislocationDisplacement<_dim>::a2=1.0;  // square of core size a
    
	template<short unsigned int _dim>
	const Eigen::Matrix<double,_dim,_dim> DislocationDisplacement<_dim>::I=Eigen::Matrix<double,_dim,_dim>::Identity();  // square of core size a

    
    /**************************************************************************/
    /**************************************************************************/
}	// close namespace
#endif
