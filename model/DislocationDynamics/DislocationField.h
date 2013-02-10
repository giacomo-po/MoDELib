/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DISLOCATIONFIELD_H_
#define model_DISLOCATIONFIELD_H_


//#include <vector>
#include <Eigen/Dense>


namespace model {
	
	struct Isotropic{};
	
	/****************************************************/
	/* LinearElasticGreensFunction: general case ********/
	/****************************************************/
	template <short unsigned int dim, typename Type>
	struct DislocationField {
		
		
	};
	
	
	/****************************************************/
	/* LinearElasticGreensFunction: isotropic case ******/
	/****************************************************/
	template <short unsigned int dim>
	struct DislocationField<dim,Isotropic> {
		
		typedef Eigen::Matrix<double,dim,dim> MatrixDim;
		typedef Eigen::Matrix<double,dim,1>   VectorDim;
		
		static const Eigen::Matrix<double,dim,dim> I;
		static const double C1;
		static const double C2;
		static const double C3;
		static const double C4;
		
		
		/* displacementKernel ****************************/
		VectorDim displacementKernel(const int & k, const VectorDim & Rfield, const VectorDim& S) const{	
			/*! The infinitesimal dispacement field generated by this segment at a field point
			 * @param[in] k			the current quadrature point
			 * @param[in] Rfield	the field point
			 *
			 * The return value is calculated according to:
			 *	\f[
			 *		\mathbf{u}^*(\alpha_k,\mathbf{r}_f) = \frac{1}{R}\left\{ \frac{2(1-\nu)}{R+\mathbf{R}\cdot\mathbf{S}}\mathbf{b} \left[\left(\mathbf{S}\times\mathbf{R}\right)\cdot\frac{d\mathbf{r}_s}{d\alpha}\right] 
			 *                                          + (1-2\nu)\left( \mathbf{b}\times \frac{d\mathbf{r}_s}{d\alpha}\right)
			 *                                          +\frac{1}{R^2}\left[\left(\frac{d\mathbf{r}_s}{d\alpha}\times\mathbf{b}\right)\cdot\mathbf{R}\right]\mathbf{R} \right\}
			 *	\f]
			 *	where:
			 *  - \f$\mathbf{r}(\alpha)\f$ is the parametrized source line segment;
			 *  - \f$\mathbf{R}(\mathbf{r}_f,\alpha) = \mathbf{r}_f-\mathbf{r}_s(\alpha)\f$ is vector connecting the source point to the field point.
			 *
			 *  The calculation of the solid angle is performed transforming the surface integral to a line integral. From Stokes theorem:
			 *	\f[
			 *  \oint\frac{\hat{\mathbf{S}}\times\mathbf{R}}{R(R-\mathbf{S}\cdot\mathbf{R})}\cdot d\mathbf{l}
			 *  \underbrace{=}_{\mbox{Stokes Th.}}\int\left(\nabla\times\frac{\hat{\mathbf{S}}\times\mathbf{R}}{R(R-\mathbf{S}\cdot\mathbf{R})}\right)\cdot\hat{\mathbf{n}}dA
			 *  \underbrace{=}_{\mbox{identity}}
			 * -\int\frac{\mathbf{R}}{R^3}\cdot\mathbf{n}dA 
			 *	\f]
			 * References:
			 * [1] Asvestas, J. Line integrals and physical optics. Part I. The transformation of the solid-angle surface integral to a line integral. J. Opt. Soc. Am. A, 2(6), 891–895.
			 */
			
			VectorDim R=Rfield-rgauss.col(k);	
			double Ra=std::pow(R.squaredNorm()+std::pow(coreL,2.0),0.5);			
			return 1.0/Ra * (+ 2.0*C1/(Ra+R.dot(S))*Burgers*(S.cross(R)).dot(rugauss.col(k))
			/*            */ + C3*rugauss.col(k).cross(Burgers) 
			/*            */ + 1.0/std::pow(Ra,2)*(rugauss.col(k).cross(Burgers)).dot(R)*R
			/*            */ );
		}
		
		
		/* stressKernel ********************************/
		static MatrixDim stressKernel(const VectorDim& Rfield, const VectorDim& Rsource, const VectorDim& Tsource, const VectorDim& Bsource, const double& aSquared){
			const VectorDim R(Rfield-Rsource);
			const double RaSquared(R.squaredNorm() + aSquared);

			return   (C1*(1.0+1.5*aSquared/RaSquared)*Tsource*(Bsource.cross(R)).transpose()
					  + R*(Tsource.cross(Bsource)).transpose() 
					  + 0.5* R.cross(Bsource).dot(Tsource) * (I*(1.0+3.0*aSquared/RaSquared) + 3.0/RaSquared*R*R.transpose())
					  )/std::pow(RaSquared,1.5);		
		}
		
		/* energyKernel ********************************/
		static double energyKernel(const VectorDim& Rfield, const VectorDim& Rsource, const VectorDim& Tsource, const VectorDim& Bsource, const double& aSquared){
			const VectorDim R(Rfield-Rsource);
			const double RaSquared(R.squaredNorm() + aSquared);
			return  (C1*(1.0+0.5*aSquared/RaSquared)*Bsource.dot(Tsource)*bf.dot(ruf)
					 +2.0*shared.material.nu*(1.0+0.5*aSquared/RaSquared)*(bf.dot(Tsource)*Bsource.dot(ruf))
					 -(Bsource.dot(bf)*(1.0+aSquared/RaSquared)+ Bsource.dot(DR)*bf.dot(DR)*1.0/RaSquared )*ruf.dot(Tsource)
					 )/std::pow(RaSquared,0.5);
		}
		

		/********************************************************/
		MatrixDim stress_straight(const VectorDim & R, const VectorDim & t) const {
			
			const VectorDim T=t.normalized();
			const double RdotT(R.dot(T));
			const double RdotT2(std::pow(RdotT,2));
			const double RaSquared(R.squaredNorm()+coreLsquared);
			const double Ra = std::pow(RaSquared,0.5);
			const double RaCubed(std::pow(Ra,3));
			const double A1 = - RdotT*(3.0*RaSquared-RdotT2)/std::pow(RaSquared-RdotT2,2)/RaCubed;
			const double A2 = 1.0/RaCubed-RdotT*A1;
			const double A6 = - RdotT/(RaSquared-RdotT2)/Ra;
			const double A3 = - RdotT/Ra+A6+RdotT2*A1;
			const double A4 = A6 + coreLsquared*A1;
			const double A5 = -C1*A6-coreLsquared*C1*A1*0.5;
			const double A7 = shared.material.nu/Ra - RdotT*A6 - coreLsquared*C1*A2*0.5;
			
			return R.cross(Burgers).dot(t)*(0.5*A1*R*R.transpose() + A2*t*R.transpose() + 0.5*A3*t*t.transpose() + 0.5*A4*I)
			/*  */ + A5*R.cross(Burgers)*t.transpose() + A6*t.cross(Burgers)*R.transpose()
			/*  */ + A7*t.cross(Burgers)*t.transpose();		// extract SYMMETRIC part of stress
		}
		
		/********************************************************/
		MatrixDim stress_straight_inf(const VectorDim & t) const {
			const VectorDim T=t.normalized();
			
			return MatrixDim::Zero();		// extract SYMMETRIC part of stress
		}


		
		
	};
	
	
	// Declare static data members
	template <short unsigned int dim, typename Type>
	const Eigen::Matrix<double,dim,dim> DislocationField<dim,Type>::I=Eigen::Matrix<double,dim,dim>::Identity();

	template <short unsigned int dim, typename Type>
	const double DislocationField<dim,Type>::C1=1.0-shared.material.nu;
	
	template <short unsigned int dim, typename Type>
	const double DislocationField<dim,Type>::C2=shared.material.mu/(4.0*M_PI*C1);
	
	template <short unsigned int dim, typename Type>
	const double DislocationField<dim,Type>::C3=1.0-2.0*shared.material.nu;
	
	template <short unsigned int dim, typename Type>
	const double DislocationField<dim,Type>::C4=1.0/(8.0*M_PI*C1);
	
	/********************************************************************/
	/********************************************************************/
} // end namespace
#endif