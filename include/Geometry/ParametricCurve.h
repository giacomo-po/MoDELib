/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2011 by Benjamin Ramirez<ramirezbrf@gmail.com>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PARAMETRICCURVE_H_
#define model_PARAMETRICCURVE_H_

#include <Eigen/Dense>

#include <Quadrature.h>


namespace model {
	

	
	
	/*!\brief Class template that computes tangent, curvature and torsion of a parametric curve \f$\mathbf{r}(u)\f$ in arbitrary dimension.
	 * 
	 * Let \f$u\in[0,1]\f$ and \f$\mathbf{r}(u): [0,1]\rightarrow \mathbb{R}^d\f$ a parametric curve. 
	 * The first derivative \f$\mathbf{r}_{,u}\f$ with respect to the parameter \f$u\f$ represents the (non unitary) 
	 * parametric tangent to the curve. The arc length \f$l\f$ along the curve can be calculated as:
	 * of the parameter \f$u\f$ at each time \f$t\f$:
	 * \f[
	 * l(u)=\int_0^u\sqrt{\mathbf{r}_{,v}\cdot\mathbf{r}_{,v}}dv
	 * \f]
	 * The length of the curve is therefore:
	 * \f[
	 * L=\int_0^1\sqrt{\mathbf{r}_{,u}\cdot\mathbf{r}_{,u}}du
	 * \f]
	 * The kernel of the integral above represents the jacobian of the transformation between \f$ u \f$ and \f$ l \f$:
	 * \f[
	 * \begin{eqnarray}
	 * J&=\frac{d l}{d u}=\sqrt{\mathbf{r}_{,u}\cdot\mathbf{r}_{,u}} \\
	 * J_{,u}&=\frac{\mathbf{r}_{,u}\cdot\mathbf{r}_{,uu}}{J}\\
	 * J_{,uu}&=\frac{\mathbf{r}_{,uu}\cdot\mathbf{r}_{,uu}+\mathbf{r}_{,u}\cdot\mathbf{r}_{uuu}-J^2_{,u}}{J}
	 * \end{eqnarray}
	 * \f]
	 * Derivatives with respect to the arc length are obtained by chain rule:
	 * \f[
	 * \begin{eqnarray}
	 * \mathbf{r}_{,l}&=\frac{\mathbf{r}_{,u}}{J}\\
	 * \mathbf{r}_{,ll}&=\frac{\mathbf{r}_{,uu}}{J^2}-J_{,u}\frac{\mathbf{r}_{,u}}{J^3}\\
	 * \mathbf{r}_{,lll}&=\frac{\mathbf{r}_{,uuu}}{J^3}-\frac{\mathbf{r}_{,u}J_{,uu}}{J^4}-3\frac{J_{,u}}{J^2}\left(\frac{\mathbf{r}_{,uu}}{J^2}-J_{,u}\frac{\mathbf{r}_{,u}}{J^3}\right)
	 * \end{eqnarray}
	 * \f]
	 * In matrix form
	 * \f[
	 * \left[
	 * \begin{array}{c}
	 * \mathbf{r}\\
	 * \mathbf{r}_{,u}\\
	 * \mathbf{r}_{,uu}\\
	 * \end{array}\right]=
	 * \left[
	 * \begin{array}{ccc}
	 * 1&0&0\\
	 * 0&J&0\\
	 * 0&J_{,u}&J^2
	 * \end{array}\right]
	 * \left[
	 * \begin{array}{c}
	 * \mathbf{r}\\
	 * \mathbf{r}_{,l}\\
	 * \mathbf{r}_{,ll}\\
	 * \end{array}\right]
	 * \f]
	 * 
	 * Notice that tangent and curvature are orthogonal:
	 * \f[
	 \mathbf{r}_{,l}\cdot\mathbf{r}_{,ll}=\frac{1}{J^3}\left(\mathbf{r}_{u}\cdot\mathbf{r}_{,uu}-J_{,u}\frac{\mathbf{r}_{,u}\cdot\mathbf{r}_{,u}}{J}\right)=
	 \frac{1}{J^3}\left(\mathbf{r}_{u}\cdot\mathbf{r}_{,uu}-\frac{\mathbf{r}_{,u}\cdot\mathbf{r}_{,uu}}{J}\frac{J^2}{J}\right)=0
	 * \f]
	 *
	 * The Frenet-Serret frame (TNB) frame is defined as:
	 * \f[
	 \begin{eqnarray}
	 \hat{\mathbf{t}}&=\mathbf{r}_{,l}&\text{tangent unit vector}\\
	 \hat{\mathbf{n}}&=\frac{\mathbf{r}_{,ll}}{\left|\mathbf{r}_{,ll}\right|}=\frac{\mathbf{r}_{,ll}}{\kappa}&\text{principal unit normal}\\
	 \hat{\mathbf{b}}&=\hat{\mathbf{t}}\times\hat{\mathbf{n}}=\frac{\mathbf{r}_{,u}\times\mathbf{r}_{,uu}}{J^3\kappa}&\text{binormal unit vector}\\
	 \end{eqnarray}
	 * \f]
	 * The curvature
	 * \f[
	 * \kappa=\left|\mathbf{r}_{,ll}\right| = \frac{1}{J^2}\sqrt{\left|\mathbf{r}_{,uu}\right|^2-J_{,u}^2}
	 * \f]
	 * \f[
	 * \begin{eqnarray}
	 * \kappa_{,u}&= -2\frac{J_{,u}}{J^3}\sqrt{\left|\mathbf{r}_{,uu}\right|^2-J_{,u}^2}+\frac{1}{J^2}\frac{\mathbf{r}_{,uu}\cdot\mathbf{r}_{uuu}-J_{,u}J_{,uu}}{\sqrt{\left|\mathbf{r}_{,uu}\right|^2-J_{,u}^2}}\\
	 * &=-2\frac{J_{,u}}{J}\kappa+\frac{\mathbf{r}_{,uu}\cdot\mathbf{r}_{uuu}-J_{,u}J_{,uu}}{J^4\kappa}
	 * \end{eqnarray}
	 * \f]
	 * Torsion measures the deviation from a planar curve:
	 * \f[
	 * \tau=\frac{\mathbf{r}_{,u}\cdot\left(\mathbf{r}_{uu}\times\mathbf{r}_{uuu}\right)}{J^6\kappa^2}
	 * \f]
	 */
	
	

	
	
	
	template <typename Derived, short unsigned int dim>
	class ParametricCurve
    {
		
		typedef Eigen::Matrix<double, dim, 1>		VectorDim;

	private:
		
		//! The CRTP pattern, renamed to avoid the dreaded diamond inheritance problem
		Derived&   derivedCurve() { return *static_cast<Derived*>(this); }
		Derived* p_derivedCurve() { return  static_cast<Derived*>(this); }
		const Derived&   derivedCurve() const { return *static_cast<const Derived*>(this); }
		const Derived* p_derivedCurve() const { return  static_cast<const Derived*>(this); }
		
	public:
		
/*************************************************************/
		double get_j(const double & u) const {
			/*! The  jacobian of the transformation r(u)->r(l):
			 * \f[
			 * J=\frac{d l}{d u}=\sqrt{\mathbf{r}_{,u}\cdot\mathbf{r}_{,u}}
			 * \f]
			 */
			return p_derivedCurve()->get_ru(u).norm();
		}
		
		/*************************************************************/
		double get_jSquare(const double & u) const {
			return p_derivedCurve()->get_ru(u).squaredNorm();
		}

		/*************************************************************/
		VectorDim get_rl(const double & u) const {
			return p_derivedCurve()->get_ru(u).normalized();
		}

		/*************************************************************/
		double get_ju(const double & u) const {
			/*! The  first derivative of the jacobian:
			 * \f[
			 * J_{,u}=\frac{\mathbf{r}_{,u}\cdot\mathbf{r}_{,uu}}{J}
			 * \f]
			 */
			return get_rl(u).dot(p_derivedCurve()->get_ruu(u));
		}

		
		/*************************************************************/
		VectorDim get_rll(const double & u) const {
			/*! The curvature vector at u:
			 * \f[
			 * \mathbf{r}_{,ll}=\frac{\mathbf{r}_{,uu}}{J^2}-J_{,u}\frac{\mathbf{r}_{,u}}{J^3}
			 * \f]
			 */
			const double j2(get_jSquare(u));
			VectorDim ru (p_derivedCurve()->get_ru (u));
			VectorDim ruu(p_derivedCurve()->get_ruu(u));
			return (ruu - ru*ru.dot(ruu)/j2)/j2;
		}
		
		
		/*************************************************************/
		//! the scalar curvature
		double get_kappa(const double & u) const {
			/*! The scalar curvature at u:
			 * \f[
			 * \kappa=\left|\mathbf{r}_{,ll}\right|
			 * \f]
			 */
			return get_rll(u).norm();
		}
		
		
		/*************************************************************/
		template <short unsigned int qOrder, template <short unsigned int, size_t> class QuadratureRule>
		double arcLength() const {
			double L(0.0);
			Quadrature<1,qOrder,QuadratureRule>::integrate(this,L,&Derived::get_j);
			return L;
		}

		/*************************************************************/
		template <short unsigned int qOrder, template <short unsigned int, size_t> class QuadratureRule>
		VectorDim rm() const {
			VectorDim temp(VectorDim::Zero());
			Quadrature<1,qOrder,QuadratureRule>::integrate(this,temp,&Derived::rm_integrand);
			return temp/arcLength<qOrder,QuadratureRule>();
		}
		
		VectorDim rm_integrand(const double& u) const {
			return p_derivedCurve()->get_r(u)*get_j(u);
		}
		
	};
	
} // namespace model
#endif


