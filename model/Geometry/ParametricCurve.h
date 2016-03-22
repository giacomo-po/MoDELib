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

#include <model/Quadrature/Quadrature.h>
//#include <model/Geometry/distance2line.h>


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
	class ParametricCurve {
		
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
		VectorDim& get_rll(const double & u) const {
			/*! The curvature vector at u:
			 * \f[
			 * \mathbf{r}_{,ll}=\frac{\mathbf{r}_{,uu}}{J^2}-J_{,u}\frac{\mathbf{r}_{,u}}{J^3}
			 * \f]
			 */
			double j2(get_jSquare(u));
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
	
	//////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////
} // namespace model
#endif



//		double arcLength_integrand(const double& u) const {
//			return get_j(u);
//		}





//		//////////////////////////////////////////////////////////////////////////
//		double distance2line(const VectorDim & lineDir, const VectorDim & lineOrg, const double & u_in){
//			//! CHANGE THIS WITH THE USE OF EIGEN3 BUIT-IN GEOMETRY FUNCTIONALITIES
//			return distance2line<dim>(lineDir,lineOrg,p_derivedCurve()->get_r(u_in));
//		}



/*	
 * A parametric curve is a function r(u): [0,1] -> R^d that described a curve in a d-dimensional
 * space. The function r(u) is implemented as the pure virtual make_r(u) and must be defined in
 * the derived classes. Given r(u) all the other geometric properties are computed here. 
 *
 * Notice that the first and second parametric derivatives make_ru() and make_ruu() are implemented 
 * by numerical differentiation but they are also declared as virtual so that more efficient and 
 * precise implementation can be included in derived classes (eg splines).
 *
 * The class is equipped with a Gauss Quadrature class to perform integration along the curve.
 *
 * The optimized calculation tree is shown below:
 *
 * \verbatim
 * 
 * make_ruu()       make_ru()      make_r()        (CRTP level: implemented in derived class)
 *     ^                ^              ^
 *     |                |              |				
 *     |            make_j()           |
 *     |                ^              |
 *     |                |              |
 *     |            make_rl()          |			
 *     |                ^              |
 *     |                |              |              
 *     +----------------+
 *     |
 * make_ju()
 *     ^
 *     |
 * make_rll()
 *     ^
 *     |
 * make_kappa()
 *
 * \endverbatim
 */







//	template <typename Derived, short unsigned int dim, short unsigned int qOrder,  template <short unsigned int, short unsigned int> class Rule>

//		typedef Eigen::Matrix<double, 1, qOrder>	VectorQuadrature;
//		typedef Eigen::Matrix<double, dim, qOrder>	MatrixDimQuadrature;


//const double& get_j(){return j;}

//////////////////////////////////////////////////////////////
// rl = ru/j
//		void make_rl(const double & u){
//			return p_derivedCurve()->make_ru(u).normalized();
//		}

//		const VectorDim& get_rl(const int & k){
//			return rlgauss.col(k);
//		}



//	double& get_ju(){return ju;}

//		//////////////////////////////////////////////////////////////
//		//! rll = (ruu-rl*ju)/j^2
//		void make_rll(const double & u){
//			
//			make_ju(u);
//			//	std::cout<<"ParametricCurve::make_rll"<<std::endl;
//			rll=(ruu-rl*ju)/std::pow(j,2);
//		}

//		const VectorDim& get_rll(const int & k) const {
//			double j2=get_jSquare();
//			return () / get;
//		}



//		
//		const double& get_kappa(const int & k) const {
//			return kappagauss(k);
//		}
//		
//		const double& get_kappa(const double & u) {
//			make_kappa(u); 
//			return kappa;
//		}
//	double& get_kappa(){return kappa;}


//VectorDim& get_rl(){return rl;}

//////////////////////////////////////////////////////////////
//! ruu(u) can be redefined in derived classes to take advantage of particular cases (eg. splines)
// CAREFUL, this is the derivative with respect to u in [0,1] not with respect to the parametrization !!!
//	virtual void make_ruu(const double & u){
//	std::cout<<"ParametricCurve::make_ruu"<<std::endl;
//		ru=(get_ru(u+du)-get_ru(u-du))/(2.0*du);
//	}

//		const VectorDim& get_ruu(const double & u){
//			p_derivedCurve()->make_ruu(u); 
//			return ruu;
//		}
//	VectorDim& get_ruu(){return ruu;}

//////////////////////////////////////////////////////////////
//! ju=rl*ruu
//		void make_ju(const double & u){
//			make_rl(u);
//			p_derivedCurve()->make_ruu(u);
//			//	std::cout<<"ParametricCurve::make_ju"<<std::endl;
//			ju=rl.dot(ruu);
//		}







//		//////////////////////////////////////////////////////////////
//		//! Compute all the parameters that need quadrature in an optimized way
//		void update(){
//			
//			
//			
//			rgauss.setZero();
//			rugauss.setZero();
//			ruugauss.setZero();
//			
//			rlgauss.setZero();
//			rllgauss.setZero();
//			
//			jgauss.setZero();
//			kappagauss.setZero();
//			
//			
//			//! Initialize average values
//			arcL=0.0;
//			kappam=0.0;
//			rm.setZero();
//			rlm.setZero();
//			rllm.setZero();
//			
//			
//			//double jw;
//			
//			for(int k=0;k<qOrder;++k){
//				make_all(k); // Calculates N, r, j, rl, rll for current value of u=abscissa(k)
//				double jw=j*this->weight(k);
//				
//				arcL+=jw;
//				rm+=r*jw;
//				rlm+=rl*jw;
//				rllm+=rll*jw;
//				kappam+=kappa*jw;
//			}
//			
//			rm/=arcL;
//			rlm/=arcL;	
//			rllm/=arcL;	
//			kappam/=arcL;
//			
//		}

//		double& get_arcL(){return arcL;}	
//		
//		// return rgauss
//		const MatrixDimQuadrature&	get_rgauss() const {return rgauss;}
//		double						get_rgauss(const size_t & i, const size_t & j) const {return rgauss(i,j);}
//		
//		// return rugauss
//		const MatrixDimQuadrature&	get_rugauss() const {return rugauss;}
//		double						get_rugauss(const size_t & i, const size_t & j) const {return rugauss(i,j);}
//		
//		// return rlgauss
//		const MatrixDimQuadrature&	get_rlgauss(){return rlgauss;}
//		
//		// return jgauss
//		const VectorQuadrature&		get_jgauss(){return jgauss;}
//		double						get_jgauss(const size_t & i) const {return jgauss(i);}
//		
//		// return rllgauss
//		const MatrixDimQuadrature& get_rllgauss(){return rllgauss;}
//		
//		
//		
//		const VectorDim& get_rm(){return rm;}
//		const VectorDim& get_rlm(){return rlm;}
//		const VectorDim& get_rllm(){return rllm;}
//		
//		double& get_kappam(){return kappam;}



//		//////////////////////////////////////////////////////////////////////////
//		double distance2line(const VectorDim & lineDir, const VectorDim & lineOrg, const int & k){
//			return distance2line<dim>(lineDir,lineOrg,rgauss.col(k));
//		}

//	private:
//		
//		double areaBetweenLineKernel(const double& u) const{
//		
//		}
//		
//	public:
//////////////////////////////////////////////////////////////////////////
//		template<short unsigned int quadratureOrder>
//		double areaBetweenLine(VectorDim lineDir, const VectorDim & lineOrg){
//			lineDir.normalize();
//			double A=0.0;
//			
//			for(int k=0;k<qOrder;++k){
//				A+=distance2line(lineDir,lineOrg,k)*rugauss.col(k).dot(lineDir)*this->weight(k);
//			}
//			return A;
//		}



//		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
//		
//		//////////////////////////////////////////////////////////////
//		// Default Constructor
//		ParametricCurve(){
//			
//			
//			
//			//	du=0.00001;
//			
//			r.setZero();
//			ru.setZero();
//			ruu.setZero();
//			rl.setZero();
//			rll.setZero();
//			
//			
////			rgauss.setZero();
////			
////			rlgauss.setZero();
////			rllgauss.setZero();
////			
////			jgauss.setZero();
////			kappagauss.setZero();
////			
////			
////			//! Initialize average values
////			arcL=0.0;
////			kappam=0.0;
////			rm.setZero();
////			rlm.setZero();
////			rllm.setZero();
//			
//		}

//		//////////////////////////////////////////////////////////////
//		//! Calculate all geometric parameters @ u in an optimized way
//		void make_all(const double & u){
//			p_derivedCurve()->make_r(u);		// Updates  r	
//			make_kappa(u);	// Updates (in order) ruu, ru, j, ju, rll, kappa
//		} 

//		//////////////////////////////////////////////////////////////
//		void make_all(const int & k){
//			
//			p_derivedCurve()->make_r(k);				// Updates  r	
//			make_kappa(this->abscissa(k));	            // Updates (in order) ruu, ru, j, ju, rll, kappa
//			
//			rgauss.col(k)=r;
//			rugauss.col(k)=ru;
//			ruugauss.col(k)=ruu;
//			rlgauss.col(k)=rl;
//			rllgauss.col(k)=rll;
//			jgauss(k)=j;
//			kappagauss(k)=kappa;
//		} 

//		//////////////////////////////////////////////////////////////
//		//! Overwrite GaussPointsD::set_QuadratureOrder
//		void set_QuadratureOrder(int goin){
//			//! If the quadrature order changes make sure to update the quantities at the gauss points
//			GaussPoints1D::set_QuadratureOrder(goin);
//			p_derivedCurve()->update();
//		}




//! r(u): pure virtual function that must be defined in derived classes
//		virtual void make_r(const double & u) = 0 ;	
//		virtual void make_r(const int & k) = 0 ;	
//	const VectorDim& get_r(const int & k){return rgauss.col(k);}
//	const VectorDim& get_r(const double & u){p_derivedCurve()->make_r(u); return r;}
//VectorDim& get_r(){return r;}

//////////////////////////////////////////////////////////////
//! ru(u): can be redefined in derived classes to take advantage of particular cases (eg. splines)
// CAREFUL, this is the derivative with respect to u in [0,1] not with respect to the parametrization !!!
//	virtual void make_ru(const double & u){
//	std::cout<<"ParametricCurve::make_ru"<<std::endl;
//		ru=(get_r(u+du)-get_r(u-du))/(2.0*du);
//	}

//	const VectorDim& get_ru(const double & u){p_derivedCurve()->make_ru(u); return ru;}		
//VectorDim& get_ru(){return ru;}

//		//////////////////////////////////////////////////////////////
//		//! j=norm(ru*ru)
//		void make_j(const double & u){
//			p_derivedCurve()->make_ru(u);
//			//	std::cout<<"ParametricCurve::make_j"<<std::endl;	
//			j=ru.norm();
//		}


//	protected:
//! Position at u
//		VectorDim r;

//! First parametric derivative dr/du
//		VectorDim ru;

//! Second parametric derivative d^2r/du^2
//		VectorDim ruu;

//! Arc length of the curve
//double arcL;

//! Average values of postition, tangent and curvature along the curve
//VectorDim rm, rlm, rllm;

//! Average value of the norm of the curvature along the curve
//double kappam;

//! Positions corrersponding to the quadrature points
//MatrixDimQuadrature rgauss;

//! Positions corrersponding to the quadrature points
//MatrixDimQuadrature rugauss;

//! Positions corrersponding to the quadrature points
//MatrixDimQuadrature ruugauss;

//! Tangents corrersponding to the quadrature points
//MatrixDimQuadrature rlgauss;

//! Curvature corrersponding to the quadrature points
//MatrixDimQuadrature rllgauss;

//! Scalar jacobian corrersponding to the quadrature points
//VectorQuadrature jgauss;

//! Scalar jacobian corrersponding to the quadrature points
//VectorQuadrature kappagauss;

//	private:

//! offset used for numerical differentiation
//	double du;

//! Unit tangent at u
//		VectorDim rl;	

//! Curvature at u
//		VectorDim rll;





//! Jacobian of the transformation dl->du : j=dl/du=sqrt(ru*ru)
//		double j;

//! Parametric derivative of the Jacobian: ju=dj/du
//		double ju;

//! Scalar curvature
//		double kappa;



//double get_Area(fglfsd)
/*
 //////////////////////////////////////////////////////////////
 // Operators
 template <typename SomeDerived, int SomeDim>	// Notice that the typename is different!
 friend std::ostream& operator<<(std::ostream & os,  ParametricCurve<SomeDerived, SomeDim> & PC){
 
 os << "///////////////////////////////////////////" << std::endl;
 os << "ParametricCurve Object:" << std::endl;
 os << "position @ u =" <<std::endl<< PC.get_r() << std::endl;
 os << "geometric (unit) tangent @ u =" <<std::endl<< PC.get_rl() << std::endl;
 os << "geometric curvature vector @ u =" <<std::endl<< PC.get_rll() << std::endl;
 os << "geometric curvature @ u =" << PC.get_kappa() << std::endl;
 os << "parametric first derivative of r @ u =" <<std::endl<< PC.get_ru() << std::endl;
 os << "parametric second derivative of r @ u =" <<std::endl<< PC.get_ruu() << std::endl;
 os << "jacobian j=dl/du @ u =" << PC.get_j() << std::endl;
 os << "derivative of the jacobian ju=dj/du @ u =" << PC.get_ju() << std::endl;
 return os;
 } */
