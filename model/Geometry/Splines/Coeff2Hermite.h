/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

/**********************************************************************************************************************/
/**********************************************************************************************************************/
#ifndef model_COEFF2HERMITE_H_
#define model_COEFF2HERMITE_H_

#include <Eigen/Dense>
#include <CTM.h>
#include <PermutationWithoutRepetition.h>

namespace model {
	
	
	/* FillC2H0 *******************************************************************************/
	// case i,j
	template<short unsigned int P,short unsigned int i, short unsigned int j>
	struct FillC2H0 : public FillC2H0<P,i,j-1>{		
		FillC2H0(Eigen::Matrix<double,(P+1)/2,(P+1)>& C2H0) : FillC2H0<P,i,j-1>::FillC2H0(C2H0){
				C2H0(i,j) = 0.0;
		}
	};
	
	// case i,i
	template<short unsigned int P,short unsigned int i>
	struct FillC2H0<P,i,i> : public FillC2H0<P,i,i-1>{
		enum {j=i};
		FillC2H0(Eigen::Matrix<double,(P+1)/2,(P+1)>& C2H0) : FillC2H0<P,i,j-1>::FillC2H0(C2H0){
				C2H0(i,j) = CTM::factorial(i);
		}
	};
	
	// case i,0, i!=0
	template<short unsigned int P,short unsigned int i>
	struct FillC2H0<P,i,0> : public FillC2H0<P,i-1,P>{
		enum {j=0};		
		FillC2H0(Eigen::Matrix<double,(P+1)/2,(P+1)>& C2H0) : FillC2H0<P,i-1,P>::FillC2H0(C2H0){
				C2H0(i,j) = 0.0;
		}
	};
	
	// case i=j=0
	template<short unsigned int P>
	struct FillC2H0<P,0,0>{
		enum {i=0};
		enum {j=0};		
		FillC2H0(Eigen::Matrix<double,(P+1)/2,(P+1)>& C2H0){
				C2H0(i,j) = CTM::factorial(i);
		}
	};
	
	/* FillC2H1 *******************************************************************************/
	// case i,j
	template<short unsigned int P,short unsigned int i, short unsigned int j>
	struct FillC2H1 : public FillC2H1<P,i,j-1>{		
		FillC2H1(Eigen::Matrix<double,(P+1)/2,(P+1)>& C2H1) : FillC2H1<P,i,j-1>::FillC2H1(C2H1){
			C2H1(i,j) = PermutationWithoutRepetition<i>::value(j)*(j>=i);
		}
	};
	
	// case i,0, i!=0
	template<short unsigned int P,short unsigned int i>
	struct FillC2H1<P,i,0> : public FillC2H1<P,i-1,P>{
		enum {j=0};		
		FillC2H1(Eigen::Matrix<double,(P+1)/2,(P+1)>& C2H1) : FillC2H1<P,i-1,P>::FillC2H1(C2H1){
			C2H1(i,j) = PermutationWithoutRepetition<i>::value(j)*(j>=i);
		}
	};
	
	// case i=j=0
	template<short unsigned int P>
	struct FillC2H1<P,0,0>{
		enum {i=0, j=0};	
		FillC2H1(Eigen::Matrix<double,(P+1)/2,(P+1)>& C2H1){
			C2H1(i,j) = PermutationWithoutRepetition<i>::value(j)*(j>=i);
		}
	};
	

	/**********************************************************************************************************************/
	/* Coeff2Hermite<polyDegree> class template ***************************************************************************/
	/**********************************************************************************************************************/
	/*! \brief Class template that performs the transformation between the polynomial coefficient matrix 
	 *  \f$ \mathbf{C}=\left[\mathbf{c}_0\ \mathbf{c}_1\ \ldots \ \mathbf{c}_P\right]\f$ and the Hermite coefficient matrix 
	 *  \f$ \mathbf{H}=\left[\mathbf{r}(0),\ \ldots \frac{d^\frac{P-1}{2}\mathbf{r}}{dt^\frac{P-1}{2}}(0),\ \mathbf{r}(1),\ \ldots \frac{d^\frac{P-1}{2}\mathbf{r}}{dt^\frac{P-1}{2}}(1)\right]\f$ 
	 *  of a spline \f$\mathbf{r}(u)=\sum_{i=0}^P\mathbf{c}_it^i\f$ of odd degree \f$P\f$, parametrized in \f$[0\ 1]\f$. 
	 *  
	 *  Consider the spline:
	 *  \[
	 *  \mathbf{r}(u)=\sum_{i=0}^P\mathbf{c}_it^i
	 *  \]
	 *  where \f$t\in[0,1]\f$ and P is odd. The matrix
	 *  \f[
	 *  \left[\mathbf{c}_0\ \mathbf{c}_1\ \ldots \ \mathbf{c}_P\right]
	 *  \f]
	 *  is the coefficient matrix of the spline, while the matrix
	 *  \f[
	 *  \mathbf{H}=\left[\mathbf{r}(0),\ \ldots \frac{d^\frac{P-1}{2}\mathbf{r}}{dt^\frac{P-1}{2}}(0),\ \mathbf{r}(1),\ \ldots \frac{d^\frac{P-1}{2}\mathbf{r}}{dt^\frac{P-1}{2}}(1)\right]
	 *  \f]
	 *  is the Hermite matrix of the spline. We now seek a relationship between the two.
	 *  \f[
	 *  \frac{d^n\mathbf{r}}{dt^n}=\sum_{i=0}^P \underbrace{i(i-1)\ldots(i-n+1)}_{n\ terms} \mathbf{c}_it^{i-n}
	 *  \f]
	 *  Since only terms having \f$i-n+1>0\f$ survive in the sum, that is \f$i\ge n\f$ then we obtain:
	 *  \f[
	 *  \frac{d^n\mathbf{r}}{dt^n}=\sum_{i=n}^P \frac{i!}{(i-n)!} \mathbf{c}_it^{i-n}
	 *  \f]
	 *  Letting \f$t=0\f$ and \f$t=1\f$, one finds that the Hermite coefficients of order n (n=0...(P-1)/2) are:
	 *  \f[
	 *  \frac{d^n\mathbf{r}}{dt^n}(0)=\sum_{i=n}^P \frac{i!}{(i-n)!} \mathbf{c}_i0^{i-n}=n!\mathbf{c}_n
	 *  \f]
	 *  \f[
	 *  \frac{d^n\mathbf{r}}{dt^n}(1)=\sum_{i=n}^P \frac{i!}{(i-n)!} \mathbf{c}_i1^{i-n}=\sum_{i=n}^P \frac{i!}{(i-n)!} \mathbf{c}_i
	 *  \f]
	 *  In matrix form this reads:
	 *  \f[
     *  \underbrace{\left[\begin{array}{l}
	 *	{d^0\mathbf{r}}/{dt^0}(0)\\
	 *	{d^1\mathbf{r}}/{dt^1}(0)\\
	 *	{d^2\mathbf{r}}/{dt^2}(0)\\
	 *	\vdots\\
	 *	{d^0\mathbf{r}}/{dt^0}(1)\\
	 *	{d^1\mathbf{r}}/{dt^1}(1)\\
	 *	{d^2\mathbf{r}}/{dt^2}(1)\\
	 *	\vdots\\
	 *	\end{array}\right]}_{\mathbf{H}^T}=
	 *	\left[\begin{array}{cccc}
	 *	0!& 0 & 0 & \ldots\\
	 *	0 & 1! & 0 & \ldots\\
	 *	0 & 0 & 2! & \ldots\\
	 *	\vdots & \vdots & \vdots & \ddots\\	 
	 *	\frac{0!}{(0-0)!}& \frac{1!}{(1-0)!} & \frac{2!}{(2-0)!} & \ldots\\
	 *	0& \frac{1!}{(1-1)!} & \frac{2!}{(2-1)!} & \ldots\\
	 *	0& 0 & \frac{2!}{(2-2)!} & \ldots\\
	 *	\vdots & \vdots & \vdots & \ddots\\	 
	 *	\end{array}\right]
	 *  \underbrace{\left[\begin{array}{l}
	 *	\mathbf{c}_0\\
	 *	\mathbf{c}_1\\
	 *	\mathbf{c}_2\\	 
	 *	\vdots\\
	 *	\end{array}\right]}_{\mathbf{C}^T}
	 *  \f]
	 *  For a linear spline (P=1, n=0) we find:
	 *  \f[
	 *  \mathbf{H}=\left[\mathbf{c}_0,\  \mathbf{c}_0+\mathbf{c}_1\right]
	 *  \f]
	 *  \f[
	 *  \mathbf{C}=\left[\mathbf{h}_0,\  \mathbf{h}_1-\mathbf{h}_0\right]
	 *  \f]
	 *  For a cubic spline (P=3, n=0,1) we find:
	 *  \f[
	 *  \mathbf{H}=\left[\mathbf{c}_0,\ \mathbf{c}_1,\  \mathbf{c}_0+\mathbf{c}_1+\mathbf{c}_2+\mathbf{c}_3,\ \mathbf{c}_1+2\mathbf{c}_2+3\mathbf{c}_3\right]
	 *  \f]
	 *  \f[
	 *  \mathbf{C}=\left[\mathbf{h}_0,\ \mathbf{h}_1,\  -3\mathbf{h}_0-2\mathbf{h}_1+3\mathbf{h}_2-\mathbf{h}_3,\ 2\mathbf{h}_0+\mathbf{h}_1-2\mathbf{h}_2+\mathbf{h}_3\right]
	 *  \f]
	 */

	template<short unsigned int polyDegree>
	class Coeff2Hermite{
				
		enum {polyCoeff=polyDegree+1};
		
		static Eigen::Matrix<double,polyCoeff,polyCoeff> get_C2H(){
			Eigen::Matrix<double,polyCoeff/2,polyCoeff> temp0;
			FillC2H0<polyDegree,(polyDegree-1)/2,polyDegree> FH0(temp0);

			Eigen::Matrix<double,polyCoeff/2,polyCoeff> temp1;
			FillC2H1<polyDegree,(polyDegree-1)/2,polyDegree> FH1(temp1);
			
			Eigen::Matrix<double,polyCoeff,polyCoeff> temp;
			temp<<temp0,temp1;
			return temp;
		}
		
	public:

		static const Eigen::Matrix<double,polyDegree+1,polyDegree+1> C2H;
		static const Eigen::Matrix<double,polyDegree+1,polyDegree+1> C2HT;

		static const Eigen::Matrix<double,polyDegree+1,polyDegree+1> H2C;
		static const Eigen::Matrix<double,polyDegree+1,polyDegree+1> H2CT;

		
		/* c2h ***********************************************************/
		template<short unsigned int dim>
		static Eigen::Matrix<double,dim,polyCoeff> c2h(const Eigen::Matrix<double,dim,polyCoeff>& C){
			/*! The Hermite matrix corresponding to the polynomial coefficients C (coefficients in columns)
			 *
			 */
			return C*C2HT;
		}
		

		/* h2c ***********************************************************/
		template<short unsigned int dim>
		static Eigen::Matrix<double,dim,polyCoeff> h2c(const Eigen::Matrix<double,dim,polyCoeff>& H){
			/*! The Hermite matrix corresponding to the polynomial coefficients C (coefficients in columns)
			 *
			 */
			return H*H2CT;
		}
	
	};

	/*************************************************************/
	// Static const data members are computed only once at runtime
	template<short unsigned int polyDegree>
	const Eigen::Matrix<double,polyDegree+1,polyDegree+1> Coeff2Hermite<polyDegree>::C2H =Coeff2Hermite<polyDegree>::get_C2H();

	template<short unsigned int polyDegree>
	const Eigen::Matrix<double,polyDegree+1,polyDegree+1> Coeff2Hermite<polyDegree>::C2HT=Coeff2Hermite<polyDegree>::get_C2H().transpose();

	template<short unsigned int polyDegree>
	const Eigen::Matrix<double,polyDegree+1,polyDegree+1> Coeff2Hermite<polyDegree>::H2C =Coeff2Hermite<polyDegree>::get_C2H().inverse();

	template<short unsigned int polyDegree>
	const Eigen::Matrix<double,polyDegree+1,polyDegree+1> Coeff2Hermite<polyDegree>::H2CT=Coeff2Hermite<polyDegree>::get_C2H().inverse().transpose();	
	
}// namespace model
#endif
/**********************************************************************************************************************/
/**********************************************************************************************************************/
