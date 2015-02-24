/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_GAUSSLEGENDRERECURSIVERELATION_H_
#define model_GAUSSLEGENDRERECURSIVERELATION_H_


#include <assert.h>

namespace model {

	/**************************************************/
	/* GaussLegendreRecursiveRelation:  */
	/**************************************************/
	/*! \brief The
	 *	P_n(x) = (2n-1)/n * x * P_(n-1) - (n-1)/n * P_(n-2)
	 *
	 *	gives:
	 *	A_n = (2*n-1)/n
	 *	B_n = 0
	 *	C_n = (n-1)/n
	 *
	 *	
	 *	
	 *	References:
	 *	[1] Golub, G.H., and Welsch, J.H. 1969, "Calculation of Gauss Quadrature Rules", 
	 *		Mathematics of Computation, vol. 23, pp. 221–230.
	 *
	 *	[2] Press, W. et al. "Numerical Recipes", 3rd edition, p.188.
	 *	
	 */
//	template<short unsigned int qOrder>
	struct GaussLegendreRecursiveRelation
    {
		
		/*! \brief Class template defining the coefficients A, B and C 
		 *  of the recursive relation for orthogonal Legendre polinomials:
		 *	P_n(x) = (A_n*x+B_n) * P_(n-1) - C_n * P_(n-2)
		 */
		
		
		//! The coefficient a
		static double a(const int& n){
			return (2.0*n-1.0)/n;
		}
		
		//! The coefficient b
		static double b(const int& n){
			return 0.0*n;
		}
		
		//! The coefficient c 
		static double c(const int& n){
			return (n-1.0)/n;
		}
		
	};
	
		
	
	//////////////////////////////////////////////////////////////
} // namespace model
#endif

