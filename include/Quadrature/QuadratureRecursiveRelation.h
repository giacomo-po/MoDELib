/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_QUADRATURERECURSIVERELATION_H_
#define model_QUADRATURERECURSIVERELATION_H_

#include <iostream>
#include <iomanip>
#include <assert.h>
#include <math.h>
#include <set>
#include <utility>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>


namespace model
{
    
    /*! \brief A static class template for generating quadrature abcsissas and
     * weights for arbitrary-order and rule in dimension 1.
     *
     * For a given quadrature rule (eg. Gauss-Legendre, Gauss-Chebyshev...),
     * abscissas and weights can be coputed from the coefficients of the recursive
     * relation that define the orthogonal polynomials of that rule.
     * Recursive relations are in the form:
     *	\f[
     *		P_n(x) = (A_n*x+B_n) * P_{n-1} - C_n * P_{n-2}
     *	\f]
     *	The actual coefficients \f$A_n\f$, \f$B_n\f$ and \f$C_n\f$
     *
     *	References:
     *	- [1] Golub, G.H., and Welsch, J.H. 1969, "Calculation of Gauss Quadrature Rules",
     *		Mathematics of Computation, vol. 23, pp. 221–230.
     *
     *	- [2] Press, W. et al. "Numerical Recipes", 3rd edition, p.188.
     */
    struct QuadratureRecursiveRelation
    {
        
        enum {dim=1};
        
        /* abcsissasAndWeights *******************************/
        template<typename PolynomialType>
        static Eigen::MatrixXd abcsissasAndWeights(const int& qOrder)
        {
            std::cout<<"computing qOrder="<<qOrder<<"..."<<std::flush;

            // Create and fill the Jacobi upper triangular matrix
            Eigen::MatrixXd jUT = Eigen::MatrixXd::Zero(qOrder,qOrder);
            for (int n=0;n<qOrder;++n)
            {
                int k = n+1;
                jUT(n,n) = -PolynomialType::b(k)/PolynomialType::a(k);
                if (n<qOrder-1)
                {
                    jUT(n,n+1)= std::sqrt(PolynomialType::c(k+1)/PolynomialType::a(k)/PolynomialType::a(k+1));
                }
            }
            
            // Create a SelfAdjointEigenSolver from the self-adjoint part of jUT
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(jUT.template selfadjointView<Eigen::Upper>());
            Eigen::VectorXd		eigenVal = es.eigenvalues();
            Eigen::MatrixXd	eigenVec = es.eigenvectors();
            
            // Insert each abscissa-weight pair in std::set for automatic sorting with increasing abscissas
            typedef std::set<std::pair<double,double> > EigenSetType;
            typedef EigenSetType::const_iterator		EigenSetIterType;
            EigenSetType eigenSet;
            for (int k=0;k<qOrder;++k)
            {
                eigenSet.insert( std::make_pair( (eigenVal(k)+1.0)*0.5, std::pow(eigenVec(0,k),2) ) );
            }
            
            // Copy the sorted values to abscissas_ and weigts_
            Eigen::MatrixXd abcsissasAndWeights_(Eigen::MatrixXd::Zero(dim+1,qOrder));
            int k=0;
            for (EigenSetIterType iter=eigenSet.begin(); iter!=eigenSet.end();++iter)
            {
                abcsissasAndWeights_(0,k)= iter->first;
                abcsissasAndWeights_(1,k)= iter->second;
                ++k;
            }
            
            std::cout<<" done."<<std::endl;

            
            return abcsissasAndWeights_;
        }
        
    };
    
} // namespace model
#endif

