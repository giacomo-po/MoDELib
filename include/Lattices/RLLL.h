/* This file is part of PIL, the Particle Interaction Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_RLLL_h_
#define model_RLLL_h_

#include <iostream>
#include <iomanip>
#include <cfloat> // FLT_EPSILON
#include <Eigen/Dense>
//#include <Eigen/QR>
//#include <RoundEigen.h>
//#include <SmithDecomposition.h>
//#include <Lattice.h>
//#include <RationalMatrix.h>
//#include <LLL.h>


namespace model
{
    
    
    class RLLL
    {/* Ying Hung Gan, Cong Ling, Complex Lattice Reduction Algorithm for Low-Complexity Full-Diversity MIMO Detection
      * IEEE TRANSACTIONS ON SIGNAL PROCESSING, VOL. 57, NO. 7, JULY 2009
      */
        
        typedef Eigen::MatrixXd MatrixType;
        typedef Eigen::VectorXd VectorType;
        
        
        MatrixType B;
        Eigen::Matrix<long int,Eigen::Dynamic,Eigen::Dynamic> U;
        
        /**********************************************************************/
        void update(VectorType& H,
                    MatrixType& M,
                    const int& k);
        /**********************************************************************/
        void size_reduce(MatrixType& M,
                         const int& k,
                         const int& j);
    public:
        
        
        /**********************************************************************/
        RLLL(const MatrixType& B0,
             const double& delta) ;
        
        /**********************************************************************/
        const MatrixType& reducedBasis() const;
        /**********************************************************************/
        const Eigen::Matrix<long int,Eigen::Dynamic,Eigen::Dynamic>& unimodularMatrix() const;
        
    };
    
    
} // end namespace
#endif


