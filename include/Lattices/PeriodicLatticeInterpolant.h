/* This file is part of PIL, the Particle Interaction Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PeriodicLatticeInterpolant_h_
#define model_PeriodicLatticeInterpolant_h_

#include <vector>
#include <array>
#include <cfloat>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <iostream>
namespace model
{

    template <int dim>
    struct PeriodicLatticeInterpolant
    {
        
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef Eigen::Matrix<size_t,dim,1> VectorDimI;
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        
        const MatrixDim A;   // column matrix of lattice vectors
        const MatrixDim B;   // column matrix of reciprocal lattice vectors
        const Eigen::Matrix<double,Eigen::Dynamic,dim> waveVectors;
        const std::vector<MatrixDim> symMatrices;
        const Eigen::Matrix<double,Eigen::Dynamic,2> sinCosCoeffs;
        
        /**********************************************************************/
        static std::vector<MatrixDim> getSymMatrices(const int& rotSymm,
                                                     const std::vector<Eigen::Matrix<double,dim,1>>& mirSymm);
        
        /**********************************************************************/
        Eigen::MatrixXd getSinCosCoeffs(const Eigen::Matrix<double,Eigen::Dynamic,dim+1>& f) const;
        
        /**********************************************************************/
        std::pair<double,double> sinKXcosKX(const int& r,const VectorDim& x) const;

        /**********************************************************************/
        PeriodicLatticeInterpolant(const MatrixDim& A_in,
                                   const Eigen::Matrix<double,Eigen::Dynamic,dim>& waveVectors_in,
                                   const Eigen::Matrix<double,Eigen::Dynamic,dim+1>& f,
                                   const int& rotSymm,
                                   const std::vector<Eigen::Matrix<double,dim,1>>& N) ;
        /**********************************************************************/
        double operator()(const VectorDim& x) const;

    };
    
} // end namespace
#endif
