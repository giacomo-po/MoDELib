/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_GammaSurface_H_
#define model_GammaSurface_H_

#include <memory>
#include <assert.h>
#include <LatticeModule.h>
//#include <LatticePlaneBase.h>
//#include <LatticeVector.h>
//#include <RationalLatticeDirection.h>
#include <PeriodicLatticeInterpolant.h>

namespace model
{
    
    struct GammaSurface : PeriodicLatticeInterpolant<2>
    {
        
        static constexpr int dim=3;
        static constexpr int lowerDim=dim-1;
        
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef Eigen::Matrix<size_t,lowerDim,1> VectorLowerDimI;
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        typedef Eigen::Matrix<double,lowerDim,lowerDim> MatrixLowerDim;
        
//        const LatticePlaneBase latticePlane;
        const Eigen::Matrix3d G2L;
        
//        GammaSurface(const LatticePlaneBase& n,
//                     const Eigen::Matrix<double,Eigen::Dynamic,lowerDim>& waveVectors,
//                     const Eigen::Matrix<double,Eigen::Dynamic,lowerDim+1>& f,
//                     const int& rotSymm,
//                     const std::vector<Eigen::Matrix<double,lowerDim,1>>& mirSymm);
        
        GammaSurface(const LatticeVector<dim>& a1,
                     const LatticeVector<dim>& a2,
                     const Eigen::Matrix<double,Eigen::Dynamic,lowerDim>& waveVectors,
                     const Eigen::Matrix<double,Eigen::Dynamic,lowerDim+1>& f,
                     const int& rotSymm,
                     const std::vector<Eigen::Matrix<double,lowerDim,1>>& mirSymm);
        
        double operator()(const VectorDim& b);
        static MatrixDim getG2L(const VectorDim& x,const VectorDim& z);
//        static MatrixLowerDim getLocalBasis(const LatticePlaneBase& n);
        static MatrixLowerDim getLocalBasis(const LatticeVector<dim>& a1,const LatticeVector<dim>& a2);

    };
    
}
#endif
