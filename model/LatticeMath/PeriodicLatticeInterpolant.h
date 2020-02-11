/* This file is part of PIL, the Particle Interaction Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PeriodicLatticeInterpolant_h_
#define model_PeriodicLatticeInterpolant_h_

#include <Eigen/Dense>
#include <vector>
#include <array>

template<size_t N>
std::vector<std::array<int N>> nestedLoopIndices(const Eigen::Matrix<int,N,1>& starts,
                                                 const Eigen::Matrix<int,N,1>& ends)
{
    for(int i=starts(0);i<ends(0);++i)
    {
        nestedLoopIndices(starts.segment<N-1>(1),end.segment<N-1>(1));
    }
        
}

template<size_t N>
struct NestedForLoop
{
    
    template<int...Is,typename T...args>
    void apply(const int& iMin,const int& iMax )
    {
        for(int i=iMin;i<iMax;++i)
        {
            NestedForLoop.apply(Is...,args...)
        }
    }
    
};

namespace model
{
    template <int dim>
    class PeriodicLatticeInterpolant
    {
        
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef Eigen::Matrix<int,dim,1> VectorDimI;
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;

        
        const MatrixDim A; // column matrix of lattice vectors
        const MatrixDim B; // column matrix of reciprocal lattice vectors
        const VectorDim  N; // column matrix of reciprocal lattice vectors
        const VectorDimI D; // column matrix of reciprocal lattice vectors
        
        std::vector<VectorDim> waveVectors;
        
    public:
        PeriodicLatticeInterpolant(const MatrixDim& A_in,
                                   const VectorDimI& nums_in,
                                   const VectorDimI& dens_in) :
        /* init */ A(A_in)
        /* init */,B(A.inverse().transpose())
        /* init */,N(nums_in.template cast<double>())
        /* init */,D(dens_in)
        {
            
        }
        
    };
    
} // end namespace
#endif

