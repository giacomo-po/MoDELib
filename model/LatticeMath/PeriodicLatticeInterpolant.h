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


namespace model
{
    
    template <int dim>
    struct WaveVectorsAssembler
    {
        typedef Eigen::Matrix<size_t,dim,1> VectorDimI;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        
        
        static Eigen::Matrix<double,Eigen::Dynamic,dim> get(const VectorDimI& N,
                                                            const VectorDim& D)
        {
            const size_t rows((N.array()+1).prod());
            const auto lowerBlock(WaveVectorsAssembler<dim-1>::get(N.template segment<dim-1>(1),D.template segment<dim-1>(1)));
            Eigen::Matrix<double,Eigen::Dynamic,dim> temp(rows,dim);
            for(size_t k=0;k<N(0)+1;++k)
            {
//                temp.block(k*lowerBlock.rows(),0,lowerBlock.rows(),1)=MatrixXd::Constant(owerBlock.rows(), 1, k/D(0));
                temp.block(k*lowerBlock.rows(),0,lowerBlock.rows(),dim)<<Eigen::MatrixXd::Constant(lowerBlock.rows(), 1, k/D(0)),lowerBlock;
            }
            return temp;
        }

    };


    template<>
    struct WaveVectorsAssembler<1>
    {
        static constexpr int dim=1;
        typedef Eigen::Matrix<size_t,dim,1> VectorDimI;
        typedef Eigen::Matrix<double,dim,1> VectorDim;

        
        static Eigen::Matrix<double,Eigen::Dynamic,dim> get(const VectorDimI& N,
                                                            const VectorDim& D)
        {
            const size_t rows((N.array()+1).prod());
            Eigen::Matrix<double,Eigen::Dynamic,dim> temp(rows,dim);
            for(size_t k=0;k<N(0)+1;++k)
            {
                temp(k,0)=k/D(0);
            }
            return temp;
        }
        
    };

    
    template <int dim>
    class PeriodicLatticeInterpolant
    {
        
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef Eigen::Matrix<size_t,dim,1> VectorDimI;
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;

        
        const MatrixDim A;   // column matrix of lattice vectors
        const MatrixDim B;   // column matrix of reciprocal lattice vectors
        const VectorDimI N; // numerators of k-point indices
        const VectorDim D;   // denumerators of k-point indices
        
        const Eigen::Matrix<double,Eigen::Dynamic,dim> waveVectors;
        
//        std::vector<VectorDim> waveVectors;
        
    public:
        
        PeriodicLatticeInterpolant(const MatrixDim& A_in,
                                   const VectorDimI& nums_in,
                                   const VectorDimI& dens_in) :
        /* init */ A(A_in)
        /* init */,B(2.0*M_PI*A.inverse().transpose())
        /* init */,N(nums_in)
        /* init */,D(dens_in.template cast<double>())
        /* init */,waveVectors(WaveVectorsAssembler<dim>::get(N,D))
        {
            
            std::cout<<"waveVectors=\n"<<waveVectors<<std::endl;
            
        }
        
    };
    
} // end namespace
#endif

