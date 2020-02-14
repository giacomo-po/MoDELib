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
            const size_t rows(N.prod());
            const auto lowerBlock(WaveVectorsAssembler<dim-1>::get(N.template segment<dim-1>(1),D.template segment<dim-1>(1)));
            Eigen::Matrix<double,Eigen::Dynamic,dim> temp(rows,dim);
            for(size_t k=0;k<N(0);++k)
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
            const size_t rows(N.prod());
            Eigen::Matrix<double,Eigen::Dynamic,dim> temp(rows,dim);
            for(size_t k=0;k<N(0);++k)
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
        Eigen::Matrix<double,Eigen::Dynamic,2> sinCosCoeffs;

//        std::vector<VectorDim> waveVectors;
        
    public:
        
        PeriodicLatticeInterpolant(const MatrixDim& A_in,
                                   const VectorDimI& nums_in,
                                   const VectorDimI& dens_in) :
        /* init */ A(A_in)
        /* init */,B(2.0*M_PI*A.inverse().transpose())
        /* init */,N(nums_in)
        /* init */,D(dens_in.template cast<double>())
        /* init */,waveVectors((B*WaveVectorsAssembler<dim>::get(N,D).transpose()).transpose())
        {/*!\param[in] A_in the lattice matrix with lattice basis in column
          * \param[in] dens_in the size of the supercell along each lattice basis
          * \param[in] nums_in the number of subdivision along each supercell side, nums_in=dens_in for 1-st Brillouin zone, nums_in>dens_in for sub-cell waves
          */
            
            std::cout<<"waveVectors=\n"<<waveVectors<<std::endl;
            
            
            
        }
        
        PeriodicLatticeInterpolant(const MatrixDim& A_in,
                                   const VectorDimI& nums_in,
                                   const VectorDimI& dens_in,
                                   const Eigen::Matrix<double,Eigen::Dynamic,dim+1>& f,
                                   const Eigen::Matrix<double,Eigen::Dynamic,2*dim+1>& df) :
        /* init */ A(A_in)
        /* init */,B(2.0*M_PI*A.inverse().transpose())
        /* init */,N(nums_in)
        /* init */,D(dens_in.template cast<double>())
        /* init */,waveVectors((B*WaveVectorsAssembler<dim>::get(N,D).transpose()).transpose())
        /* init */,sinCosCoeffs(getSinCosCoeffs(f,df).transpose())
        {/*!\param[in] A_in the lattice matrix with lattice basis in column
          * \param[in] dens_in the size of the supercell along each lattice basis
          * \param[in] nums_in the number of subdivision along each supercell side, nums_in=dens_in for 1-st Brillouin zone, nums_in>dens_in for sub-cell waves
          */
            
            std::cout<<"waveVectors=\n"<<waveVectors<<std::endl;
            std::cout<<"sinCosCoeffs=\n"<<sinCosCoeffs<<std::endl;

        }

        void setConditions(const Eigen::Matrix<double,Eigen::Dynamic,dim+1>& f,
                           const Eigen::Matrix<double,Eigen::Dynamic,2*dim+1>& df)
        {
            sinCosCoeffs=getSinCosCoeffs(f,df).transpose();
            std::cout<<"sinCosCoeffs=\n"<<sinCosCoeffs<<std::endl;
        }
        
        Eigen::Map<Eigen::MatrixXd> getSinCosCoeffs(const Eigen::Matrix<double,Eigen::Dynamic,dim+1>& f,
                           const Eigen::Matrix<double,Eigen::Dynamic,2*dim+1>& df) const
        {/*!\param[in] f Matrix of function values in each row.
          * Each row has dim+1 entries with format
          * x1,x2,... f(x1,x2,...)
          * \param[in] df Matrix of function derivatives values in each row.
          * Each row has 2*dim+1 entries with format
          * x1,x2,...d1,d2,... dot(grad(f(x1,x2,...)),d)
          */
            Eigen::MatrixXd M1(Eigen::MatrixXd::Zero(f.rows(),2*waveVectors.rows())); // sine and cosine coeffs for each row
            Eigen::MatrixXd V1(Eigen::MatrixXd::Zero(f.rows(),1));
            for(int i=0;i<f.rows();i++)
            {
                for(int j=0;j<waveVectors.rows();j++)
                {
                    const double KdotX(waveVectors.row(j).dot(f.template block<1,dim>(i,0)));
                    M1(i,2*j)=sin(KdotX);
                    M1(i,2*j+1)=cos(KdotX);
                }
                V1(i)=f(i,2);
            }
            
            Eigen::MatrixXd M2(Eigen::MatrixXd::Zero(df.rows(),2*waveVectors.rows()));
            Eigen::MatrixXd V2(Eigen::MatrixXd::Zero(df.rows(),1));
            for(int i=0;i<df.rows();i++)
            {
                for(int j=0;j<waveVectors.rows();j++)
                {
                    const double dirNorm(df.template block<1,dim>(i,dim).norm());
                    const double KdotD(dirNorm>FLT_EPSILON? waveVectors.row(j).dot(df.template block<1,dim>(i,dim))/dirNorm : 0.0);
                    const double KdotX(waveVectors.row(j).dot(df.template block<1,dim>(i,0)));
                    M2(i,2*j)=KdotD*cos(KdotX);
                    M2(i,2*j+1)=-KdotD*sin(KdotX);
                }
                V2(i)=df(i,2*dim);
            }
            
            Eigen::MatrixXd M(M1.rows()+M2.rows(),2*waveVectors.rows()-1);
            M<<M1.block(0,1,M1.rows(),M1.cols()-1),M2.block(0,1,M2.rows(),M2.cols()-1);
            Eigen::MatrixXd V(V1.rows()+V2.rows(),1);
            V<<V1,V2;
            Eigen::MatrixXd x((M.transpose()*M).llt().solve(M.transpose()*V));
            Eigen::MatrixXd x1(x.rows()+1,1);
            x1<<0.0,x;
            std::cout<<"x1=\n"<<x1<<std::endl;
            return Eigen::Map<Eigen::MatrixXd>(x1.data(),2,x1.rows()/2);
        }
        
    };
    
} // end namespace
#endif

