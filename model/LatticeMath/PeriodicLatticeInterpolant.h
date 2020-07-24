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
    struct PeriodicLatticeInterpolant
    {
        
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef Eigen::Matrix<size_t,dim,1> VectorDimI;
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;

        
        const MatrixDim A;   // column matrix of lattice vectors
        const MatrixDim B;   // column matrix of reciprocal lattice vectors
//        const VectorDimI N; // numerators of k-point indices
//        const VectorDim D;   // denumerators of k-point indices
        
        Eigen::MatrixXd getSinCosCoeffs(const Eigen::Matrix<double,Eigen::Dynamic,dim+1>& f,
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
            
            
            assert(M.rows()>=M.cols() && "MUST HAVE AT LEAST AS MANY EQUATIONS AS UNKNOWNS");
            
            // Solve Yx=z;
            Eigen::MatrixXd x;
            bool useLeastSquares(false);
            if(useLeastSquares)
            {
                const Eigen::MatrixXd Y(M.transpose()*M);
                const Eigen::JacobiSVD<Eigen::MatrixXd> svd(Y);
                const double cond(svd.singularValues()(0)/ svd.singularValues()(svd.singularValues().size()-1));
                if(cond>1.0e3)
                {
                    std::cout<<"Y=\n"<<Y<<std::endl;
                    std::cout<<"cond="<<cond<<std::endl;
                                                    assert(cond<1.0e3 && "Bad condition number");
                }
                x=Y.llt().solve(M.transpose()*V);
            }
            else
            {// Y=M
                const Eigen::JacobiSVD<Eigen::MatrixXd> svd(M);
                const double cond(svd.singularValues()(0)/ svd.singularValues()(svd.singularValues().size()-1));
                if(cond>1.0e3)
                {
                    std::cout<<"M=\n"<<M<<std::endl;
                    std::cout<<"cond="<<cond<<std::endl;
                                assert(cond<1.0e3 && "Bad condition number");
                }
                x=M.lu().solve(V);
            }
//            Eigen::MatrixXd x(Y.llt().solve(M.transpose()*V));
            Eigen::MatrixXd x1(x.rows()+1,1);
            x1<<0.0,x;
            return Eigen::Map<Eigen::MatrixXd>(x1.data(),2,x1.rows()/2);
        }

//        Eigen::MatrixXd getWaveVecFromIntCoord(const Eigen::Matrix<double,Eigen::Dynamic,dim>& waveVecInt) const
//        {
//            std::cout<<"\n Inside pli constructor"<<std::endl;
//            MatrixDim Bmod;
//            Bmod.col(0)=B.col(0)/D(0);
//            Bmod.col(1)=B.col(1)/D(1);
//            //std::cout<<"\n B:"<<Bmod;
//            //std::cout<<"\n waveVecInt_in:"<<waveVecInt_in;
//            //waveVectors.resize(waveVecInt_in.rows(),dim);
//            return waveVecInt * Bmod.transpose();
//        }
        
        
        const Eigen::Matrix<double,Eigen::Dynamic,dim> waveVectors;
        const Eigen::Matrix<double,Eigen::Dynamic,2> sinCosCoeffs;

        
        
        PeriodicLatticeInterpolant(const MatrixDim& A_in,
                                   const VectorDimI& N,
                                   const VectorDimI& D,
                                   const Eigen::Matrix<double,Eigen::Dynamic,dim+1>& f,
                                   const Eigen::Matrix<double,Eigen::Dynamic,2*dim+1>& df) :
        /* init */ A(A_in)
        /* init */,B(2.0*M_PI*A.inverse().transpose())
        /* init */,waveVectors((B*WaveVectorsAssembler<dim>::get(N,D.template cast<double>()).transpose()).transpose())
        /* init */,sinCosCoeffs(getSinCosCoeffs(f,df).transpose()) 
        {/*!\param[in] A_in the lattice matrix with lattice basis in column
          * \param[in] dens_in the size of the supercell along each lattice basis
          * \param[in] nums_in the number of subdivision along each supercell side, nums_in=dens_in for 1-st Brillouin zone, nums_in>dens_in for sub-cell waves
          */
            
        }

        PeriodicLatticeInterpolant(const MatrixDim& A_in,
                                   const Eigen::Matrix<double,Eigen::Dynamic,dim>& waveVectors_in,
                                   const Eigen::Matrix<double,Eigen::Dynamic,dim+1>& f,
                                   const Eigen::Matrix<double,Eigen::Dynamic,2*dim+1>& df) :
        /* init */ A(A_in)
        /* init */,B(2.0*M_PI*A.inverse().transpose())
        /* init */,waveVectors((B*waveVectors_in.transpose()).transpose())
        /* init */,sinCosCoeffs(getSinCosCoeffs(f,df).transpose())
        { /*!\param[in] A_in the lattice matrix with lattice basis in column
          * \param[in] 
          */

        }
        
        
        double operator()(const VectorDim& x)
        {
            double temp(0);
            for(int r=0;r<waveVectors.rows();++r)
            {
                const double KdotX(waveVectors.row(r).dot(x));
                temp+=sinCosCoeffs(r,0)*sin(KdotX)+sinCosCoeffs(r,1)*cos(KdotX);
            }
            return temp;
        }

    };
    
} // end namespace
#endif

