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
                                                     const std::vector<Eigen::Matrix<double,dim,1>>& mirSymm)
        {
            std::vector<MatrixDim> temp;
            for(int k=0;k<rotSymm;++k)
            {// apply rotations
                const double theta(k*2.0*M_PI/rotSymm);
                const double c(cos(theta));
                const double s(sin(theta));
                temp.push_back((MatrixDim()<<c,-s,s,c).finished());
            }
            
            for(const auto& N : mirSymm)
            {// apply mirSymm to each existing entry in temp
                const double nNorm(N.norm());
                if(nNorm>FLT_EPSILON)
                {
                    const VectorDim n(N/nNorm);
                    const auto temp1=temp;
                    for(const auto& m : temp1)
                    {
                        temp.push_back(m*(MatrixDim::Identity()-2.0*n*n.transpose())); // mirror symm
                    }
                }
            }
            return temp;
        }
        
        /**********************************************************************/
        Eigen::MatrixXd getSinCosCoeffs(const Eigen::Matrix<double,Eigen::Dynamic,dim+1>& f) const
        {/*!\param[in] f Matrix of function values in each row.
          * Each row has dim+1 entries with format
          * x1,x2,... f(x1,x2,...)
          * \param[in] df Matrix of function derivatives values in each row.
          * Each row has 2*dim+1 entries with format
          * x1,x2,...d1,d2,... dot(grad(f(x1,x2,...)),d)
          */
            Eigen::MatrixXd M(Eigen::MatrixXd::Zero(f.rows(),2*waveVectors.rows())); // sine and cosine coeffs for each row
            Eigen::MatrixXd V(Eigen::MatrixXd::Zero(f.rows(),1));
            for(int i=0;i<f.rows();i++)
            {
                for(int j=0;j<waveVectors.rows();j++)
                {
                    const auto sc(sinKXcosKX(j,f.template block<1,dim>(i,0)));
                    M(i,2*j)  =sc.first;
                    M(i,2*j+1)=sc.second;
                }
                V(i)=f(i,2);
            }
            
            // Some of the S and C coefficients may be identically zero. This means that some columns of M may be vanishing, we need to remove those colums
            std::vector<int> nonZcols;
            for(int j=0;j<M.cols();++j)
            {
                if(M.col(j).squaredNorm()>FLT_EPSILON)
                {
                    nonZcols.push_back(j);
                }
            }
            
            Eigen::MatrixXd M1(M.rows(),nonZcols.size());
            for(size_t j=0;j<nonZcols.size();++j)
            {
                M1.col(j)=M.col(nonZcols[j]);
            }
            
            const Eigen::JacobiSVD<Eigen::MatrixXd> svd(M1);
            
            // Solve
            if(M1.rows()==M1.cols() && svd.rank()==M1.rows())
            {
                const double cond(svd.singularValues()(0)/ svd.singularValues()(svd.singularValues().size()-1));
                if(cond>1.0e3)
                {
                    std::cout<<"M1=\n"<<M1<<std::endl;
                    std::cout<<"rank="<<svd.rank()<<std::endl;
                    std::cout<<"cond="<<cond<<std::endl;
                    assert(cond<1.0e3 && "Bad condition number");
                }
                const Eigen::VectorXd x(M1.lu().solve(V));

                Eigen::VectorXd x1(Eigen::VectorXd::Zero(M.cols()));
                for(size_t j=0;j<nonZcols.size();++j)
                {
                    x1(nonZcols[j])=x(j);
                }
                return Eigen::Map<Eigen::MatrixXd>(x1.data(),2,x1.rows()/2); // reshape
            }
            else
            {
                std::cout<<"M=\n"<<M<<std::endl;
                std::cout<<"M1=\n"<<M1<<std::endl;
                std::cout<<"M1 has size "<<M1.rows()<<"x"<<M1.cols()<<",and rank="<<svd.rank()<<std::endl;
                assert(false && "M1 MUST BE SQUARE AND FULL RANK");
                return Eigen::MatrixXd::Zero(1,1);
            }
        }
        
        /**********************************************************************/
        std::pair<double,double> sinKXcosKX(const int& r,const VectorDim& x) const
        {
            std::pair<double,double> sc(0.0,0.0);
            for(const auto& R : symMatrices)
            {
                const double KdotX(waveVectors.row(r).dot(R*x));
                sc.first +=sin(KdotX);
                sc.second+=cos(KdotX);
            }
            return sc;
        }

        /**********************************************************************/
        PeriodicLatticeInterpolant(const MatrixDim& A_in,
                                   const Eigen::Matrix<double,Eigen::Dynamic,dim>& waveVectors_in,
                                   const Eigen::Matrix<double,Eigen::Dynamic,dim+1>& f,
                                   const int& rotSymm,
                                   const std::vector<Eigen::Matrix<double,dim,1>>& N) : // rotSymm only works for 2d, in 3d we need to prescride axis and n-fold symm
        /* init */ A(A_in)
        /* init */,B(2.0*M_PI*A.inverse().transpose())
        /* init */,waveVectors((B*waveVectors_in.transpose()).transpose())
        /* init */,symMatrices(getSymMatrices(rotSymm,N))
        /* init */,sinCosCoeffs(getSinCosCoeffs(f).transpose())
        { /*!\param[in] A_in the lattice matrix with lattice basis in column
          * \param[in] 
          */

        }
        
        /**********************************************************************/
        double operator()(const VectorDim& x) const
        {
            double temp(0.0);
            for(int r=0;r<waveVectors.rows();++r)
            {
                const auto sc(sinKXcosKX(r,x));
                temp+=sinCosCoeffs(r,0)*sc.first+sinCosCoeffs(r,1)*sc.second;
            }
            return temp;
        }

    };
    
} // end namespace
#endif


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





//        PeriodicLatticeInterpolant(const MatrixDim& A_in,
//                                   const VectorDimI& N,
//                                   const VectorDimI& D,
//                                   const Eigen::Matrix<double,Eigen::Dynamic,dim+1>& f,
//                                   const Eigen::Matrix<double,Eigen::Dynamic,2*dim+1>& df) :
//        /* init */ A(A_in)
//        /* init */,B(2.0*M_PI*A.inverse().transpose())
//        /* init */,waveVectors((B*WaveVectorsAssembler<dim>::get(N,D.template cast<double>()).transpose()).transpose())
//        /* init */,sinCosCoeffs(getSinCosCoeffs(f,df).transpose())
//        {/*!\param[in] A_in the lattice matrix with lattice basis in column
//          * \param[in] dens_in the size of the supercell along each lattice basis
//          * \param[in] nums_in the number of subdivision along each supercell side, nums_in=dens_in for 1-st Brillouin zone, nums_in>dens_in for sub-cell waves
//          */
//
//        }


//
//    template <int dim>
//    struct WaveVectorsAssembler
//    {
//        typedef Eigen::Matrix<size_t,dim,1> VectorDimI;
//        typedef Eigen::Matrix<double,dim,1> VectorDim;
//
//
//        static Eigen::Matrix<double,Eigen::Dynamic,dim> get(const VectorDimI& N,
//                                                            const VectorDim& D)
//        {
//            const size_t rows(N.prod());
//            const auto lowerBlock(WaveVectorsAssembler<dim-1>::get(N.template segment<dim-1>(1),D.template segment<dim-1>(1)));
//            Eigen::Matrix<double,Eigen::Dynamic,dim> temp(rows,dim);
//            for(size_t k=0;k<N(0);++k)
//            {
//                temp.block(k*lowerBlock.rows(),0,lowerBlock.rows(),dim)<<Eigen::MatrixXd::Constant(lowerBlock.rows(), 1, k/D(0)),lowerBlock;
//            }
//            return temp;
//        }
//
//    };
//
//
//    template<>
//    struct WaveVectorsAssembler<1>
//    {
//        static constexpr int dim=1;
//        typedef Eigen::Matrix<size_t,dim,1> VectorDimI;
//        typedef Eigen::Matrix<double,dim,1> VectorDim;
//
//
//        static Eigen::Matrix<double,Eigen::Dynamic,dim> get(const VectorDimI& N,
//                                                            const VectorDim& D)
//        {
//            const size_t rows(N.prod());
//            Eigen::Matrix<double,Eigen::Dynamic,dim> temp(rows,dim);
//            for(size_t k=0;k<N(0);++k)
//            {
//                temp(k,0)=k/D(0);
//            }
//            return temp;
//        }
//
//    };
