/* This file is part of model.
 *
 * model is distributed without any warranty under the MIT License.
 */


#ifndef model_RLLL_cpp_
#define model_RLLL_cpp_

#include <RLLL.h>

namespace model
{


    /**********************************************************************/
    void RLLL::update(VectorType &H,
                MatrixType &M,
                const int &k)
    {
        const int n = B.cols();

        //            B1=B;
        //            B1(:,k-1)=B(:,k);
        //            B1(:,k)=B(:,k-1);
        B.col(k).swap(B.col(k - 1));

        MatrixType H1 = H;
        MatrixType M1 = M;
        H1(k - 1) = H(k) + std::pow(M(k, k - 1), 2) * H(k - 1);
        //        M1(k,k-1)=M(k,k-1)*H(k-1)/H1(k-1);
        M1(k, k - 1) *= H(k - 1) / H1(k - 1);
        H1(k) = H(k - 1) - std::pow(M1(k, k - 1), 2) * H1(k - 1);
        //        for i=k+1:n
        for (int i = k + 1; i < n; ++i)
        {
            M1(i, k - 1) = M(i, k - 1) * M1(k, k - 1) + M(i, k) * H(k) / H1(k - 1);
            M1(i, k) = M(i, k - 1) - M(i, k) * M(k, k - 1);
        }
        //            for j=1:k-2
        for (int j = 0; j <= k - 2; ++j)
        {
            M1(k - 1, j) = M(k, j);
            M1(k, j) = M(k - 1, j);
        }
        H = H1;
        M = M1;
    }

    /**********************************************************************/
    void RLLL::size_reduce(MatrixType &M,
                     const int &k,
                     const int &j)
    {

        const long int c = std::round(M(k, j));
        B.col(k) -= c * B.col(j);
        U.col(k) -= c * U.col(j);
        for (int l = 0; l <= j; ++l) // j is already 0-based, so use <=
        {
            M(k, l) -= c * M(j, l);
        }
    }

    /**********************************************************************/
    RLLL::RLLL(const MatrixType &B0,
         const double &delta) : /* init */ B(B0),
                                /* init */ U(Eigen::Matrix<long int, Eigen::Dynamic, Eigen::Dynamic>::Identity(B0.cols(), B0.cols()))
    {
        assert(delta >= 0.5 && delta <= 1.0 && "delta must be in [0.5 1]");

        const int n = B0.cols();
        assert(n <= B0.rows() && "B.rows() must be >= B.cols()");

        //            //std::cout<<"here 0"<<std::endl;
        VectorType H(VectorType::Zero(n));
        for (int j = 0; j < n; ++j)
        {
            H(j) = B.col(j).squaredNorm();
        }

        //std::cout<<H.transpose()<<std::endl;

        //std::cout<<"here 1"<<std::endl;
        MatrixType M(MatrixType::Identity(n, n));
        for (int j = 0; j < n; ++j)
        {
            for (int i = j + 1; i < n; ++i)
            {
                double temp = 0.0;
                //std::cout<<i<<" "<<j<<std::endl;

                for (int k = 0; k <= j - 1; ++k) // j is already 0-based, so use <=
                {
                    //std::cout<<i<<" "<<j<<" "<<k<<std::endl;

                    temp += M(j, k) * M(i, k) * H(k);
                }
                M(i, j) = (B.col(i).dot(B.col(j)) - temp) / H(j);
                H(i) -= std::pow(M(i, j), 2) * H(j);
            }
        }

        //std::cout<<std::setprecision(15)<<std::scientific<<M<<std::endl;
        //std::cout<<std::setprecision(15)<<std::scientific<<H<<std::endl;

        //            MatrixType M(MatrixType::Identity(n,n));
        int k = 1;
        while (k < n)
        {
            if (fabs(M(k, k - 1)) > 0.5)
            {
                size_reduce(M, k, k - 1);
            }

            if (H(k) < (delta - std::pow(M(k, k - 1), 2)) * H(k - 1))
            {
                update(H, M, k);
                U.col(k).swap(U.col(k - 1));
                k = std::max(1, k - 1);
            }
            else
            {
                //                    for j=k-2:-1:1
                for (int j = k - 2; j >= 0; --j)
                {
                    if (fabs(M(k, j)) > 0.5)
                    {
                        size_reduce(M, k, j);
                    }
                }
                k = k + 1;
            }
        }

        const double err = (B0.lu().solve(B) - U.cast<double>()).norm() / U.cast<double>().norm();
        if (err > FLT_EPSILON)
        {
            std::cout << "RLLL relative error= " << std::setprecision(15) << std::scientific << err << " > " << FLT_EPSILON << std::endl;
            assert(false && "Relative error too large. RLLL failed.");
        }

        const double absDetU(fabs(U.cast<double>().determinant()));
        if (fabs(absDetU - 1.0) > FLT_EPSILON)
        {
            std::cout << "|det(U)|= " << std::setprecision(15) << std::scientific << absDetU << std::endl;
            assert(false && "U is not unimodular. RLLL failed.");
        }
    }

    /**********************************************************************/
    const typename RLLL::MatrixType& RLLL::reducedBasis() const
    {
        return B;
    }

    /**********************************************************************/
    const Eigen::Matrix<long int, Eigen::Dynamic, Eigen::Dynamic>& RLLL::unimodularMatrix() const
    {
        return U;
    }

} // end namespace
#endif
