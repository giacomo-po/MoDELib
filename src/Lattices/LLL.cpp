/* This file is part of PIL, the Particle Interaction Library.
 *
 * Copyright (C) 2016 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LLL_cpp_
#define model_LLL_cpp_

#include <LLL.h>
// http://www.arageli.org/download
// https://www.mathworks.com/matlabcentral/fileexchange/49457-lattice-reduction-mimo?focused=3859922&tab=function

namespace model
{
    /*! An implementation of the  Lenstra–Lenstra–Lovász (LLL) lattice basis
     *  reduction algorithm for integers.
     */

    /**********************************************************************/
    void LLL::lll_gram_schmidt_int(const int &k)
    {

        for (int j = 0; j <= k; j++)
        {
            int u = B.col(k).dot(B.col(j));
            for (int i = 0; i < j; i++)
            {
                u = (d(i + 1) * u - Lambda(k, i) * Lambda(j, i)) / d(i);
            }
            if (j < k)
            {
                Lambda(k, j) = u;
            }
            else
            {
                d(k + 1) = u;
            }
        }
    }

    /**********************************************************************/
    void LLL::lll_size_reduction_int(const int &k,
                                     const int &l)
    {

        int q = (2 * Lambda(k, l) + d[l + 1]) / (2 * d(l + 1)); //quotient
        B.col(k) -= B.col(l) * q;
        H.col(k) -= H.col(l) * q;
        Lambda(k, l) -= q * d(l + 1);

        for (int i = 0; i < l; i++)
        {
            Lambda(k, i) -= q * Lambda(l, i);
        }
    }

    /**********************************************************************/
    void LLL::lll_interchange_int(const int &k,
                                  const int &k_max)
    {
        //            typedef typename Lambda_type::value_type T;
        std::cout << "1" << std::endl;

        Eigen::VectorXi tempCol = B.col(k);
        B.col(k) = B.col(k - 1);
        B.col(k - 1) = tempCol;

        std::cout << "2" << std::endl;

        tempCol = H.col(k);
        H.col(k) = H.col(k - 1);
        H.col(k - 1) = tempCol;

        //            B.swap_cols(k, k - 1);
        //            H.swap_cols(k, k - 1);

        std::cout << "3" << std::endl;

        for (int j = 0; j < k - 1; j++)
        {
            std::swap(Lambda(k, j), Lambda(k - 1, j));
        }

        std::cout << "4" << std::endl;

        int lambda = Lambda(k, k - 1);
        int b = (d(k - 1) * d(k + 1) + pow(lambda, 2)) / d(k);

        std::cout << "5" << std::endl;
        std::cout << "d(k+1)=" << d(k + 1) << std::endl;

        for (int i = k + 1; i <= k_max; i++)
        {
            int t = Lambda(i, k);
            Lambda(i, k) = (d(k + 1) * Lambda(i, k - 1) - lambda * t) / d(k);
            Lambda(i, k - 1) = (b * t + lambda * Lambda(i, k)) / d(k + 1);
        }
        std::cout << "6" << std::endl;

        d(k) = b;
        std::cout << "7" << std::endl;
    }

    template <int m, int n>
    LLL::LLL(const Eigen::Matrix<int, m, n> &B_in) : /* init */ B(B_in),
                                                     /* init */ d(Eigen::Matrix<int, 1, n + 1>::Ones()),
                                                     /* init */ H(Eigen::Matrix<int, n, n>::Identity()),
                                                     /* init */ Lambda(Eigen::Matrix<int, n, n>::Identity())
    {

        d(0) = 1;
        d(1) = B.col(0).squaredNorm();

        int k_max = 0;
        for (int k = 1; k < n;)
        {
            std::cout << "k=" << k << std::endl;
            if (k > k_max)
            {
                std::cout << "0" << std::endl;
                k_max = k;

                lll_gram_schmidt_int(k);

                //                                        if (is_null(d[k + 1]))
                if (d(k + 1) == 0)
                {
                    assert(0 && "B(i) did not form the basis");
                    //                                            return false; // B(i) did not form the basis
                }
            }

            if (2 * std::abs(Lambda(k, k - 1)) > d(k))
            {
                std::cout << "A" << std::endl;
                lll_size_reduction_int(k, k - 1);
            }

            if (4 * d(k + 1) * d(k - 1) < 3 * pow(d(k), 2) - 4 * pow(Lambda(k, k - 1), 2))
            {
                std::cout << "B" << std::endl;
                lll_interchange_int(k, k_max);
                k = std::max(1, k - 1);
                std::cout << "B1" << std::endl;
            }
            else
            {
                std::cout << "C" << std::endl;
                for (int l = k - 1; l > 0; l--)
                    if (2 * std::abs(Lambda(k, l - 1)) > d(l))
                        lll_size_reduction_int(k, l - 1);
                k++;
            }
            std::cout << "B=" << B << std::endl;
            std::cout << "d=" << d.transpose() << std::endl;
        }

        std::cout << B << std::endl;
    }

} // end namespace
#endif
