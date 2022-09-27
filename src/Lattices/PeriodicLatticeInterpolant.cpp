/* This file is part of PIL, the Particle Interaction Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PeriodicLatticeInterpolant_cpp_
#define model_PeriodicLatticeInterpolant_cpp_

#include <numbers>
#include <PeriodicLatticeInterpolant.h>

namespace model
{


   

    /**********************************************************************/
    template <int dim>
    std::vector<typename PeriodicLatticeInterpolant<dim>::MatrixDim> PeriodicLatticeInterpolant<dim>::getSymMatrices(const int &rotSymm,
                                                 const std::vector<Eigen::Matrix<double, dim, 1>> &mirSymm)
    {
        std::vector<MatrixDim> temp;
        for (int k = 0; k < rotSymm; ++k)
        { // apply rotations
            const double theta(k * 2.0 * std::numbers::pi / rotSymm);
            const double c(cos(theta));
            const double s(sin(theta));
            temp.push_back((MatrixDim() << c, -s, s, c).finished());
        }

        for (const auto &N : mirSymm)
        { // apply mirSymm to each existing entry in temp
            const double nNorm(N.norm());
            if (nNorm > FLT_EPSILON)
            {
                const VectorDim n(N / nNorm);
                const auto temp1 = temp;
                for (const auto &m : temp1)
                {
                    temp.push_back(m * (MatrixDim::Identity() - 2.0 * n * n.transpose())); // mirror symm
                }
            }
        }
        return temp;
    }

    /**********************************************************************/
    template <int dim>
    Eigen::MatrixXd PeriodicLatticeInterpolant<dim>::getSinCosCoeffs(const Eigen::Matrix<double, Eigen::Dynamic, dim + 1> &f) const
    {                                                                               /*!\param[in] f Matrix of function values in each row.
          * Each row has dim+1 entries with format
          * x1,x2,... f(x1,x2,...)
          * \param[in] df Matrix of function derivatives values in each row.
          * Each row has 2*dim+1 entries with format
          * x1,x2,...d1,d2,... dot(grad(f(x1,x2,...)),d)
          */
        Eigen::MatrixXd M(Eigen::MatrixXd::Zero(f.rows(), 2 * waveVectors.rows())); // sine and cosine coeffs for each row
        Eigen::MatrixXd V(Eigen::MatrixXd::Zero(f.rows(), 1));
        for (int i = 0; i < f.rows(); i++)
        {
            for (int j = 0; j < waveVectors.rows(); j++)
            {
                const auto sc(sinKXcosKX(j, f.template block<1, dim>(i, 0)));
                M(i, 2 * j) = sc.first;
                M(i, 2 * j + 1) = sc.second;
            }
            V(i) = f(i, 2);
        }

        // Some of the S and C coefficients may be identically zero. This means that some columns of M may be vanishing, we need to remove those colums
        std::vector<int> nonZcols;
        for (int j = 0; j < M.cols(); ++j)
        {
            if (M.col(j).squaredNorm() > FLT_EPSILON)
            {
                nonZcols.push_back(j);
            }
        }

        Eigen::MatrixXd M1(M.rows(), nonZcols.size());
        for (size_t j = 0; j < nonZcols.size(); ++j)
        {
            M1.col(j) = M.col(nonZcols[j]);
        }

        const Eigen::JacobiSVD<Eigen::MatrixXd> svd(M1);

        // Solve
        if (M1.rows() == M1.cols() && svd.rank() == M1.rows())
        {
            const double cond(svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size() - 1));
            if (cond > 1.0e3)
            {
                std::cout << "M1=\n"
                          << M1 << std::endl;
                std::cout << "rank=" << svd.rank() << std::endl;
                std::cout << "cond=" << cond << std::endl;
                throw std::runtime_error("Bad condition number");
            }
            const Eigen::VectorXd x(M1.lu().solve(V));

            Eigen::VectorXd x1(Eigen::VectorXd::Zero(M.cols()));
            for (size_t j = 0; j < nonZcols.size(); ++j)
            {
                x1(nonZcols[j]) = x(j);
            }
            return Eigen::Map<Eigen::MatrixXd>(x1.data(), 2, x1.rows() / 2); // reshape
        }
        else
        {
            std::cout << "M=\n"
                      << M << std::endl;
            std::cout << "M1=\n"
                      << M1 << std::endl;
            std::cout << "M1 has size " << M1.rows() << "x" << M1.cols() << ",and rank=" << svd.rank() << std::endl;
            throw std::runtime_error("M1 MUST BE SQUARE AND FULL RANK");
            return Eigen::MatrixXd::Zero(1, 1);
        }
    }

    /**********************************************************************/
    template <int dim>
    std::pair<double, double> PeriodicLatticeInterpolant<dim>::sinKXcosKX(const int &r, const VectorDim &x) const
    {
        std::pair<double, double> sc(0.0, 0.0);
        for (const auto &R : symMatrices)
        {
            const double KdotX(waveVectors.row(r).dot(R * x));
            sc.first += sin(KdotX);
            sc.second += cos(KdotX);
        }
        return sc;
    }

    /**********************************************************************/
    template <int dim>
    PeriodicLatticeInterpolant<dim>::PeriodicLatticeInterpolant(const MatrixDim &A_in,
                               const Eigen::Matrix<double, Eigen::Dynamic, dim> &waveVectors_in,
                               const Eigen::Matrix<double, Eigen::Dynamic, dim + 1> &f,
                               const int &rotSymm,
                               const std::vector<Eigen::Matrix<double, dim, 1>> &N) : // rotSymm only works for 2d, in 3d we need to prescride axis and n-fold symm
                                                                                      /* init */ A(A_in)
                                                                                      /* init */,
                                                                                      B(2.0 * std::numbers::pi * A.inverse().transpose())
                                                                                      /* init */,
                                                                                      waveVectors((B * waveVectors_in.transpose()).transpose())
                                                                                      /* init */,
                                                                                      symMatrices(getSymMatrices(rotSymm, N))
                                                                                      /* init */,
                                                                                      sinCosCoeffs(getSinCosCoeffs(f).transpose())
    { /*!\param[in] A_in the lattice matrix with lattice basis in column
          * \param[in] 
          */
    }

    /**********************************************************************/
    template <int dim>
    double PeriodicLatticeInterpolant<dim>::operator()(const VectorDim &x) const
    {
        double temp(0.0);
        for (int r = 0; r < waveVectors.rows(); ++r)
        {
            const auto sc(sinKXcosKX(r, x));
            temp += sinCosCoeffs(r, 0) * sc.first + sinCosCoeffs(r, 1) * sc.second;
        }
        return temp;
    }
    
template struct PeriodicLatticeInterpolant<2>;
} // end namespace
#endif
