/* This file is part of model.
 *
 * model is distributed without any warranty under the MIT License.
 */


#ifndef model_RationalMatrix_cpp_
#define model_RationalMatrix_cpp_

#include <iomanip>
#include <RationalMatrix.h>

namespace model
{

    /**********************************************************************/
    template <int dim>
    std::pair<typename RationalMatrix<dim>::MatrixDimI, typename RationalMatrix<dim>::IntScalarType> RationalMatrix<dim>::compute(const MatrixDimD &R)
    {

        // Find the BestRationalApproximation of each entry
        MatrixDimI nums(MatrixDimI::Zero());
        MatrixDimI dens(MatrixDimI::Ones());

        IntScalarType sigma = 1;
        for (int i = 0; i < dim; ++i)
        {
            for (int j = 0; j < dim; ++j)
            {
                BestRationalApproximation bra(R(i, j), maxDen);
                nums(i, j) = bra.num;
                dens(i, j) = bra.den;
                sigma = IntegerMath<IntScalarType>::lcm(sigma, bra.den);
            }
        }

        MatrixDimI im(MatrixDimI::Zero());
        for (int i = 0; i < dim; ++i)
        {
            for (int j = 0; j < dim; ++j)
            {
                im(i, j) = nums(i, j) * sigma / dens(i, j);
            }
        }

        const double error = (im.template cast<double>() / sigma - R).norm() / (dim * dim);
        if (error > 100.0 * DBL_EPSILON)
        {
            std::cout << "error=" << error << std::endl;
            std::cout << "maxDen=" << maxDen << std::endl;
            std::cout << "im=\n"
                      << std::setprecision(15) << std::scientific << im.template cast<double>() / sigma << std::endl;
            std::cout << "= 1/" << sigma << "*\n"
                      << std::setprecision(15) << std::scientific << im << std::endl;
            std::cout << "R=\n"
                      << std::setprecision(15) << std::scientific << R << std::endl;
            throw std::runtime_error("Rational Matrix failed, check maxDen");
        }

        return std::make_pair(im, sigma);
    }

    /**********************************************************************/
    template <int dim>
    RationalMatrix<dim>::RationalMatrix(const MatrixDimD &R) :
    /* init */ returnPair(compute(R)),
    /* init */ integerMatrix(returnPair.first),
    /* init */ mu(returnPair.second)
    {
    }

    template <int dim>
    RationalMatrix<dim>::RationalMatrix(const MatrixDimI& Rn,const IntScalarType& Rd) :
    /* init */ returnPair(std::make_pair(Rn, Rd)),
    /* init */ integerMatrix(returnPair.first),
    /* init */ mu(returnPair.second)
    {
    }

    template <int dim>
    typename RationalMatrix<dim>::MatrixDimD RationalMatrix<dim>::asMatrix() const
    {
        return integerMatrix.template cast<double>()/mu;
    }



    template class RationalMatrix<1>;
    template class RationalMatrix<2>;
    template class RationalMatrix<3>;

} // end namespace
#endif
