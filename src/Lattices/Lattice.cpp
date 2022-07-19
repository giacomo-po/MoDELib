/* This file is part of PIL, the Particle Interaction Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_Lattice_cpp_
#define model_Lattice_cpp_

#include <LatticeModule.h>

namespace model
{


    /**********************************************************************/
    template <int dim>
    typename Lattice<dim>::MatrixDimD Lattice<dim>::getLatticeBasis(const MatrixDimD &A, const MatrixDimD &Q)
    {

        // Check that Q is orthogonal
        const MatrixDimD QQT(Q * Q.transpose());
        if ((QQT - MatrixDimD::Identity()).norm() > 2.0 * DBL_EPSILON * dim * dim)
        {
            std::cout << "Q=\n"
                      << Q << std::endl;
            std::cout << "Q*Q^T=\n"
                      << QQT << std::endl;
            throw std::runtime_error("ROTATION MATRIX IS NOT ORTHOGONAL.");
//            assert(false && "ROTATION MATRIX IS NOT ORTHOGONAL.");
        }

        // Check sure that C2G is proper
        if(std::fabs(Q.determinant() - 1.0) > FLT_EPSILON)
        {
            throw std::runtime_error("ROTATION MATRIX IS NOT PROPER.");

        }
//        assert(std::fabs(Q.determinant() - 1.0) < FLT_EPSILON && "ROTATION MATRIX IS NOT PROPER.");

        // Check that A is full rank
        if(std::fabs(A.determinant()) < FLT_EPSILON)
        {
            throw std::runtime_error("A matrix is singular");

        }
//        assert(std::fabs(A.determinant()) > FLT_EPSILON && "A matrix is singular");

        return Q * A;
        //            const MatrixDimD QA(Q*A);
        //            return std::make_tuple(QA,QA.inverse().transpose(),Q);
    }

    //        /**********************************************************************/
    //        Lattice(const MatrixDimD& A,const MatrixDimD& invAT) :
    //        /* init */ latticeBases(std::make_pair(A,invAT)),
    //        /* init */ latticeBasis(latticeBases.first),
    //        /* init */ reciprocalBasis(latticeBases.second)
    //        {
    //
    //        }

    //        std::tuple<MatrixDimD,MatrixDimD,MatrixDimD> latticeBases;


    /**********************************************************************/
        /**********************************************************************/
    template <int dim>
    Lattice<dim>::Lattice(const MatrixDimD& A,const MatrixDimD& Q) :
//  /* init */ latticeBases(getLatticeBases(A,Q))
    /* init */ latticeBasis(getLatticeBasis(A,Q))
    /* init */,reciprocalBasis(latticeBasis.inverse().transpose())
    /* init */,C2G(Q)
    {

    }



    /**********************************************************************/
    template <int dim>
    Eigen::Matrix<long int, dim, 1> Lattice<dim>::rationalApproximation(VectorDimD nd)
    {
        Eigen::Array<long int, dim, 1> nums = Eigen::Matrix<long int, dim, 1>::Zero();

        if (nd.squaredNorm() > 0.0)
        {
            const Eigen::Array<double, dim, 1> nda(nd.array().abs()); // vector of close-to-integer numbers corresponding to lattice coordinates
            size_t maxID = 0;
            const double maxVal(nda.maxCoeff(&maxID));
            nd /= maxVal; // make each value of nd in [-1:1]

            nums = Eigen::Matrix<long int, dim, 1>::Ones();
            Eigen::Array<long int, dim, 1> dens = Eigen::Matrix<long int, dim, 1>::Ones();
            long int denProd = 1;

            for (int k = 0; k < dim; ++k)
            {
                BestRationalApproximation bra(nd(k), 10000);

                nums(k) = bra.num;
                dens(k) = bra.den;
                denProd *= bra.den;
            }

            for (int k = 0; k < dim; ++k)
            {
                nums(k) *= (denProd / dens(k));
            }
        }

        return nums.matrix();
    }

    /**********************************************************************/
    template <int dim>
    typename Lattice<dim>::LatticeVectorType Lattice<dim>::snapToLattice(const VectorDimD &d) const
    {
        VectorDimD nd(reciprocalBasis.transpose() * d);
        return LatticeVectorType(nd.array().round().matrix().template cast<long int>(), *this);
    }

    /**********************************************************************/
    template <int dim>
    typename Lattice<dim>::LatticeDirectionType Lattice<dim>::latticeDirection(const VectorDimD &d) const
    {

        const VectorDimD nd(reciprocalBasis.transpose() * d);
        const LatticeVectorType temp(rationalApproximation(nd), *this);

        const double crossNorm(temp.cartesian().normalized().cross(d.normalized()).norm());
        if (crossNorm > FLT_EPSILON)
        {
            std::cout << "input direction=" << d.normalized().transpose() << std::endl;
            std::cout << "lattice direction=" << temp.cartesian().normalized().transpose() << std::endl;
            std::cout << "cross product norm=" << std::setprecision(15) << std::scientific << crossNorm << std::endl;
            throw std::runtime_error("LATTICE DIRECTION NOT FOUND");
        }

        return LatticeDirectionType(temp);
    }

    /**********************************************************************/
    template <int dim>
    typename Lattice<dim>::ReciprocalLatticeDirectionType Lattice<dim>::reciprocalLatticeDirection(const VectorDimD &d) const
    {

        const VectorDimD nd(latticeBasis.transpose() * d);
        const ReciprocalLatticeVectorType temp(rationalApproximation(nd), *this);

        const double crossNorm(temp.cartesian().normalized().cross(d.normalized()).norm());
        if (crossNorm > FLT_EPSILON)
        {
            std::cout << "input direction=" << std::setprecision(15) << std::scientific << d.normalized().transpose() << std::endl;
            std::cout << "reciprocal lattice direction=" << std::setprecision(15) << std::scientific << temp.cartesian().normalized().transpose() << std::endl;
            std::cout << "cross product norm=" << std::setprecision(15) << std::scientific << crossNorm << std::endl;
//            assert(0 && "RECIPROCAL LATTICE DIRECTION NOT FOUND");
            throw std::runtime_error("RECIPROCAL LATTICE DIRECTION NOT FOUND");

        }

        return ReciprocalLatticeDirectionType(temp);
    }

    /**********************************************************************/
    template <int dim>
    typename Lattice<dim>::RationalLatticeDirectionType Lattice<dim>::rationalLatticeDirection(const VectorDimD &d,
                                                          const typename BestRationalApproximation::LongIntType &maxDen ) const
    {

        const LatticeDirectionType ld(latticeDirection(d));
        const BestRationalApproximation bra(d.norm() / ld.cartesian().norm(), maxDen);
        Rational rat(bra.num, bra.den);
        RationalLatticeDirectionType rld(rat, ld);
        if ((rld.cartesian() - d).squaredNorm() > FLT_EPSILON)
        {
            std::cout << "input vector=" << d.transpose() << std::endl;
            std::cout << "lattice direction=" << ld.cartesian().transpose() << std::endl;
            std::cout << "rational=" << rat << std::endl;
            std::cout << "d.norm()/ld.cartesian().norm()=" << d.norm() / ld.norm() << std::endl;
//            assert(0 && "RationalLatticeDirectionType NOT FOUND");
            throw std::runtime_error("RationalLatticeDirectionType NOT FOUND");
        }
        return rld;
    }

    /**********************************************************************/
    template <int dim>
    typename Lattice<dim>::LatticeVectorType Lattice<dim>::latticeVector(const VectorDimD &p) const
    {
        return LatticeVectorType(p, *this);
    }

    /**********************************************************************/
    template <int dim>
    typename Lattice<dim>::ReciprocalLatticeVectorType Lattice<dim>::reciprocalLatticeVector(const VectorDimD &p) const
    {
        return ReciprocalLatticeVectorType(p, *this);
    }
    template class Lattice<3>;
}
#endif
