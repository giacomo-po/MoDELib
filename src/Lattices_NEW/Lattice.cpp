/* This file is part of model.
 *
 * model is distributed without any warranty under the MIT License.
 */


#ifndef model_Lattice_cpp_
#define model_Lattice_cpp_

#include <Eigen/Eigenvalues>

#include <LatticeModule.h>
#include <GramMatrix.h>
#include <iomanip>

namespace model
{


    /**********************************************************************/
    template <int dim>
    typename Lattice<dim>::MatrixDimD Lattice<dim>::getLatticeBasis(const MatrixDimD &A, const MatrixDimD &Q)
    {

        // Check that Q is orthogonal
        const MatrixDimD QQT(Q * Q.transpose());
        const double QQTerror((QQT - MatrixDimD::Identity()).norm());
        if (QQTerror > 10.0 * DBL_EPSILON * dim * dim)
        {
            throw std::runtime_error("The rotation matrix Q is not orthogonal: norm(Q*Q^T-I)="+std::to_string(QQTerror)+"\n");
        }

        const double Qdet(Q.determinant());
        if (std::fabs(Q.determinant() - 1.0)>FLT_EPSILON)
        {
            throw std::runtime_error("The rotation matrix is not proper: det(Q)="+std::to_string(Qdet)+"\n");
        }
        
        return Q * A;
    }


    /**********************************************************************/
        /**********************************************************************/
    template <int dim>
    Lattice<dim>::Lattice(const MatrixDimD& A,const MatrixDimD& Q) :
    /* init */ latticeBasis(getLatticeBasis(A,Q))
    /* init */,reciprocalBasis(latticeBasis.inverse().transpose())
    /* init */,C2G(Q)
    {

    }

    /**********************************************************************/
    template <int dim>
    LatticeDirection<dim> Lattice<dim>::latticeDirection(const VectorDimD &d) const
    {
        const VectorDimD nd(reciprocalBasis.transpose()*d);
        const LatticeVector<dim> temp(LatticeCore<dim>::rationalApproximation(nd),*this);
        const GramMatrix<double,2> G(std::array<VectorDimD,2>{temp.cartesian().normalized(),d.normalized()});
        const double crossNorm(sqrt(G.determinant()));
        if(crossNorm>FLT_EPSILON)
        {
            std::cout<<"input direction="<<d.normalized().transpose()<<std::endl;
            std::cout<<"lattice direction="<<temp.cartesian().normalized().transpose()<<std::endl;
            std::cout<<"cross product norm="<<std::setprecision(15)<<std::scientific<<crossNorm<<std::endl;
            throw std::runtime_error("LATTICE DIRECTION NOT FOUND\n");
        }
        return LatticeDirection<dim>(temp);
    }

    /**********************************************************************/
    template <int dim>
    ReciprocalLatticeDirection<dim> Lattice<dim>::reciprocalLatticeDirection(const VectorDimD &d) const
    {
        const VectorDimD nd(latticeBasis.transpose()*d);
        const ReciprocalLatticeVector<dim> temp(LatticeCore<dim>::rationalApproximation(nd),*this);
        const GramMatrix<double,2> G(std::array<VectorDimD,2>{temp.cartesian().normalized(),d.normalized()});
        const double crossNorm(sqrt(G.determinant()));
        if(crossNorm>FLT_EPSILON)
        {
            std::cout<<"input direction="<<std::setprecision(15)<<std::scientific<<d.normalized().transpose()<<std::endl;
            std::cout<<"reciprocal lattice direction="<<std::setprecision(15)<<std::scientific<<temp.cartesian().normalized().transpose()<<std::endl;
            std::cout<<"cross product norm="<<std::setprecision(15)<<std::scientific<<crossNorm<<std::endl;
            throw std::runtime_error("RECIPROCAL LATTICE DIRECTION NOT FOUND\n");
        }
        return ReciprocalLatticeDirection<dim>(temp);
    }

    /**********************************************************************/
    template <int dim>
    RationalLatticeDirection<dim> Lattice<dim>::rationalLatticeDirection(const VectorDimD &d,
                                                          const typename BestRationalApproximation::LongIntType &maxDen ) const
    {
        const LatticeDirection<dim> ld(latticeDirection(d));
        const BestRationalApproximation bra(d.norm() / ld.cartesian().norm(), maxDen);
        const Rational rat(bra.num, bra.den);
        const RationalLatticeDirection<dim> rld(rat, ld);
        if ((rld.cartesian() - d).squaredNorm() > FLT_EPSILON)
        {
            std::cout << "input vector=" << d.transpose() << std::endl;
            std::cout << "lattice direction=" << ld.cartesian().transpose() << std::endl;
            std::cout << "rational=" << rat << std::endl;
            std::cout << "d.norm()/ld.cartesian().norm()=" << d.norm() / ld.norm() << std::endl;
            throw std::runtime_error("Rational Lattice DirectionType NOT FOUND\n");
        }
        return rld;
    }

    /**********************************************************************/
    template <int dim>
    LatticeVector<dim> Lattice<dim>::latticeVector(const VectorDimD &p) const
    {
        return LatticeVector<dim>(p, *this);
    }

    /**********************************************************************/
    template <int dim>
    ReciprocalLatticeVector<dim> Lattice<dim>::reciprocalLatticeVector(const VectorDimD &p) const
    {
        return ReciprocalLatticeVector<dim>(p, *this);
    }



    template class Lattice<1>;
    template class Lattice<2>;
    template class Lattice<3>;
}
#endif
