/* This file is part of PIL, the Particle Interaction Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_Lattice_cpp_
#define model_Lattice_cpp_

#include <iomanip>
#include <cfloat> // FLT_EPSILON
#include <Eigen/Dense>

#include <LatticeMath.h>

namespace model
{
    /**********************************************************************/
    template <int dim>
    typename Lattice<dim>::MatrixDimD Lattice<dim>::getLatticeBasis(const typename Lattice<dim>::MatrixDimD& A,const typename Lattice<dim>::MatrixDimD& Q)
    {
        
        // Check that Q is orthogonal
        const typename Lattice<dim>::MatrixDimD QQT(Q*Q.transpose());
        if((QQT-typename Lattice<dim>::MatrixDimD::Identity()).norm()>2.0*DBL_EPSILON*dim*dim)
        {
            std::cout<<"Q=\n"<<Q<<std::endl;
            std::cout<<"Q*Q^T=\n"<<QQT<<std::endl;
            assert(false && "ROTATION MATRIX IS NOT ORTHOGONAL.");
        }
        
        // Check sure that C2G is proper
        assert(std::fabs(Q.determinant()-1.0) < FLT_EPSILON && "ROTATION MATRIX IS NOT PROPER.");
        
        // Check that A is full rank
        assert(std::fabs(A.determinant())>FLT_EPSILON && "A matrix is singular");
        
        return Q*A;
    }
    
    /**********************************************************************/
    template <int dim>
    Lattice<dim>::Lattice(const typename Lattice<dim>::MatrixDimD& A,const typename Lattice<dim>::MatrixDimD& Q) :
    //        /* init */ latticeBases(getLatticeBases(A,Q))
    /* init */ latticeBasis(getLatticeBasis(A,Q))
    /* init */,reciprocalBasis(latticeBasis.inverse().transpose())
    /* init */,C2G(Q)
    {
        
    }
    
    /**********************************************************************/
    template <int dim>
    Eigen::Matrix<long int,dim,1> Lattice<dim>::rationalApproximation(VectorDimD nd)
    {
        Eigen::Array<long int,dim,1> nums=Eigen::Matrix<long int,dim,1>::Zero();
        
        if(nd.squaredNorm()>0.0)
        {
            const Eigen::Array<double,dim,1> nda(nd.array().abs()); // vector of close-to-integer numbers corresponding to lattice coordinates
            size_t maxID=0;
            const double maxVal(nda.maxCoeff(&maxID));
            nd/=maxVal; // make each value of nd in [-1:1]
            
            nums=Eigen::Matrix<long int,dim,1>::Ones();
            Eigen::Array<long int,dim,1> dens=Eigen::Matrix<long int,dim,1>::Ones();
            long int denProd=1;
            
            for(int k=0;k<dim;++k)
            {
                BestRationalApproximation bra(nd(k),10000);
                
                nums(k)=bra.num;
                dens(k)=bra.den;
                denProd*=bra.den;
            }
            
            for(int k=0;k<dim;++k)
            {
                nums(k)*=(denProd/dens(k));
            }
        }
        
        return nums.matrix();
    }
    
    /**********************************************************************/
    template <int dim>
    typename Lattice<dim>:: LatticeVectorType Lattice<dim>::snapToLattice(const VectorDimD& d) const
    {
        VectorDimD nd(reciprocalBasis.transpose()*d);
        return typename Lattice<dim>:: LatticeVectorType(nd.array().round().matrix().template cast<long int>(),*this);
    }
    
    /**********************************************************************/
    template <int dim>
    typename Lattice<dim>:: LatticeDirectionType Lattice<dim>::latticeDirection(const VectorDimD& d) const
    {
        
        const VectorDimD nd(reciprocalBasis.transpose()*d);
        const typename Lattice<dim>:: LatticeVectorType temp(rationalApproximation(nd),*this);
        
        const double crossNorm(temp.cartesian().normalized().cross(d.normalized()).norm());
        if(crossNorm>FLT_EPSILON)
        {
            std::cout<<"input direction="<<d.normalized().transpose()<<std::endl;
            std::cout<<"lattice direction="<<temp.cartesian().normalized().transpose()<<std::endl;
            std::cout<<"cross product norm="<<std::setprecision(15)<<std::scientific<<crossNorm<<std::endl;
            assert(0 && "LATTICE DIRECTION NOT FOUND");
        }
        
        return typename Lattice<dim>:: LatticeDirectionType(temp);
    }
    
    /**********************************************************************/
    template <int dim>
    typename Lattice<dim>::ReciprocalLatticeDirectionType Lattice<dim>::reciprocalLatticeDirection(const VectorDimD& d) const
    {
        
        const VectorDimD nd(latticeBasis.transpose()*d);
        const typename Lattice<dim>::ReciprocalLatticeVectorType temp(rationalApproximation(nd),*this);
        
        const double crossNorm(temp.cartesian().normalized().cross(d.normalized()).norm());
        if(crossNorm>FLT_EPSILON)
        {
            std::cout<<"input direction="<<std::setprecision(15)<<std::scientific<<d.normalized().transpose()<<std::endl;
            std::cout<<"reciprocal lattice direction="<<std::setprecision(15)<<std::scientific<<temp.cartesian().normalized().transpose()<<std::endl;
            std::cout<<"cross product norm="<<std::setprecision(15)<<std::scientific<<crossNorm<<std::endl;
            assert(0 && "RECIPROCAL LATTICE DIRECTION NOT FOUND");
        }
        
        return typename Lattice<dim>::ReciprocalLatticeDirectionType(temp);
    }
    
    /**********************************************************************/
    template <int dim>
    typename Lattice<dim>::RationalLatticeDirectionType Lattice<dim>::rationalLatticeDirection(const VectorDimD& d,
                                                                        const typename BestRationalApproximation::LongIntType& maxDen) const
    {
        
        const typename Lattice<dim>:: LatticeDirectionType ld(latticeDirection(d));
        const BestRationalApproximation bra(d.norm()/ld.cartesian().norm(),maxDen);
        Rational rat(bra.num,bra.den);
        typename Lattice<dim>:: RationalLatticeDirectionType rld(rat,ld);
        if((rld.cartesian()-d).squaredNorm()>FLT_EPSILON)
        {
            std::cout<<"input vector="<<d.transpose()<<std::endl;
            std::cout<<"lattice direction="<<ld.cartesian().transpose()<<std::endl;
            std::cout<<"rational="<<rat<<std::endl;
            std::cout<<"d.norm()/ld.cartesian().norm()="<<d.norm()/ld.norm()<<std::endl;
            assert(0 && "typename Lattice<dim>:: RationalLatticeDirectionType NOT FOUND");
        }
        return rld;
        
    }
    
    /**********************************************************************/
    template <int dim>
    typename Lattice<dim>:: LatticeVectorType Lattice<dim>::latticeVector(const VectorDimD& p) const
    {
        return typename Lattice<dim>:: LatticeVectorType(p,*this);
    }
    
    /**********************************************************************/
    template <int dim>
    typename Lattice<dim>::ReciprocalLatticeVectorType Lattice<dim>::reciprocalLatticeVector(const VectorDimD& p) const
    {
        return typename Lattice<dim>::ReciprocalLatticeVectorType(p,*this);
    }
}
#endif
