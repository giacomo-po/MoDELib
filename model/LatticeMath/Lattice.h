/* This file is part of PIL, the Particle Interaction Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_Lattice_h_
#define model_Lattice_h_

#include <iomanip> // FLT_EPSILON
#include <cfloat> // FLT_EPSILON
#include <Eigen/Dense>
#include <StaticID.h>
//#include <RoundEigen.h>
//#include <LatticeBase.h>
//#include <ReciprocalLatticeVector.h>
#include <LatticeVector.h>
#include <ReciprocalLatticeVector.h>
#include <LatticeDirection.h>
#include <ReciprocalLatticeDirection.h>
#include <BestRationalApproximation.h>

#include <RationalLatticeDirection.h>

namespace model
{

    /**************************************************************************/
    /**************************************************************************/
    template <int dim>
    class Lattice
    {
        static constexpr double roundTol=FLT_EPSILON;
        typedef Eigen::Matrix<  double,dim,1> VectorDimD;
        typedef Eigen::Matrix<long int,dim,1> VectorDimI;
        typedef Eigen::Matrix<double,dim,dim> MatrixDimD;
        typedef LatticeVector<dim> LatticeVectorType;
        typedef LatticeDirection<dim> LatticeDirectionType;
        typedef RationalLatticeDirection<dim> RationalLatticeDirectionType;
        typedef ReciprocalLatticeVector<dim> ReciprocalLatticeVectorType;
        typedef ReciprocalLatticeDirection<dim> ReciprocalLatticeDirectionType;

//        typedef LatticeBase<dim,dim> LatticeBaseType;
        
        //! The static column matrix of lattice vectors
 
        /**********************************************************************/
        std::tuple<MatrixDimD,MatrixDimD,MatrixDimD> getLatticeBases(const MatrixDimD& A,const MatrixDimD& Q)
        {

            // Check that Q is orthogonal
            const MatrixDimD QQT(Q*Q.transpose());
            if((QQT-MatrixDimD::Identity()).norm()>2.0*DBL_EPSILON*dim*dim)
            {
                std::cout<<"Q=\n"<<Q<<std::endl;
                std::cout<<"Q*Q^T=\n"<<QQT<<std::endl;
                assert(false && "ROTATION MATRIX IS NOT ORTHOGONAL.");
            }
            
            // Check sure that C2G is proper
            assert(std::fabs(Q.determinant()-1.0) < FLT_EPSILON && "ROTATION MATRIX IS NOT PROPER.");

            // Check that A is full rank
            assert(std::fabs(A.determinant())>FLT_EPSILON && "A matrix is singular");
            
            const MatrixDimD QA(Q*A);
            return std::make_tuple(QA,QA.inverse().transpose(),Q);
        }

        

//        /**********************************************************************/
//        Lattice(const MatrixDimD& A,const MatrixDimD& invAT) :
//        /* init */ latticeBases(std::make_pair(A,invAT)),
//        /* init */ latticeBasis(latticeBases.first),
//        /* init */ reciprocalBasis(latticeBases.second)
//        {
//            
//        }
        
        std::tuple<MatrixDimD,MatrixDimD,MatrixDimD> latticeBases;

        
    public:
        
        const MatrixDimD&    latticeBasis;
        const MatrixDimD& reciprocalBasis;
        const MatrixDimD& C2G;
        
        /**********************************************************************/
        Lattice(const MatrixDimD& A,const MatrixDimD& Q) :
        /* init */ latticeBases(getLatticeBases(A,Q))
        /* init */,latticeBasis(std::get<0>(latticeBases))
        /* init */,reciprocalBasis(std::get<1>(latticeBases))
        /* init */,C2G(std::get<2>(latticeBases))
        {

        }
        
        /**********************************************************************/
        Lattice(const Lattice& other) :
        /* init */ latticeBases(other.latticeBases)
        /* init */,latticeBasis(std::get<0>(latticeBases))
        /* init */,reciprocalBasis(std::get<1>(latticeBases))
        /* init */,C2G(std::get<2>(latticeBases))
        {/*!The copy contructor initializes latticeBases to other.latticeBases,
          * and then the references latticeBasis and reciprocalBasis to the local matrices
          * Note that this allows to wite:
          * Lattice L1(A); // ok
          * Lattice L2(L1); // ok, copy contructor
          * Lattice L3=L1;  // ok, assignemt operator here works as copy contructor
          * Lattice L4(A);  
          * L4=L1;          // ERROR: const references latticeBasis and reciprocalBasis cannot be assigned
          */
            
        }
        
//        /**********************************************************************/
//        void rotate(const MatrixDimD& Q)
//        {/*! Z is atomic number
//          */
//
//            const MatrixDimD QQT(Q*Q.transpose());
//            if((QQT-MatrixDimD::Identity()).norm()>2.0*DBL_EPSILON*dim*dim)
//            {
//                std::cout<<"Q=\n"<<Q<<std::endl;
//                std::cout<<"Q*Q^T=\n"<<QQT<<std::endl;
//                assert(false && "ROTATION MATRIX IS NOT ORTHOGONAL.");
//            }
//            // make sure that C2G is proper
//            assert(std::fabs(Q.determinant()-1.0) < FLT_EPSILON && "ROTATION MATRIX IS NOT PROPER.");
////            setLatticeBasis();
//            latticeBases=getLatticeBases((Q*latticeBasis).eval());
//
//        }

        /**********************************************************************/
        static Eigen::Matrix<long int,dim,1> rationalApproximation(VectorDimD nd)
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
        LatticeVectorType snapToLattice(const VectorDimD& d) const
        {
            VectorDimD nd(reciprocalBasis.transpose()*d);
            return LatticeVectorType(nd.array().round().matrix().template cast<long int>(),*this);
        }
        
        /**********************************************************************/
        LatticeDirectionType latticeDirection(const VectorDimD& d) const
        {
            
            const VectorDimD nd(reciprocalBasis.transpose()*d);
            const LatticeVectorType temp(rationalApproximation(nd),*this);
            
            const double crossNorm(temp.cartesian().normalized().cross(d.normalized()).norm());
            if(crossNorm>FLT_EPSILON)
            {
                model::cout<<"input direction="<<d.normalized().transpose()<<std::endl;
                model::cout<<"lattice direction="<<temp.cartesian().normalized().transpose()<<std::endl;
                model::cout<<"cross product norm="<<std::setprecision(15)<<std::scientific<<crossNorm<<std::endl;
                assert(0 && "LATTICE DIRECTION NOT FOUND");
            }
            
            return LatticeDirectionType(temp);
        }
      
        /**********************************************************************/
        ReciprocalLatticeDirectionType reciprocalLatticeDirection(const VectorDimD& d) const
        {
            
            const VectorDimD nd(latticeBasis.transpose()*d);
            const ReciprocalLatticeVectorType temp(rationalApproximation(nd),*this);
            
            const double crossNorm(temp.cartesian().normalized().cross(d.normalized()).norm());
            if(crossNorm>FLT_EPSILON)
            {
                model::cout<<"input direction="<<std::setprecision(15)<<std::scientific<<d.normalized().transpose()<<std::endl;
                model::cout<<"reciprocal lattice direction="<<std::setprecision(15)<<std::scientific<<temp.cartesian().normalized().transpose()<<std::endl;
                model::cout<<"cross product norm="<<std::setprecision(15)<<std::scientific<<crossNorm<<std::endl;
                assert(0 && "RECIPROCAL LATTICE DIRECTION NOT FOUND");
            }
            
            return ReciprocalLatticeDirectionType(temp);
        }
        
        /**********************************************************************/
        RationalLatticeDirectionType rationalLatticeDirection(const VectorDimD& d,
                                                              const typename BestRationalApproximation::LongIntType& maxDen=1000) const
        {
            
            const LatticeDirectionType ld(latticeDirection(d));
            const BestRationalApproximation bra(d.norm()/ld.cartesian().norm(),maxDen);
            Rational rat(bra.num,bra.den);
            RationalLatticeDirectionType rld(rat,ld);
            if((rld.cartesian()-d).squaredNorm()>FLT_EPSILON)
            {
                model::cout<<"input vector="<<d.transpose()<<std::endl;
                model::cout<<"lattice direction="<<ld.cartesian().transpose()<<std::endl;
                model::cout<<"rational="<<rat<<std::endl;
                model::cout<<"d.norm()/ld.cartesian().norm()="<<d.norm()/ld.norm()<<std::endl;
                assert(0 && "RationalLatticeDirectionType NOT FOUND");
            }
            return rld;

        }
        
        /**********************************************************************/
        LatticeVectorType latticeVector(const VectorDimD& p) const
        {
            return LatticeVectorType(p,*this);
        }
        
        /**********************************************************************/
        ReciprocalLatticeVectorType reciprocalLatticeVector(const VectorDimD& p) const
        {
            return ReciprocalLatticeVectorType(p,*this);
        }
        
    };

    
}
#endif

//        /**********************************************************************/
//        Lattice<dim,dim> reciprocal() const
//        {
//            return Lattice<dim,dim>(reciprocalBasis(),latticeBasis());
//        }

//        /**********************************************************************/
//        const MatrixDimD& latticeBasis() const __attribute__ ((deprecated)) // make latticeBasis public since it is const
//        {
//            return latticeBasis;
//        }
//
//        /**********************************************************************/
//        const MatrixDimD& reciprocalBasis() const __attribute__ ((deprecated)) // make reciprocalBasis public since it is const
//        {
//            return reciprocalBasis;
//        }
//
//        const VectorDimD basisVector(const size_t& c) const
//        {
//            assert(c<dim);
//            return latticeBasis().col(c);
//        }

//        /**********************************************************************/
//        void setLatticeBasis(const MatrixDimD& A)
//        {
//            latticeBases=getLatticeBases(A);
////            latticeBasis=A;
////            LatticeBaseType::setLatticeBasis(A);
////            reciprocalBasis=latticeBasis().inverse().transpose();
////            reciprocalBasis=latticeBasis.inverse().transpose();
////            std::cout<<"Lattice basis (in columns) =\n"<<latticeBasis<<std::endl;
////            std::cout<<"Lattice reciprocal basis (in columns) =\n"<<reciprocalBasis<<std::endl;
//        }

