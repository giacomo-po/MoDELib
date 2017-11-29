/* This file is part of PIL, the Particle Interaction Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_Lattice_h_
#define model_Lattice_h_

#include <cfloat> // FLT_EPSILON
#include <Eigen/Dense>
#include <model/Utilities/StaticID.h>
#include <model/Math/RoundEigen.h>
//#include <model/LatticeMath/LatticeBase.h>
//#include <model/LatticeMath/ReciprocalLatticeVector.h>
#include <model/LatticeMath/LatticeVector.h>
#include <model/LatticeMath/ReciprocalLatticeVector.h>
#include <model/LatticeMath/LatticeDirection.h>
#include <model/LatticeMath/ReciprocalLatticeDirection.h>
#include <model/Math/BestRationalApproximation.h>


namespace model
{

//    /**************************************************************************/
//    /**************************************************************************/
//    template <int dim,int k>
//    class LatticeBase : public StaticID<LatticeBase<dim,k>>
//    {
//        static_assert(dim>0,"dim must be > 0.");
//        static_assert(k>0,"k must be > 0.");
//        static_assert(k<=dim,"k must be <= dim.");
//        
//        typedef Eigen::Matrix<double,dim,dim> MatrixDimK;
//
//        MatrixDimK    _covBasis;
//        
//    public:
//        
//        /**********************************************************************/
//        LatticeBase() :
//        /* init */ _covBasis(MatrixDimK::Identity())
//        {
//            
//        }
//        
//        /**********************************************************************/
//        LatticeBase(const MatrixDimK& cov) :
//        /* init */ _covBasis(cov)
//        {
//            
//        }
//        
//        /**********************************************************************/
//        void setLatticeBasis(const MatrixDimK& A)
//        {
//            _covBasis=A;
//            std::cout<<"Lattice basis (in columns) =\n"<<_covBasis<<std::endl;
//        }
//        
//        /**********************************************************************/
//        const MatrixDimK& covBasis() const
//        {
//            return _covBasis;
//        }
//    };
//    
//    /**************************************************************************/
//    /**************************************************************************/
//    template <int dim,int k>
//    class Lattice : public LatticeBase<dim,k>
//    {
//
//        public:
//    };

    /**************************************************************************/
    /**************************************************************************/
    template <int dim>
//    class Lattice<dim,dim> : public LatticeBase<dim,dim>
    class Lattice
    {
        static constexpr double roundTol=FLT_EPSILON;
        typedef Eigen::Matrix<  double,dim,1> VectorDimD;
        typedef Eigen::Matrix<long int,dim,1> VectorDimI;
        typedef Eigen::Matrix<double,dim,dim> MatrixDimD;
        typedef LatticeVector<dim> LatticeVectorType;
        typedef LatticeDirection<dim> LatticeDirectionType;
        typedef ReciprocalLatticeVector<dim> ReciprocalLatticeVectorType;
        typedef ReciprocalLatticeDirection<dim> ReciprocalLatticeDirectionType;
//        typedef LatticeBase<dim,dim> LatticeBaseType;
        
        //! The static column matrix of lattice vectors
 
        /**********************************************************************/
        std::pair<MatrixDimD,MatrixDimD> getLatticeBases(const MatrixDimD& A)
        {
            assert(std::fabs(A.determinant())>FLT_EPSILON && "A matrix is singular");
            
            return std::pair<MatrixDimD,MatrixDimD>(A,A.inverse().transpose());
        }

        
        std::pair<MatrixDimD,MatrixDimD> latticeBases;
        MatrixDimD&    _covBasis;
        MatrixDimD& _contraBasis;

        /**********************************************************************/
        Lattice(const MatrixDimD& A,const MatrixDimD& invAT) :
        /* init */ latticeBases(std::make_pair(A,invAT)),
        /* init */ _covBasis(latticeBases.first),
        /* init */ _contraBasis(latticeBases.second)
        {
            
        }
        
    public:
        
        /**********************************************************************/
        Lattice() :
        /* init */ latticeBases(getLatticeBases(MatrixDimD::Identity())),
        /* init */ _covBasis(latticeBases.first),
        /* init */ _contraBasis(latticeBases.second)
        {

        }
        
        /**********************************************************************/
        Lattice(const MatrixDimD& A) :
        /* init */ latticeBases(getLatticeBases(A)),
        /* init */ _covBasis(latticeBases.first),
        /* init */ _contraBasis(latticeBases.second)
        {
            
        }
        
        /**********************************************************************/
        void setLatticeBasis(const MatrixDimD& A)
        {
            latticeBases=getLatticeBases(A);
//            _covBasis=A;
//            LatticeBaseType::setLatticeBasis(A);
//            _contraBasis=covBasis().inverse().transpose();
//            _contraBasis=_covBasis.inverse().transpose();
//            std::cout<<"Lattice basis (in columns) =\n"<<_covBasis<<std::endl;
//            std::cout<<"Lattice reciprocal basis (in columns) =\n"<<_contraBasis<<std::endl;
        }
        
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
            VectorDimD nd(_contraBasis.transpose()*d);
            return LatticeVectorType(RoundEigen<double,dim>::round(nd).template cast<long int>(),*this);
        }
        
        /**********************************************************************/
        LatticeDirectionType latticeDirection(const VectorDimD& d) const
        {
            
            const VectorDimD nd(_contraBasis.transpose()*d);
//            const LatticeVectorType temp(rationalApproximation(nd),covBasis(),contraBasis());
            const LatticeVectorType temp(rationalApproximation(nd),*this);
            
            if(temp.cartesian().normalized().cross(d.normalized()).norm()>FLT_EPSILON)
            {
                model::cout<<"input direction="<<d.normalized().transpose()<<std::endl;
                model::cout<<"lattice direction="<<temp.cartesian().normalized().transpose()<<std::endl;
                assert(0 && "LATTICE DIRECTION NOT FOUND");
            }
            
            return LatticeDirectionType(temp);
        }
      
        /**********************************************************************/
        ReciprocalLatticeDirectionType reciprocalLatticeDirection(const VectorDimD& d) const
        {
            
            const VectorDimD nd(covBasis().transpose()*d);
            const ReciprocalLatticeVectorType temp(rationalApproximation(nd),*this);
            
            if(temp.cartesian().normalized().cross(d.normalized()).norm()>FLT_EPSILON)
            {
                model::cout<<"input direction="<<d.normalized().transpose()<<std::endl;
                model::cout<<"reciprocal lattice direction="<<temp.cartesian().normalized().transpose()<<std::endl;
                assert(0 && "RECIPROCAL LATTICE DIRECTION NOT FOUND");
            }
            
            return ReciprocalLatticeDirectionType(temp);
        }
        
        /**********************************************************************/
        const MatrixDimD& covBasis() const
        {
            return _covBasis;
        }
        
        const VectorDimD basisVector(const size_t& c) const
        {
            assert(c<dim);
            return covBasis().col(c);
        }
        
        /**********************************************************************/
        const MatrixDimD& contraBasis() const
        {
            return _contraBasis;
        }
        
        /**********************************************************************/
        LatticeVectorType latticeVector(const VectorDimD& p) const
        {
//            return LatticeVectorType(p,covBasis(),contraBasis());
            return LatticeVectorType(p,*this);

        }
        
        /**********************************************************************/
        ReciprocalLatticeVectorType reciprocalLatticeVector(const VectorDimD& p) const
        {
//            return ReciprocalLatticeVectorType(p,covBasis(),contraBasis());
            return ReciprocalLatticeVectorType(p,*this);

        }
        
//        /**********************************************************************/
//        Lattice<dim,dim> reciprocal() const
//        {
//            return Lattice<dim,dim>(contraBasis(),covBasis());
//        }
        
    };

    
} // end namespace
#endif


