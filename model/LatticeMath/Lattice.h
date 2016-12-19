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
    
    template <int dim>
    class Lattice
    {
        static_assert(dim>0,"dim must be > 0.");
        static constexpr double roundTol=FLT_EPSILON;
        typedef Eigen::Matrix<  double,dim,1> VectorDimD;
        typedef Eigen::Matrix<long int,dim,1> VectorDimI;
        typedef Eigen::Matrix<double,dim,dim> MatrixDimD;
        typedef LatticeVector<dim> LatticeVectorType;
        typedef LatticeDirection<dim> LatticeDirectionType;
        typedef ReciprocalLatticeVector<dim> ReciprocalLatticeVectorType;
        typedef ReciprocalLatticeDirection<dim> ReciprocalLatticeDirectionType;

        
        //! The static column matrix of lattice vectors
        MatrixDimD    _covBasis;
        MatrixDimD _contraBasis;

        
    public:
        

//        Eigen::Matrix<double,dim,dim> C2G;

        /**********************************************************************/
        Lattice() :
        /* init */ _covBasis(MatrixDimD::Identity()),
        /* init */ _contraBasis(MatrixDimD::Identity())
 //       /* init */ C2G(Eigen::Matrix<double,dim,dim>::Identity())
        {
//            model::cout<<"Creating Grain "<<grainID<<std::endl;
//            selectMaterial(materialZ);
        }
        
        /**********************************************************************/
        void setLatticeBasis(const MatrixDimD& A)
        {
            _covBasis=A;
            _contraBasis=_covBasis.inverse().transpose();
            std::cout<<"Lattice basis (in columns) =\n"<<_covBasis<<std::endl;
            std::cout<<"Lattice reciprocal basis (in columns) =\n"<<_contraBasis<<std::endl;
        }
        
        /**********************************************************************/
        static Eigen::Matrix<long int,dim,1> rationalApproximation(VectorDimD nd)
        {
            const Eigen::Array<double,dim,1> nda(nd.array().abs()); // vector of close-to-integer numbers corresponding to lattice coordinates
            size_t maxID=0;
            const double maxVal(nda.maxCoeff(&maxID));
            nd/=maxVal; // make each value of nd in [-1:1]
            
            Eigen::Array<long int,dim,1> nums=Eigen::Matrix<long int,dim,1>::Ones();
            Eigen::Array<long int,dim,1> dens=Eigen::Matrix<long int,dim,1>::Ones();
            long int denProd=1;
            
            for(int k=0;k<dim;++k)
            {
                BestRationalApproximation bra(nd(k),100);
                
                nums(k)=bra.num;
                dens(k)=bra.den;
                denProd*=bra.den;
            }
            
            for(int k=0;k<dim;++k)
            {
                nums(k)*=(denProd/dens(k));
            }
            return nums.matrix();
        }
        
        /**********************************************************************/
        LatticeVectorType snapToLattice(const VectorDimD& d) const
        {
            VectorDimD nd(_contraBasis.transpose()*d);
            return LatticeVectorType(RoundEigen<double,dim>::round(nd).template cast<long int>(),_covBasis,_contraBasis);
        }
        
        /**********************************************************************/
        LatticeDirectionType latticeDirection(const VectorDimD& d) const
        {
            
            const VectorDimD nd(_contraBasis.transpose()*d);
            const LatticeVectorType temp(rationalApproximation(nd),_covBasis,_contraBasis);
            
            if(temp.cartesian().normalized().cross(d.normalized()).norm()>FLT_EPSILON)
            {
                std::cout<<"input direction="<<d.normalized().transpose()<<std::endl;
                std::cout<<"lattice direction="<<temp.cartesian().normalized().transpose()<<std::endl;
                assert(0 && "LATTICE DIRECTION NOT FOUND");
            }
            
            return LatticeDirectionType(temp);
        }
      
        /**********************************************************************/
        ReciprocalLatticeDirectionType reciprocalLatticeDirection(const VectorDimD& d) const
        {
            
            const VectorDimD nd(_covBasis.transpose()*d);
            const ReciprocalLatticeVectorType temp(rationalApproximation(nd),_covBasis,_contraBasis);
            
            if(temp.cartesian().normalized().cross(d.normalized()).norm()>FLT_EPSILON)
            {
                std::cout<<"input direction="<<d.normalized().transpose()<<std::endl;
                std::cout<<"reciprocal lattice direction="<<temp.cartesian().normalized().transpose()<<std::endl;
                assert(0 && "RECIPROCAL LATTICE DIRECTION NOT FOUND");
            }
            
            return ReciprocalLatticeDirectionType(temp);
        }
        
        /**********************************************************************/
        const MatrixDimD& covBasis() const
        {
            return _covBasis;
        }
        
        /**********************************************************************/
        const MatrixDimD& contraBasis() const
        {
            return _contraBasis;
        }
        
        /**********************************************************************/
        LatticeVectorType latticeVector(const VectorDimD& p) const
        {
            return LatticeVectorType(p,_covBasis,_contraBasis);
        }
        
        /**********************************************************************/
        ReciprocalLatticeVectorType reciprocalLatticeVector(const VectorDimD& p) const
        {
            return ReciprocalLatticeVectorType(p,_covBasis,_contraBasis);
        }
        
    };

    
} // end namespace
#endif


//        /**********************************************************************/
//        LatticeVector(const LatticeVectorType& other) :
//        /* base copy */ BaseType(static_cast<VectorDimI>(other))
//        {
//            assert(&covBasis==&other.covBasis && "LatticeVectors have different bases.");
//            assert(&contraBasis==&other.contraBasis && "LatticeVectors have different bases.");
//        }
//
//        /**********************************************************************/
//        LatticeVector(LatticeVectorType&& other) :
//        /* base copy */ BaseType(std::move(static_cast<VectorDimI>(other)))
//        {
//            assert(&covBasis==&other.covBasis && "LatticeVectors have different bases.");
//            assert(&contraBasis==&other.contraBasis && "LatticeVectors have different bases.");
//        }


//        template<typename OtherDerived>
//        LatticeVector(const Eigen::MatrixBase<OtherDerived>& other) :
//        /* base init */ BaseType(other)
//        {// This constructor allows  to construct LatticeVector from Eigen expressions
//
//
//        }

//        template<typename OtherDerived>
//        LatticeVector& operator=(const Eigen::MatrixBase <OtherDerived>& other)
//        {// This method allows to assign Eigen expressions to LatticeVector
//            BaseType::operator=(other);
//            return *this;
//        }

//        VectorDimD cartesian() const
//        {
//            return LatticeBaseType::covBasis()*this->template cast<double>();
//        }
