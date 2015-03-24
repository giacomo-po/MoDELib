/* This file is part of PIL, the Particle Interaction Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LatticeVector_h_
#define model_LatticeVector_h_

#include <cfloat> // FLT_EPSILON
#include <Eigen/Dense>
#include <model/Math/RoundEigen.h>
//#include <model/DislocationDynamics/Materials/CrystalOrientation.h>
#include <model/DislocationDynamics/LatticeMath/LatticeBase.h>
#include <model/DislocationDynamics/LatticeMath/ContravariantCoordinate.h>


namespace model
{
    template <int dim>
    struct LatticeVector //: private ContraVariantVector<long int,dim,1>,
                        //  private     CoVariantVector<long int,dim,1>
    //Eigen::Matrix<long int,dim,1>
    {
        static_assert(dim>0,"dim must be > 0.");
        
        typedef LatticeBase<dim> LatticeBaseType;
        typedef LatticeVector<dim> LatticeVectorType;
        typedef Eigen::Matrix<long int,dim,1> BaseType;
        typedef Eigen::Matrix<long int,dim,1> VectorDimI;
        typedef Eigen::Matrix<  double,dim,1> VectorDimD;
        typedef ContravariantCoordinate<long int, dim> ContravariantCoordinateType;
        
//        //! The static column matrix of lattice vectors
//        static Eigen::Matrix<double,dim,dim>    A;
//        static Eigen::Matrix<double,dim,dim> invA;
//        static Eigen::Matrix<double,dim,dim> cofA;
//        
//        static double roundTol;
    private:
        /**********************************************************************/
        static VectorDimI d2contra(const VectorDimD& d)
        {
            VectorDimD nd(LatticeBaseType::invCovBasis()*d);
            VectorDimD rd(RoundEigen<double,dim>::round(nd));
            std::cout<<"d2contra, nd="<<nd.transpose()<<std::endl;
            std::cout<<"d2contra, rd="<<rd.transpose()<<std::endl;
            assert((nd-rd).norm()<roundTol && "Input vector is not a lattice vector");
            return rd.template cast<long int>();
        }
        
        /**********************************************************************/
        static VectorDimI d2cov(const VectorDimD& d)
        {
            VectorDimD nd(LatticeBaseType::invContraBasis()*d);
            VectorDimD rd(RoundEigen<double,dim>::round(nd));
            std::cout<<"d2cov, nd="<<nd.transpose()<<std::endl;
            std::cout<<"d2cov, rd="<<rd.transpose()<<std::endl;
            assert((nd-rd).norm()<roundTol && "Input vector is not a lattice vector");
            return rd.template cast<long int>();
        }

        VectorDimI contraVariant; // contravariant components
        VectorDimI coVariant;     //     covariant components
        
    public:
        
        static double roundTol;

        
        // Empty constructor
        LatticeVector() : //BaseType(VectorDimI::Zero())
        contraVariant(VectorDimI::Zero()),
            coVariant(VectorDimI::Zero())
        {
            
        }
        
        LatticeVector(const VectorDimD& d) :
        //BaseType(d2i(d))
        /* init */ contraVariant(d2contra(d)),
        /* init */ coVariant(d2cov(d))
        {/*!@param[in] d vector in real space
          * Constructs *this by spapping d to the lattice
          */
            
        }
        
        LatticeVector(const ContravariantCoordinateType& cont) :
        /* init */ contraVariant(cont),
        /* init */ coVariant(???)
        {
        
        }
        
//        LatticeVector(const LatticeVector& n1,const LatticeVector& n2) :
//        BaseType(d2i(cofA.transpose()*(gcd(n1.cross(n2))).template cast<double>()))
//        {/*!@param[in] n1 Lattice vector defiing the first plane
//          *!@param[in] n2 Lattice vector defiing the scond plane
//          * Constructs *this as the lattice vector direction at the intersection of the two planes
//          */
//            
//        }
//        
//        // This constructor allows to construct LatticeVector from Eigen expressions
//        template<typename OtherDerived>
//        LatticeVector(const Eigen::MatrixBase<OtherDerived>& other)
//        : BaseType(other)
//        { }
//        
//        // This method allows to assign Eigen expressions to LatticeVector
//        template<typename OtherDerived>
//        LatticeVector& operator= (const Eigen::MatrixBase<OtherDerived>& other)
//        {
//            this->BaseType::operator=(other);
//            return *this;
//        }
        
//        VectorDimD global() const
//        {
//            return A*(*this);
//        }
        

        const VectorDimI& contra() const
        {
            return contraVariant;
        }
        
        const VectorDimI& cov() const
        {
            return coVariant;
        }
        
        long int dot(const LatticeVectorType& other) const
        {
            return contra().dot(other.cov());
        }
        
        
    };
    
    
    template <int dim>
    double LatticeVector<dim>::roundTol=FLT_EPSILON;
    
//    template <int dim>
//    Eigen::Matrix<double,dim,dim> LatticeVector<dim>::A=Eigen::Matrix<double,dim,dim,1>::Identity();
//    
//    template <int dim>
//    Eigen::Matrix<double,dim,dim> LatticeVector<dim>::invA=A.inverse();
//    
//    template <int dim>
//    Eigen::Matrix<double,dim,dim> LatticeVector<dim>::cofA=invA*A.determinant();
//    
//    
//    template <int dim>
//    double LatticeVector<dim>::roundTol=FLT_EPSILON;
    
} // end namespace
#endif
