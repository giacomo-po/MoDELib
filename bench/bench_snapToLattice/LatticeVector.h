/* This file is part of PIL, the Particle Interaction Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LatticeVector_h_
#define model_LatticeVector_h_

#include <Eigen/Dense>
#include <model/Math/RoundEigen.h>

namespace model
{
    template <int dim>
    struct LatticeVector : public Eigen::Matrix<long int,dim,1>
    {
        static_assert(dim>0,"dim must be > 0.");
        
        typedef LatticeVector<dim> LatticeVectorType;
        typedef Eigen::Matrix<long int,dim,1> BaseType;
        typedef Eigen::Matrix<long int,dim,1> VectorDimI;
        typedef Eigen::Matrix<  double,dim,1> VectorDimD;
        
        //! The static column matrix of lattice vectors
        static Eigen::Matrix<double,dim,dim>    A;
        static Eigen::Matrix<double,dim,dim> invA;
        static Eigen::Matrix<double,dim,dim> cofA;
        
        static double roundTol;
        
        /**********************************************************************/
        static VectorDimI d2i(const VectorDimD& d)
        {
            VectorDimD nd(invA*d);
            VectorDimD rd(RoundEigen<double,dim>::round(nd));
            assert((nd-rd).norm()<roundTol && "Input vector is not a lattice vector");
            return rd.template cast<long int>();
        }
        
        /**********************************************************************/
        static int gcd(const size_t& a,const size_t& b)
        {
            return b>0? gcd(b, a % b) : a;
        }
        
        /**********************************************************************/
        static int gcd(const size_t& a,const size_t& b,const size_t& c)
        {
            return gcd(a,gcd(b,c));
        }
        
    public:
        
        /**********************************************************************/
        static void setLatticeMatrix(const Eigen::Matrix<double,dim,dim,1>& A_in)
        {
            
            A=A_in;
            invA=A.inverse();
            
        }
        
        LatticeVector() : BaseType(VectorDimI::Zero())
        {
            
        }
        
        LatticeVector(const VectorDimD& d) : BaseType(d2i(d))
        {
            
        }
        
        // This constructor allows to construct LatticeVector from Eigen expressions
        template<typename OtherDerived>
        LatticeVector(const Eigen::MatrixBase<OtherDerived>& other)
        : BaseType(other)
        { }
        
        // This method allows to assign Eigen expressions to LatticeVector
        template<typename OtherDerived>
        LatticeVector& operator= (const Eigen::MatrixBase<OtherDerived>& other)
        {
            this->BaseType::operator=(other);
            return *this;
        }
        
        LatticeVector cross(const LatticeVector& other) const
        {
            std::cout<<cofA * (dynamic_cast<const BaseType*>(this)->cross(other)).template cast<double>();
//            return LatticeVector(cofA * (dynamic_cast<const BaseType*>(this)->cross(other)).template cast<double>());
            return LatticeVector(dynamic_cast<const BaseType*>(this)->cross(other));
        }
        
//        LatticeVector asDirection() const
//        {
//        
//        }
        
    };
    
    template <int dim>
    Eigen::Matrix<double,dim,dim> LatticeVector<dim>::A=Eigen::Matrix<double,dim,dim,1>::Identity();
    
    template <int dim>
    Eigen::Matrix<double,dim,dim> LatticeVector<dim>::invA=A.inverse();
    
    template <int dim>
    Eigen::Matrix<double,dim,dim> LatticeVector<dim>::cofA=invA*A.determinant();
    
    
    template <int dim>
    double LatticeVector<dim>::roundTol=FLT_EPSILON;
    
} // end namespace
#endif
