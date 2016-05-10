/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LatticeBase_h_
#define model_LatticeBase_h_

#include <iostream>
#include <cfloat> // FLT_EPSILON
#include <Eigen/Dense>
#include <model/Math/RoundEigen.h>

namespace model
{
    
    template <int dim>
    class LatticeVector;
    
    template <int dim>
    class ReciprocalLatticeVector;
    
    class LatticePlane;
    class LatticeLine;
    
    
    template <int dim>
    class LatticeBase
    {
        static_assert(dim>0,"dim must be > 0.");
        
        typedef Eigen::Matrix<long int,dim,1> VectorDimI;
        typedef Eigen::Matrix<  double,dim,1> VectorDimD;
        typedef Eigen::Matrix<double,dim,dim> MatrixDimD;
        
        //! The static column matrix of lattice vectors
        static MatrixDimD    A;
        static MatrixDimD    AT;
        static MatrixDimD invA;
        static MatrixDimD invAT;
        //        static MatrixDimD cofA;
        
        //        static double roundTol;
        
        
    public:
        static const double roundTol;
        
        
        /**********************************************************************/
        static void setLatticeBasis(const Eigen::Matrix<double,dim,dim,1>& A_in)
        {
            A=A_in;
            AT=A.transpose();
            invA=A.inverse();
            invAT=invA.transpose();
            //            cofA=invA*A.determinant();
            
            std::cout<<"Lattice basis (in columns) =\n"<<A<<std::endl;
            std::cout<<"Lattice reciprocal basis (in columns) =\n"<<invAT<<std::endl;
            
        }
        
        /**********************************************************************/
        static VectorDimI d2contra(const VectorDimD& d)
        {
            const VectorDimD nd(invA*d);
            const VectorDimD rd(RoundEigen<double,dim>::round(nd));
            if((nd-rd).norm()>roundTol)
            {
                std::cout<<"d2contra, nd="<<nd.transpose()<<std::endl;
                std::cout<<"d2contra, rd="<<rd.transpose()<<std::endl;
                assert(0 && "Input vector is not a lattice vector");
            }
            return rd.template cast<long int>();
        }
        
        /**********************************************************************/
        static VectorDimI latticeDirection(const VectorDimD& d)
        {
            bool found=false;
            VectorDimD rdk(VectorDimD::Zero());

            const VectorDimD nd(invA*d);
            
            
            for(int k=0;k<dim;++k)
            {
                const VectorDimD ndk(nd/nd(k));
                rdk=RoundEigen<double,dim>::round(ndk);
                if((ndk-rdk).norm()<roundTol)
                {
                    found=true;
                    break;
                }

            }
            assert(found && "Input vector is not on a lattice direction");
            return rdk.template cast<long int>();
        }
        
        /**********************************************************************/
        static VectorDimI snapToLattice(const VectorDimD& d)
        {
            VectorDimD nd(invA*d);
            return RoundEigen<double,dim>::round(nd).template cast<long int>();
        }
        
        /**********************************************************************/
        static VectorDimI d2cov(const VectorDimD& d)
        {
            const VectorDimD nd(AT*d);
            const VectorDimD rd(RoundEigen<double,dim>::round(nd));
            if((nd-rd).norm()>roundTol)
            {
                std::cout<<"d2cov, nd="<<nd.transpose()<<std::endl;
                std::cout<<"d2cov, rd="<<rd.transpose()<<std::endl;
                assert(0 && "Input vector is not a reciprocal lattice vector");
            }
            //            assert((nd-rd).norm()<roundTol && "Input vector is not a lattice vector");
            return rd.template cast<long int>();
        }
        
        /**********************************************************************/
        static VectorDimI reciprocalLatticeDirection(const VectorDimD& d)
        {
            bool found=false;
            VectorDimD rdk(VectorDimD::Zero());
            
            const VectorDimD nd(AT*d);
            
            
            for(int k=0;k<dim;++k)
            {
                const VectorDimD ndk(nd/nd(k));
                rdk=RoundEigen<double,dim>::round(ndk);
                if((ndk-rdk).norm()<roundTol)
                {
                    found=true;
                    break;
                }
                
            }
            assert(found && "Input vector is not on a lattice direction");
            return rdk.template cast<long int>();
        }
        
        static const MatrixDimD& covBasis()
        {
            return A;
        }
        
        static const MatrixDimD& invCovBasis()
        {
            return invA;
        }
        
        static const MatrixDimD& contraBasis()
        {
            return invAT;
        }
        
        static const MatrixDimD& invContraBasis()
        {
            return AT;
        }
        
        //        /**********************************************************************/
        //        static int gcd(const size_t& a,const size_t& b)
        //        {
        //            return b>0? gcd(b, a % b) : a;
        //        }
        
        //        static int signOfFirstNonzero(const VectorDimI& v)
        //        {
        //            int sgn=1;
        //            for(int d=0;d<dim;++d)
        //            {
        //                if(v(d)!=0)
        //                {
        //                    sgn = ((v(d)>0)? 1 : -1);
        //                    break;
        //                }
        //            }
        //            return sgn;
        //        }
        
    };
    
    template <int dim>
    Eigen::Matrix<double,dim,dim> LatticeBase<dim>::A=Eigen::Matrix<double,dim,dim,1>::Identity();
    
    template <int dim>
    Eigen::Matrix<double,dim,dim> LatticeBase<dim>::AT=A.transpose();
    
    
    template <int dim>
    Eigen::Matrix<double,dim,dim> LatticeBase<dim>::invA=A.inverse();
    
    template <int dim>
    Eigen::Matrix<double,dim,dim> LatticeBase<dim>::invAT=invA.transpose();
    
    //    template <int dim>
    //    Eigen::Matrix<double,dim,dim> LatticeVector<dim>::cofA=invA*A.determinant();
    
    
    template <int dim>
    const double LatticeBase<dim>::roundTol=FLT_EPSILON;
    
} // end namespace
#endif
