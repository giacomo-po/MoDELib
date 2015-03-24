/* This file is part of PIL, the Particle Interaction Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LatticeBase_h_
#define model_LatticeBase_h_

#include <Eigen/Dense>

namespace model
{
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
        
        /**********************************************************************/
        static void setLatticeMatrix(const Eigen::Matrix<double,dim,dim,1>& A_in)
        {
            A=A_in;
            AT=A.transpose();
            invA=A.inverse();
            invAT=invA.transpose();
//            cofA=invA*A.determinant();
            
            std::cout<<"Lattice basis (in columns) =\n"<<A<<std::endl;
            std::cout<<"Lattice reciprocal basis (in columns) =\n"<<invAT<<std::endl;
            
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
    
//    
//    template <int dim>
//    double LatticeVector<dim>::roundTol=FLT_EPSILON;
    
} // end namespace
#endif
