/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2014 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef  model_BarycentricTraits_H_
#define  model_BarycentricTraits_H_

#include <set>
#include <Eigen/Dense>

namespace model {
	
    
    /**************************************************************************/
	/**************************************************************************/
    /*!Class template that performs auxiliary barycentric-to-standard
     * coordinate transformation in a Simplex.
     */
	template<short int dim>
	struct BarycentricTraits
    {
        static_assert(dim>0,"dim must be > 0");

        
        typedef Eigen::Matrix<double,dim,dim> MatrixDimD;
        typedef Eigen::Matrix<double,dim+1,dim+1> MatrixHigherDimD;
        typedef Eigen::Matrix<double,dim+1,dim> MatrixHigherDimDimD;
        typedef Eigen::Matrix<double,dim,dim+1> MatrixDimHigherDimD;

        typedef Eigen::Matrix<double,dim  ,1> VectorDimD;
        typedef Eigen::Matrix<double,dim+1,1> VectorHigherDimD;

        static const MatrixHigherDimD L2X;
        static const MatrixHigherDimD X2L;
        static const MatrixHigherDimDimD dLdX;
        static const MatrixDimHigherDimD NdA;
        
        /**********************************************************************/
        static VectorHigherDimD x2l(const VectorDimD& X);
        
        /**********************************************************************/
        static VectorHigherDimD x2l(const double& X);
        
        /**********************************************************************/
        static VectorDimD l2x(const VectorHigherDimD& L);

        /**********************************************************************/
        static Eigen::Matrix<double,dim+1,1> face2domainBary(const Eigen::Matrix<double,dim,1>& b1,const int& boundaryFace);
        
    private:
        
        /**********************************************************************/
        static MatrixHigherDimD get_L2X();
        
        /**********************************************************************/
        static MatrixHigherDimD get_X2L();
        
        /**********************************************************************/
        static MatrixDimHigherDimD get_NdA();
        
	};
    
//    // Declare static data member
//    template<short int dim>
//    const typename BarycentricTraits<dim>::MatrixHigherDimD BarycentricTraits<dim>::L2X=BarycentricTraits<dim>::get_L2X();
//    
//	template<short int dim>
//    const typename BarycentricTraits<dim>::MatrixHigherDimD BarycentricTraits<dim>::X2L=BarycentricTraits<dim>::get_X2L();
//
//    template<short int dim>
//    const typename BarycentricTraits<dim>::MatrixHigherDimDimD BarycentricTraits<dim>::dLdX=BarycentricTraits<dim>::X2L. template block<dim+1,dim>(0,0);
//    
//    template<short int dim>
//    const typename BarycentricTraits<dim>::MatrixDimHigherDimD BarycentricTraits<dim>::NdA=BarycentricTraits<dim>::get_NdA();

	
}
#endif
