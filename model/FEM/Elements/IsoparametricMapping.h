/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2014 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_IsoparametricMapping_H_
#define model_IsoparametricMapping_H_

#include <Eigen/Dense>

namespace model
{
    
    /**************************************************************************/
	/**************************************************************************/
    /*!\brief Class template that defines the isoparametric coordinate
     * transformation between actual and standared domain.
     */
    template<typename ElementType>
	class IsoparametricMapping
    {
        
        constexpr static int dim=ElementType::dim;
        
        const ElementType& ele;
        
    public:
        
        /**********************************************************************/
        IsoparametricMapping(const ElementType& ele_in) : ele(ele_in)
        {
            
        }
        
        /**********************************************************************/
        Eigen::Matrix<double,dim,1> position(const Eigen::Matrix<double,1,dim+1>& bary) const
        {/*!@param[in] bary the vector of baricentric coordinates
          * \returns the position in the actual domain corresponding to bary.
          */
            return ele.Xe*ele.sf(bary).transpose();
        }
        
        /**********************************************************************/
        Eigen::Matrix<double,dim,dim> Fs(const Eigen::Matrix<double,1,dim+1>& bary) const
        {/*!@param bary the vector of baricentric coordinates
          * \returns the Jacobian matrix dx_i / ds_j, evaulated at bary.
          */
            return ele.Xe*ele.gradS(bary).transpose();
        }
        
        /**********************************************************************/
        double J(const Eigen::Matrix<double,1,dim+1>& bary) const
        {/*!@param bary the vector of baricentric coordinates
          * \returns the absolute value of the determinant of the Jacobian
          * matrix dx_i / ds_j, evaulated at bary.
          */
            return Fs(bary).determinant();
        }
        
        /**********************************************************************/
        double absJ(const Eigen::Matrix<double,1,dim+1>& bary) const
        {/*!@param bary the vector of baricentric coordinates
          * \returns the absolute value of the determinant of the Jacobian
          * matrix dx_i / ds_j, evaulated at bary.
          */
            return std::fabs(Fs(bary).determinant());
        }
        
        /**********************************************************************/
        Eigen::Matrix<double,dim,dim> Gs(const Eigen::Matrix<double,1,dim+1>& bary) const
        {/*!@param bary the vector of baricentric coordinates
          * \returns the inverse Jacobian matrix ds_i / dx_j, evaulated at bary.
          */
            return Fs(bary).inverse();
        }
        
        /**********************************************************************/
        Eigen::Matrix<double,dim,1> jGN(const Eigen::Matrix<double,1,dim+1>& bary, const int& n) const
        {
            assert(n>=0);
            assert(n<dim+1);
            assert(bary(n)==0.0);
            
            return Gs(bary).transpose()*BarycentricTraits<dim>::NdA.col(n)*absJ(bary);
        }
        
    };
    
}	// close namespace
#endif

//        /**********************************************************************/
//        Eigen::Matrix<double,dim,1> outNormal(const Eigen::Matrix<double,1,dim+1>& bary, const int& n) const
//        {
//            assert(n>=0);
//            assert(n<dim+1);
//            assert(bary(n)==0.0);
//
//            /* a matrix having in colums the derivatives of the position vector
//             * with respect to barycentric coordinates */
//            const Eigen::Matrix<double,dim,dim+1> dxdL(Xe*diff(bary).transpose());
//
//            // remove the n-th derivative
//            Eigen::Matrix<double,dim,dim> dxdLt;
//            int c(0);
//            for (int k=0;k<dim+1;++k)
//            {
//                if (k!=n)
//                {
//                    dxdLt.col(c)=dxdL.col(k);
//                    c++;
//                }
//            }
//
//            // compute dim-1 tangent vectors
//            Eigen::Matrix<double,dim,dim-1> tan;
//            for (int k=0;k<dim-1;++k)
//            {
//                tan.col(k)=dxdLt.col(k)-dxdLt.col(k+1);
//            }
//
//            //std::cout<<"dxdL=\n"<<dxdL<<std::endl;
//
//            //            Eigen::Matrix<double,dim,1> temp(-dxdL.col(n));
//            Eigen::Matrix<double,dim,1> temp(position(bary)-position(Eigen::Matrix<double,dim+1,1>::Constant(1.0/(dim+1))));
//
//            //const double jF(Fs(bary).determinant());
//            //temp*=(jF/std::fabs(jF));
//            for (int k=0;k<dim-1;++k)
//            {
//                const double n2(tan.col(k).squaredNorm());
//                if(n2>0.0) // this should not happen
//                {
//                    //temp += dxdL.col(n).dot(dxdL.col(k))*dxdL.col(k)/dxdL.col(k).squaredNorm();
//                    temp -= (temp.dot(tan.col(k))*tan.col(k)/n2).eval();
//                }
//            }
//            return temp.normalized();
//        }

