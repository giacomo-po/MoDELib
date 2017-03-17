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
