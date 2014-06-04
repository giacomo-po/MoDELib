/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


/**********************************************************************/
template <typename BaryType>
Eigen::Matrix<double,Eigen::Dynamic,1> operator()(const ElementType& ele, const BaryType& bary) const
{/*!@param[in] ele the element
  * @param[in] bary the vector of barycentric coordinates
  * \returns the value of the Derived expression at bary.
  */
    return sfm(ele,bary)*trial().dofs(ele);
}

//        /**********************************************************************/
//        template <int dim>
//        Eigen::Matrix<double,Eigen::Dynamic,1> operator()(const Eigen::Matrix<double,dim,1>& P, const Simplex<dim,dim>* guess) const
//        {/*!@param[in] ele the element
//          * @param[in] bary the vector of barycentric coordinates
//          * \returns the value of the Derived expression at bary.
//          */
//            const std::pair<bool,const ElementType*> temp=derived().fe.searchWithGuess(P,guess);
//
//
//            Eigen::Matrix<double,Eigen::Dynamic,1> temp1(derived().sfm(*temp.second,bary)*derived().trial().dofs(*temp.second));
//
//            return temp.first?  : ;
//        }





}	// close namespace
#endif