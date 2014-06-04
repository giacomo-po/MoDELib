/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationNegativeStress_H_
#define model_DislocationNegativeStress_H_



namespace model
{
	template <typename DislocationNetworkType>
    struct DislocationNegativeStress :
    /* inheritance */ public EvalExpression<DislocationNegativeStress<DislocationNetworkType> >
    {
        constexpr static int dim=DislocationNetworkType::dim;
        constexpr static int rows=dim;
        constexpr static int cols=dim;
        
        const DislocationNetworkType& DN;
        
        /**********************************************************************/
        DislocationNegativeStress(const DislocationNetworkType& DN_in) : DN(DN_in)
        {
            
        }
        
        /**********************************************************************/
        template<typename ElementType, typename BaryType>
        Eigen::Matrix<double,dim,dim> operator() (const ElementType& ele, const BaryType& bary) const
        {/*!@param[in] elem the element
          * @param[in] bary the barycentric cooridinate
          *\returns the stress
          */
            
            //Eigen::Matrix<double,dim,1> outNormal;
            
            ///* inheritance */ public FieldPoint<DislocationNegativeStress,dim,DislocationStress<dim> >,
            // Let the DislocationNetwork compute the stress at the field points
            
            //    DN.computeField<SimpleFieldPoint,StressField>(fieldPoints);
            
            // Ouput results to file
            
            return -DN.stress(ele.position(bary));
        }
        
    };

} // namespace model
#endif


