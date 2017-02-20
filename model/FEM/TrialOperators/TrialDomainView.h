/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_TrialDomainView_H_
#define model_TrialDomainView_H_

//#include <model/Utilities/TypeTraits.h>
//#include <model/FEM/TrialOperators/TrialExpressionBase.h>


namespace model
{
    
    template<typename Derived>
    struct TrialExpressionBase;
    
    /**************************************************************************/
	/**************************************************************************/
    /*!\brief A class template that provides the base for all the expressions
     * involving a TrialFunction. The expression-template mechanism is based
     * on the CRTP pattern.
     */
	template<typename Derived>
	struct TrialDomainView
    {
        
        const TrialExpressionBase<Derived>& trialBase;
        
        
        TrialDomainView(const TrialExpressionBase<Derived>& te) :
        trialBase(te)
        {
        
        }
        

    };
    
    
    /**************************************************************************/
    template <typename T, typename Derived>
    T& operator << (T& os, const TrialDomainView<Derived>& dv)
    {/*!@param[in] os the stream object
      * @param[in] expr the expression
      *
      * Outputs the value of Derived on the faces of the mesh
      */
        constexpr int dim=Derived::TrialFunctionType::dim;
        const Eigen::Matrix<double,dim+1,dim+1> vertexBary(Eigen::Matrix<double,dim+1,dim+1>::Identity());
        
        
        for (typename Derived::FiniteElementType::ElementContainerType::const_iterator eIter =dv.trialBase.derived().trial().fe.elementBegin();
             /*                                                                     */ eIter!=dv.trialBase.derived().trial().fe.elementEnd();
             /*                                                                   */ ++eIter)
        {
//            if(eIter->second.isBoundaryElement())
//            {
//                const std::vector<int> boundaryFaces=eIter->second.boundaryFaces();
//                for (size_t f=0;f<boundaryFaces.size();++f)
//                {
                    for (int v=0;v<dim+1;++v)
                    {
//                        if (v!=boundaryFaces[f])
//                        {
                            os<<eIter->second.position(vertexBary.col(v)).transpose()<<" "
                            <<dv.trialBase.derived()(eIter->second,vertexBary.col(v)).transpose()<<"\n";
//                        }
//                    }
                }
//            }
        }
        return os;
    }
    
    
}	// close namespace
#endif
