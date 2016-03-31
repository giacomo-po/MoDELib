/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_TrialBase_H_
#define model_TrialBase_H_

//#include <model/FEM/TrialOperators/TestExpression.h>

namespace model
{
    	
    /**************************************************************************/
	/**************************************************************************/
    /*!\brief A class template that provides the base for all the expressions
     * involving a TrialFunction. The expression-template mechanism is based 
     * on the CRTP pattern.
     */
	template<typename _TrialFunctionType>
	class TrialBase
    {

    public:
        typedef typename TypeTraits<_TrialFunctionType>::FiniteElementType FiniteElementType;
        typedef typename TypeTraits<_TrialFunctionType>::ElementType ElementType;
        typedef typename TypeTraits<_TrialFunctionType>::NodeType NodeType;
        typedef _TrialFunctionType TrialFunctionType;

        //        const FiniteElementType& fe;
    private:
        
        const _TrialFunctionType& tf;

        
    public:
        
//        const TrialFunctionType& tf;

        
        /**********************************************************************/
        TrialBase(const TrialFunctionType& tf_in) : tf(tf_in)
        {/*!\param[in] tf_in a const reference to a TrialFunction
          * Constructor initializes tf.
          */
            
        }
        
        /**********************************************************************/
        const ElementType& element(const size_t& n) const
        {/*!@param[in] n the node ID
          * \returns a const reference to the n-th node in the FiniteElement
          */
            return tf.fe.element(n);
        }

        /**********************************************************************/
        const TrialFunctionType& trial() const
        {
            return tf;
        }
        
        /**********************************************************************/
        size_t elementSize() const
        {/*!\returns the number of elements in the FiniteElement
          */
            return tf.fe.elementSize();
        }
        
        /**********************************************************************/
        typename FiniteElementType::ElementContainerType::const_iterator elementBegin() const
        {
            return tf.fe.elementBegin();
        }
        
        /**********************************************************************/
        typename FiniteElementType::ElementContainerType::const_iterator elementEnd() const
        {
            return tf.fe.elementEnd();
        }
        
        /**********************************************************************/
        size_t nodeSize() const
        {
            return tf.fe.nodeSize();
        }

        /**********************************************************************/
        const NodeType& node(const size_t& n) const
        {
            return tf.fe.node(n);
        }

        
    };
    
    
    
}	// close namespace
#endif

//        /**********************************************************************/
//        void operator()(const Eigen::Matrix<double,TypeTraits<TrialFunctionType>::dim,1>& p) const
//        {
//            std::cout<<"Finish implementation here"<<std::endl;
//        }


//        template <typename ElementType>
//        Eigen::Matrix<double,rows,derived::nodesPerElement> at(const typename Derived::ElementType& ele, const Eigen::Matrix<double,Derived::dim+1,1>& bary) const
//        {
//
//            return derived().template sfm(ele,bary);
//
//        }