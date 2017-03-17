/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_TrialGrad_H_
#define model_TrialGrad_H_

#include <utility> // for std::move
#include <model/FEM/TrialOperators/ExpressionRef.h>
#include <type_traits> //static_assert
#include <model/Utilities/TypeTraits.h>
#include <model/Mesh/Simplex.h>
#include <model/FEM/TrialOperators/TrialBase.h>
#include <model/FEM/TrialOperators/TrialExpressionBase.h>

namespace model
{
    
    /**************************************************************************/
	/**************************************************************************/
    template <typename T>
    class TrialGrad : public TrialBase<typename T::TrialFunctionType>,
    /*             */ public TrialExpressionBase<TrialGrad<T> >
    {
        
        
        
        //const T& op;
        internal::ExpressionRef<T> op;

    public:
        
        typedef typename T::TrialFunctionType _TrialFunctionType;
        typedef TrialBase<_TrialFunctionType> TrialBaseType;
        typedef typename TypeTraits<_TrialFunctionType>::ElementType ElementType;
        typedef typename TypeTraits<_TrialFunctionType>::BaryType BaryType;
        constexpr static int dim=TypeTraits<_TrialFunctionType>::dim;
        constexpr static int dofPerElement=TypeTraits<_TrialFunctionType>::dofPerElement;
        constexpr static int rows=_TrialFunctionType::dim*T::rows;
        constexpr static int cols=T::cols;
        typedef Eigen::Matrix<double,rows,dofPerElement> ShapeFunctionMatrixType;

        /**********************************************************************/
        TrialGrad(const TrialExpressionBase<T>& x) :
        /* base init */ TrialBaseType(x.derived().trial()),
        /* init list */ op(x.derived())
        {
            std::cout<<"TrialGrad constructor 1"<<std::endl;
        }
        
        /**********************************************************************/
        TrialGrad(TrialExpressionBase<T>&& x) :
        /* base init */ TrialBaseType(x.derived().trial()),
        /* init list */ op(std::move(x.derived()))
        {
            std::cout<<"TrialGrad constructor 2"<<std::endl;
        }
        
        /**********************************************************************/
        ShapeFunctionMatrixType sfm(const ElementType& ele,
                                           const BaryType& bary) const
        {
            return op().sfmGrad(ele,bary);
        }
        
        //        /**********************************************************************/
        //        void sfmGrad(const typename TF::ElementType& ele, const Eigen::Matrix<double,TF::dim+1,1>& bary) const
        //        {
        //            static_assert(0,"SECOND GRADIENTS ARE NOT SUPPORTED YET.");
        //        }
        
        /**********************************************************************/
        Eigen::Matrix<double,rows,1> operator()(const ElementType& ele,
                                                const BaryType& bary) const
        {/*!@param[in] ele the element
          * @param[in] bary the vector of barycentric coordinates
          * \returns the value of the Derived expression at bary.
          *
          * \todo: in order to be optimized, this function should be Derived-specific
          */
            return sfm(ele,bary)*this->trial().dofs(ele);
        }
        
        /**********************************************************************/
        Eigen::Matrix<double,rows,1> operator()(const Eigen::Matrix<double,dim,1>& P,
                                                const Simplex<dim,dim>* guess) const
        {/*!@param[in] P the position vector
          * @param[in] guess the Simplex where the search starts
          * \returns the value of the Derived expression at P.
          */
            const std::pair<bool,const ElementType*> temp=this->trial().fe.searchWithGuess(P,guess);
            //            return (temp.first? sfm(*(temp.second),temp.second->simplex.pos2bary(P))*this->trial().dofs(*(temp.second)) : Eigen::Matrix<double,rows,1>::Zero());
            Eigen::Matrix<double,rows,1> val(Eigen::Matrix<double,rows,1>::Zero());
            if(temp.first)
            {
                val=sfm(*(temp.second),temp.second->simplex.pos2bary(P))*this->trial().dofs(*(temp.second));
            }
            return val;
        }

    };
    
    /**************************************************************************/
    template <typename T>
    TrialGrad<T> grad(const TrialExpressionBase<T>& x)
    {
        return TrialGrad<T>(x);
    }
    
    /**************************************************************************/
    template <typename T>
    TrialGrad<T> grad(TrialExpressionBase<T>&& x)
    {
        return TrialGrad<T>(std::move(x));
    }
    
}
#endif


//        /**********************************************************************/
//        Eigen::Matrix<double,rows,1> operator()(const Eigen::Matrix<double,dim,1>& P) const
//        {/*!@param[in] P the position vector
//          * @param[in] guess the Simplex where the search starts
//          * \returns the value of the Derived expression at P.
//          */
//            return operator()(P,&(op.trial.fe.mesh.begin()->second));
//        }
