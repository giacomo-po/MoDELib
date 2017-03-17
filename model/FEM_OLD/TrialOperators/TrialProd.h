/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_TrialProd_H_
#define model_TrialProd_H_

#include <utility> // for std::move
#include <model/FEM/TrialOperators/ExpressionRef.h>
#include <model/Utilities/TypeTraits.h>
#include <model/Mesh/Simplex.h>
#include <model/FEM/TrialOperators/TrialBase.h>
#include <model/FEM/TrialOperators/TrialExpressionBase.h>
#include <model/FEM/Constant.h>
#include <model/FEM/TrialOperators/EvalExpression.h>


namespace model
{
    
    template<int rows1,int cols1, int rows2>
    struct TrialProdRows
    {
        static_assert(cols1==rows2,"YOU ARE MULTIPLYING TrialExpressionBaseS WITH DIFFERENT SIZE.");
        constexpr static int rows=rows1;
    };
    
    template<int rows2>
    struct TrialProdRows<1,1,rows2>
    {
        constexpr static int rows=rows2;
    };
    
    template<>
    struct TrialProdRows<1,1,1>
    {
        constexpr static int rows=1;
    };
    
    
    
    /**************************************************************************/
	/**************************************************************************/
    template <typename T1, typename T2>
	class TrialProd : public TrialBase<typename T2::TrialFunctionType>,
    /*             */ public TrialExpressionBase<TrialProd<T1,T2> >
    {
        typedef typename T2::TrialFunctionType _TrialFunctionType;
        typedef TrialBase<_TrialFunctionType> TrialBaseType;
        typedef typename TypeTraits<_TrialFunctionType>::ElementType _ElementType;
        typedef typename TypeTraits<_TrialFunctionType>::BaryType _BaryType;
        constexpr static int dim=TypeTraits<_TrialFunctionType>::dim;
        constexpr static int dofPerElement=TypeTraits<_TrialFunctionType>::dofPerElement;


//        const T1& op1;  // first  operand (the Constant).
//        const T2& op2;
        
        internal::ExpressionRef<T1> op1;
        internal::ExpressionRef<T2> op2;
        
    public:

        constexpr static int rows=TrialProdRows<T1::rows,T1::cols,T2::rows>::rows;
        constexpr static int cols=T2::cols;

        /**********************************************************************/
        TrialProd(const EvalExpression<T1>& x,
                  const TrialExpressionBase<T2>& y) :
        /* base initialization */ TrialBaseType(y.derived().trial()),
        /* init list           */ op1(x.derived()),
        /* init list           */ op2(y.derived())
        {/*!
          */
                        std::cout<<"TrialProd constructor 1"<<std::endl;
        }
        
        /**********************************************************************/
        TrialProd(EvalExpression<T1>&& x,
                  const TrialExpressionBase<T2>& y) :
        /* base initialization */ TrialBaseType(y.derived().trial()),
        /* init list           */ op1(std::move(x.derived())),
        /* init list           */ op2(y.derived())
        {/*!
          */
                        std::cout<<"TrialProd constructor 2"<<std::endl;
        }
        
        /**********************************************************************/
        TrialProd(const EvalExpression<T1>& x,
                  TrialExpressionBase<T2>&& y) :
        /* base initialization */ TrialBaseType(y.derived().trial()),
        /* init list           */ op1(std::move(x.derived())),
        /* init list           */ op2(std::move(y.derived()))
        {/*!
          */
                        std::cout<<"TrialProd constructor 3"<<std::endl;
        }
        
        /**********************************************************************/
        TrialProd(EvalExpression<T1>&& x,
                  TrialExpressionBase<T2>&& y) :
        /* base initialization */ TrialBaseType(y.derived().trial()),
        /* init list           */ op1(std::move(x.derived())),
        /* init list           */ op2(std::move(y.derived()))
        {/*!
          */
            std::cout<<"TrialProd constructor 4"<<std::endl;

        }
        
        
        /**********************************************************************/
        Eigen::Matrix<double,rows,dofPerElement> sfm(const _ElementType& ele,
                                                     const _BaryType& bary) const
        {
//            return op1.c*op2.sfm(ele,bary);
            return op1()(ele,bary)*op2().sfm(ele,bary);

        }
        
//        /**********************************************************************/
//        Eigen::Matrix<double,TF::dim,TypeTraits<TF>::dofPerElement> sfmGrad(const typename TF::ElementType& ele, const BaryType& bary) const
//        {
//              what is the number of rows in this case?
//            return op1.c*op2.sfmGrad(ele,bary);
//        }
        
        /**********************************************************************/
        Eigen::Matrix<double,rows,1> operator()(const _ElementType& ele,
                                                const _BaryType& bary) const
        {/*!@param[in] ele the element
          * @param[in] bary the vector of barycentric coordinates
          * \returns the value of the Derived expression at bary.
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
            const std::pair<bool,const _ElementType*> temp=this->trial().fe.searchWithGuess(P,guess);
//            return (temp.first? sfm(*(temp.second),temp.second->simplex.pos2bary(P))*this->trial().dofs(*(temp.second)) : Eigen::Matrix<double,rows,1>::Zero());
            Eigen::Matrix<double,rows,1> val(Eigen::Matrix<double,rows,1>::Zero());
            if(temp.first)
            {
                val=sfm(*(temp.second),temp.second->simplex.pos2bary(P))*this->trial().dofs(*(temp.second));
            }
            return val;
        }
        
        /**********************************************************************/
        Eigen::Matrix<double,rows,1> operator()(const Eigen::Matrix<double,dim,1>& P) const
        {/*!@param[in] P the position vector
          * @param[in] guess the Simplex where the search starts
          * \returns the value of the Derived expression at P.
          */
            return operator()(P,&(this->mesh().begin()->second));
        }
        
    };
    
    
//    /**************************************************************************/
//    /**************************************************************************/
//    template <typename T1, typename T2>
//    class TrialProd : public TrialBase<typename T2::TrialFunctionType>,
//    /*             */ public TrialExpressionBase<TrialProd<T1,T2> >
//    {
//        typedef typename T2::TrialFunctionType _TrialFunctionType;
//        typedef TrialBase<_TrialFunctionType> TrialBaseType;
//        typedef typename TypeTraits<_TrialFunctionType>::ElementType _ElementType;
//        typedef typename TypeTraits<_TrialFunctionType>::BaryType _BaryType;
//        constexpr static int dim=TypeTraits<_TrialFunctionType>::dim;
//        constexpr static int dofPerElement=TypeTraits<_TrialFunctionType>::dofPerElement;
//        
//        
//        //        const T1& op1;  // first  operand (the Constant).
//        //        const T2& op2;
//        
//        internal::ExpressionRef<T1> op1;
//        internal::ExpressionRef<T2> op2;
//        
//    public:
//        
//        constexpr static int rows=TrialProdRows<T1::rows,T1::cols,T2::rows>::rows;
//        constexpr static int cols=T2::cols;
//        
//        /**********************************************************************/
//        TrialProd(const T1& x,
//                  const T2& y) :
//        /* base initialization */ TrialBaseType(y.derived().trial()),
//        /* init list           */ op1(x),
//        /* init list           */ op2(y)
//        {/*!
//          */
//            std::cout<<"TrialProd constructor 1"<<std::endl;
//        }
//        
//        /**********************************************************************/
//        TrialProd(T1&& x,
//                  const T2& y) :
//        /* base initialization */ TrialBaseType(y.derived().trial()),
//        /* init list           */ op1(std::move(x)),
//        /* init list           */ op2(y)
//        {/*!
//          */
//            std::cout<<"TrialProd constructor 2"<<std::endl;
//        }
//        
//        /**********************************************************************/
//        TrialProd(const T1& x,
//                  T2&& y) :
//        /* base initialization */ TrialBaseType(y.derived().trial()),
//        /* init list           */ op1(x),
//        /* init list           */ op2(std::move(y))
//        {/*!
//          */
//            std::cout<<"TrialProd constructor 3"<<std::endl;
//        }
//        
//        /**********************************************************************/
//        TrialProd(T1&& x,
//                  T2&& y) :
//        /* base initialization */ TrialBaseType(y.derived().trial()),
//        /* init list           */ op1(std::move(x)),
//        /* init list           */ op2(std::move(y))
//        {/*!
//          */
//            std::cout<<"TrialProd constructor 4"<<std::endl;
//            
//        }
//        
//        
//        /**********************************************************************/
//        Eigen::Matrix<double,rows,dofPerElement> sfm(const _ElementType& ele,
//                                                     const _BaryType& bary) const
//        {
//            //            return op1.c*op2.sfm(ele,bary);
//            return op1()(ele,bary)*op2.derived().sfm(ele,bary);
//            
//        }
//        
//        //        /**********************************************************************/
//        //        Eigen::Matrix<double,TF::dim,TypeTraits<TF>::dofPerElement> sfmGrad(const typename TF::ElementType& ele, const BaryType& bary) const
//        //        {
//        //              what is the number of rows in this case?
//        //            return op1.c*op2.sfmGrad(ele,bary);
//        //        }
//        
//        /**********************************************************************/
//        Eigen::Matrix<double,rows,1> operator()(const _ElementType& ele,
//                                                const _BaryType& bary) const
//        {/*!@param[in] ele the element
//          * @param[in] bary the vector of barycentric coordinates
//          * \returns the value of the Derived expression at bary.
//          */
//            return sfm(ele,bary)*this->trial().dofs(ele);
//        }
//        
//        /**********************************************************************/
//        Eigen::Matrix<double,rows,1> operator()(const Eigen::Matrix<double,dim,1>& P,
//                                                const Simplex<dim,dim>* guess) const
//        {/*!@param[in] P the position vector
//          * @param[in] guess the Simplex where the search starts
//          * \returns the value of the Derived expression at P.
//          */
//            const std::pair<bool,const _ElementType*> temp=this->trial().fe.searchWithGuess(P,guess);
//            //            return (temp.first? sfm(*(temp.second),temp.second->simplex.pos2bary(P))*this->trial().dofs(*(temp.second)) : Eigen::Matrix<double,rows,1>::Zero());
//            Eigen::Matrix<double,rows,1> val(Eigen::Matrix<double,rows,1>::Zero());
//            if(temp.first)
//            {
//                val=sfm(*(temp.second),temp.second->simplex.pos2bary(P))*this->trial().dofs(*(temp.second));
//            }
//            return val;
//        }
//        
//        /**********************************************************************/
//        Eigen::Matrix<double,rows,1> operator()(const Eigen::Matrix<double,dim,1>& P) const
//        {/*!@param[in] P the position vector
//          * @param[in] guess the Simplex where the search starts
//          * \returns the value of the Derived expression at P.
//          */
//            return operator()(P,&(this->mesh().begin()->second));
//        }
//        
//    };

    
    /**************************************************************************/
    template <typename DerivedEval,typename DerivedTrial>
    TrialProd<DerivedEval,DerivedTrial> operator*(const EvalExpression<DerivedEval>& evalExp,
                                                     const TrialExpressionBase<DerivedTrial>& trialExp)
    {
        std::cout<<"TrialProd operator* 1d"<<std::endl;
        return TrialProd<DerivedEval,DerivedTrial>(evalExp,trialExp);
    }
    
    template <typename DerivedEval,typename DerivedTrial>
    TrialProd<DerivedEval,DerivedTrial> operator*(EvalExpression<DerivedEval>&& evalExp,
                                                                             const TrialExpressionBase<DerivedTrial>& trialExp)
    {
        std::cout<<"TrialProd operator* 2d"<<std::endl;
        return TrialProd<DerivedEval,DerivedTrial>(std::move(evalExp),trialExp);
    }
    
    template <typename DerivedEval,typename DerivedTrial>
    TrialProd<DerivedEval,DerivedTrial> operator*(const EvalExpression<DerivedEval>& evalExp,
                                                                            TrialExpressionBase<DerivedTrial>&& trialExp)
    {
        std::cout<<"TrialProd operator* 3d"<<std::endl;
        return TrialProd<DerivedEval,DerivedTrial>(evalExp,std::move(trialExp));
    }
    
    template <typename DerivedEval,typename DerivedTrial>
    TrialProd<DerivedEval,DerivedTrial> operator*(EvalExpression<DerivedEval>&& evalExp,
                                                                            TrialExpressionBase<DerivedTrial>&& trialExp)
    {
        std::cout<<"TrialProd operator* 1d"<<std::endl;
        return TrialProd<DerivedEval,DerivedTrial>(std::move(evalExp),std::move(trialExp));
    }
    
//    /**************************************************************************/
//    template <typename T>
//    TrialProd<Constant<double,1,1>,T> operator*(const double& op1,
//                                                const TrialExpressionBase<T>& op2)
//    {
//        std::cout<<"TrialProd operator* 1b"<<std::endl;
//        return operator*(make_constant(op1),op2);
//    }
//    
//    template <typename T>
//    TrialProd<Constant<double,1,1>,T> operator*(double&& op1,
//                                                const TrialExpressionBase<T>& op2)
//    {
//                std::cout<<"TrialProd operator* 2b"<<std::endl;
//        return operator*(make_constant(std::move(op1)),op2);
//    }
//    
//    template <typename T>
//    TrialProd<Constant<double,1,1>,T> operator*(const double& op1,
//                                                TrialExpressionBase<T>&& op2)
//    {
//                std::cout<<"TrialProd operator* 3b"<<std::endl;
//        return operator*(make_constant(op1),std::move(op2));
//    }
//    
//    template <typename T>
//    TrialProd<Constant<double,1,1>,T> operator*(double&& op1,
//                                                TrialExpressionBase<T>&& op2)
//    {
//                std::cout<<"TrialProd operator* 4b"<<std::endl;
//        return operator*(make_constant(std::move(op1)),std::move(op2));
// //       return TrialProd<Constant<double,1,1>,T>(make_constant(std::move(op1)),std::move(op2));
//    }
    
//    /**************************************************************************/
//    template <int rows1, int cols1, typename T2>
//    TrialProd<Constant<Eigen::Matrix<double,rows1,cols1>,rows1,cols1>,T2> operator*(const Eigen::Matrix<double,rows1,cols1>& c,
//                                                                                    const TrialExpressionBase<T2>& op2)
//    {
//        return operator*(make_constant(c),op2);
//    }
//    
//    template <int rows1, int cols1, typename T2>
//    TrialProd<Constant<Eigen::Matrix<double,rows1,cols1>,rows1,cols1>,T2> operator*(Eigen::Matrix<double,rows1,cols1>&& c,
//                                                                                    const TrialExpressionBase<T2>& op2)
//    {
//        return operator*(make_constant(std::move(c)),op2);
//    }
//    
//    template <int rows1, int cols1, typename T2>
//    TrialProd<Constant<Eigen::Matrix<double,rows1,cols1>,rows1,cols1>,T2> operator*(const Eigen::Matrix<double,rows1,cols1>& c,
//                                                                                    TrialExpressionBase<T2>&& op2)
//    {
//        return operator*(make_constant(c),std::move(op2));
//    }
//    
//    template <int rows1, int cols1, typename T2>
//    TrialProd<Constant<Eigen::Matrix<double,rows1,cols1>,rows1,cols1>,T2> operator*(Eigen::Matrix<double,rows1,cols1>&& c,
//                                                                                    TrialExpressionBase<T2>&& op2)
//    {
//        return operator*(make_constant(std::move(c)),std::move(op2));
//    }
    
}	// close namespace
#endif


//    /**************************************************************************/
//    template <typename T1,int rows1, int cols1, typename T2>
//    TrialProd<Constant<T1,rows1,cols1>,T2> operator*(const Constant<T1,rows1,cols1>& c,
//                                                     const TrialExpressionBase<T2>& op2)
//    {
//                std::cout<<"TrialProd operator* 1a"<<std::endl;
//        return TrialProd<Constant<T1,rows1,cols1>,T2>(c,op2);
//    }
//
//    template <typename T1,int rows1, int cols1, typename T2>
//    TrialProd<Constant<T1,rows1,cols1>,T2> operator*(Constant<T1,rows1,cols1>&& c,
//                                                     const TrialExpressionBase<T2>& op2)
//    {
//                std::cout<<"TrialProd operator* 2a"<<std::endl;
//        return TrialProd<Constant<T1,rows1,cols1>,T2>(std::move(c),op2);
//    }
//
//    template <typename T1,int rows1, int cols1, typename T2>
//    TrialProd<Constant<T1,rows1,cols1>,T2> operator*(const Constant<T1,rows1,cols1>& c,
//                                                     TrialExpressionBase<T2>&& op2)
//    {
//                std::cout<<"TrialProd operator* 3a"<<std::endl;
//        return TrialProd<Constant<T1,rows1,cols1>,T2>(c,std::move(op2));
//    }
//
//    template <typename T1,int rows1, int cols1, typename T2>
//    TrialProd<Constant<T1,rows1,cols1>,T2> operator*(Constant<T1,rows1,cols1>&& c,
//                                                     TrialExpressionBase<T2>&& op2)
//    {
//                std::cout<<"TrialProd operator* 4a"<<std::endl;
//        return TrialProd<Constant<T1,rows1,cols1>,T2>(std::move(c),std::move(op2));
//
//    }
