/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_TrialProd_H_
#define model_TrialProd_H_

#include <model/FEM/TrialOperators/TrialBase.h>
#include <model/FEM/TrialOperators/TrialExpressionBase.h>
#include <model/FEM/Constant.h>

namespace model
{
    
    template<int rows1,int cols1, int rows2>
    struct TrialProdRows
    {
        static_assert(cols1==rows2,"YOU ARE MULTIPLYING TrialExpressionBaseS WITH DIFFERENT SIZE.");
//        enum {rows=rows1};
        constexpr static int rows=rows1;
    };

//    template<int rows1,int cols1>
//    struct TrialProdRows<rows1,cols1,cols1>
//    {
//        enum {rows=rows1};
//    };
    
    template<int rows2>
    struct TrialProdRows<1,1,rows2>
    {
//        enum {rows=rows2};
        constexpr static int rows=rows2;

    };
    
    template<>
    struct TrialProdRows<1,1,1>
    {
        constexpr static int rows=1;
//        enum {rows=1};
    };
    
    
    
    /**************************************************************************/
	/**************************************************************************/
    template <typename T1, typename T2>
	class TrialProd : public TrialBase<typename T2::TrialFunctionType>,
    /*             */ public TrialExpressionBase<TrialProd<T1,T2> >
    {
        typedef typename T2::TrialFunctionType _TrialFunctionType;
        typedef TrialBase<_TrialFunctionType> TrialBaseType;


        const T1 op1;  // first  operand (the Constant). This is copied by value because the input Constant coudl be a temporary.
        const T2 op2;  // Input expression could be a temporary, so copy by value
        
    public:

        constexpr static int rows=TrialProdRows<T1::rows,T1::cols,T2::rows>::rows;

        /**********************************************************************/
        TrialProd(const T1& x, const TrialExpressionBase<T2>& y) :
        /* base initialization */ TrialBaseType(y.derived().trial()),
        /* init list           */ op1(x),
        /* init list           */ op2(y.derived())
        {/*!
          */
            
        }
        
        /**********************************************************************/
        Eigen::Matrix<double,rows,TypeTraits<_TrialFunctionType>::dofPerElement> sfm(const typename TypeTraits<_TrialFunctionType>::ElementType& ele,
                                                                                     const typename TypeTraits<_TrialFunctionType>::BaryType& bary) const
        {
            return op1.c*op2.sfm(ele,bary);
        }
        
//        /**********************************************************************/
//        Eigen::Matrix<double,TF::dim,TypeTraits<TF>::dofPerElement> sfmGrad(const typename TF::ElementType& ele, const BaryType& bary) const
//        {
//              what is the number of rows in this case?
//            return op1.c*op2.sfmGrad(ele,bary);
//        }
        
    };
    
    
    /**************************************************************************/
    template <typename T1,int rows1, int cols1, typename T2>
    TrialProd<Constant<T1,rows1,cols1>,T2> operator*(const Constant<T1,rows1,cols1>& c, const TrialExpressionBase<T2>& op2)
    {
        return TrialProd<Constant<T1,rows1,cols1>,T2>(c,op2);
    }
    
    /**************************************************************************/
    template <typename T>
    TrialProd<Constant<double,1,1>,T> operator*(const double& op1, const TrialExpressionBase<T>& op2)
    {
        return operator*(make_constant(op1),op2);
    }
    
    /**************************************************************************/
    template <int rows1, int cols1, typename T2>
    TrialProd<Constant<Eigen::Matrix<double,rows1,cols1>,rows1,cols1>,T2> operator*(const Eigen::Matrix<double,rows1,cols1>& c, const TrialExpressionBase<T2>& op2)
    {
        return operator*(make_constant(c),op2);
    }

}	// close namespace
#endif



//    /**************************************************************************/
//    template <typename Derived, typename T2, typename TF, int rows2>
//    TrialProd<Constant<Derived,Derived::rowsAtCompileTime,Derived::colsAtCompileTime>,T2> operator*(const Eigen::DenseBase<Derived>& op1, const TrialExpressionBase<T2,TF,rows2>& op2)
//    {
//        return operator*(make_constant(op1),op2);
//    }
//
//    /**************************************************************************/
//    template <typename Derived, int N, typename FE>
//    TrialProd<Constant<Derived,Derived::rowsAtCompileTime,Derived::colsAtCompileTime>,TrialExpression<TrialFunction<N,FE>,N> > operator*(const Eigen::DenseBase<Derived>& op1, const TrialFunction<N,FE>& op2)
//    {
//        return operator*(make_constant(op1),TrialExpression<TrialFunction<N,FE>,N>(op2));
//    }
//
//    /**************************************************************************/
//    template <int rows1, int cols1, int N, typename FE>
//    TrialProd<Constant<Eigen::Matrix<double,rows1,cols1>,rows1,cols1>,TrialExpression<TrialFunction<N,FE> > > operator*(const Eigen::Matrix<double,rows1,cols1>& op1, const TrialFunction<N,FE>& op2)
//    {
//        return operator*(make_constant(op1),TrialExpression<TrialFunction<N,FE> >(op2));
//    }
//
//
//
//    /**************************************************************************/
//    template <int N, typename FE>
//    TrialProd<Constant<double,1,1>,TrialExpression<TrialFunction<N,FE> > > operator*(const double& op1, const TrialFunction<N,FE>& op2)
//    {
//        return operator*(make_constant(op1),TrialExpression<TrialFunction<N,FE> >(op2));
//    }