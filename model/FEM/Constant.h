/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_Constatnt_H_
#define model_Constatnt_H_

#include <utility> // for std::move

#include <Eigen/Dense>
//#include <EvalExpression.h>
#include <EvalFunction.h>


namespace model
{

    /**************************************************************************/
	/**************************************************************************/
    template <typename T, int _rows, int _cols>
	struct Constant : public EvalFunction<Constant<T,_rows,_cols> >
    {

        constexpr static int rows=_rows;
        constexpr static int cols=_cols;

        
        const T c;
//        internal::ExpressionRef<T> c;

        
        /**********************************************************************/
        Constant(const T& c_in) : c(c_in)
        {
//            std::cout<<"Constant Constructor 1"<<std::endl;
//            std::cout<<"c="<<c<<std::endl;
        }
        
//        /**********************************************************************/
//        Constant(T&& c_in) : c(std::move(c_in))
//        {
//            std::cout<<"Constant Constructor 2"<<std::endl;
//        }
        
        /**********************************************************************/
        template<typename ElementType, typename BaryType>
        const T& operator() (const ElementType&, const BaryType&) const
        {/*!@param[in] elem the element
          * @param[in] bary the barycentric cooridinate
          *\returns the constant c.
          */

            return c;
        }
        
//        /**********************************************************************/
//        EvalExpression<Constant<T,_rows,_cols>> eval() const
//        {
//            return EvalExpression<Constant<T,_rows,_cols>>(*this);
//        }
        
    };
    
    // Operators
    template <int rows, int cols>
    Constant<Eigen::Matrix<double,rows,cols>,rows,cols> make_constant(const Eigen::Matrix<double,rows,cols>& c)
    {
        return Constant<Eigen::Matrix<double,rows,cols>,rows,cols>(c);
    }
    
//    template <int rows, int cols>
//    Constant<Eigen::Matrix<double,rows,cols>,rows,cols> make_constant(Eigen::Matrix<double,rows,cols>&& c)
//    {
//        return Constant<Eigen::Matrix<double,rows,cols>,rows,cols>(std::move(c));
//    }

    static inline Constant<double,1,1> make_constant(const double& c)
    {
        return Constant<double,1,1>(c);
    }
    
//    Constant<double,1,1> make_constant(double&& c)
//    {
//        return Constant<double,1,1>(std::move(c));
//    }
    
//    /**************************************************************************/
//    template <typename T,int rows,int cols>
//    EvalExpression<Constant<T,rows,cols>> eval(const Constant<T,rows,cols>& con)
//    {
//        return EvalExpression<Constant<T,rows,cols>>(con);
//    }
//
//    template <typename T,int rows,int cols>
//    EvalExpression<Constant<T,rows,cols>> eval(Constant<T,rows,cols>&& con)
//    {
//        return EvalExpression<Constant<T,rows,cols>>(con);
//    }
    
//    Constant<double,1,1> operator()(const double& c)
//    {
//        return Constant<double,1,1>(c);
//
//    }
    
    
}	// close namespace
#endif


//    template <typename Derived>
//    int make_constant(const Eigen::DenseBase<Derived>& c)
//    {
//        return 1;
//    }

//    template <typename Derived>
//    Constant<Derived,Derived::rows,2> make_constant(const Eigen::DenseBase<Derived>& c)
//    {
//        return Constant<Derived,Derived::rows,2>(c);
//    }
