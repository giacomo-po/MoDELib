/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_Constatnt_H_
#define model_Constatnt_H_

#include <Eigen/Dense>
#include <model/FEM/TrialOperators/EvalExpression.h>


namespace model
{
    
 
    /**************************************************************************/
	/**************************************************************************/
    template <typename T, int _rows, int _cols>
	struct Constant : public EvalExpression<Constant<T,_rows,_cols> >
    {
//        enum {rows=_rows};
//        enum {cols=_cols};
        constexpr static int rows=_rows;
        constexpr static int cols=_cols;

        
        const T c;
        
        /**********************************************************************/
        Constant(const T& c_in) : c(c_in)
        {
        
        }
        
        /**********************************************************************/
        template<typename ElementType, typename BaryType>
        const T& operator() (const ElementType& ele, const BaryType& bary) const
        {/*!@param[in] elem the element
          * @param[in] bary the barycentric cooridinate
          *\returns the constant c.
          */
            return c;
        }
        
    };
    


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
    
    template <int rows, int cols>
    Constant<Eigen::Matrix<double,rows,cols>,rows,cols> make_constant(const Eigen::Matrix<double,rows,cols>& c)
    {
        return Constant<Eigen::Matrix<double,rows,cols>,rows,cols>(c);
    }

    Constant<double,1,1> make_constant(const double& c)
    {
        return Constant<double,1,1>(c);
    }
    
    
    
}	// close namespace
#endif