/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_BilinearWeakSum_H_
#define model_BilinearWeakSum_H_

#include <iterator>
#include <model/Utilities/AreSameType.h>
#include <model/FEM/WeakForms/BilinearWeakExpression.h>

namespace model
{
 	/**************************************************************************/
	/**************************************************************************/
	/*! \brief A class template that performs the sum of two BilinearWeakExpression.
	 */
    template <typename T1,typename T2>
	struct  BilinearWeakSum : public BilinearWeakExpression<BilinearWeakSum<T1,T2> >
    {
        
        static_assert(AreSameType<typename T1::TrialFunctionType,typename T2::TrialFunctionType>::value,"YOU ARE SUMMING BilinearWeakForms OF DIFFERENT TRIALFUNCTIONS.");
        typedef typename T1::TrialFunctionType TrialFunctionType;
        
        //! first operand
        const T1& op1; // expression could be a temporary, so copy by value
        //! second operand
        const T2& op2; // expression could be a temporary, so copy by value
        
        /**********************************************************************/
        BilinearWeakSum(const BilinearWeakExpression<T1>& x, const BilinearWeakExpression<T2>& y) :
        /* init list           */ op1(x.derived()),
        /* init list           */ op2(y.derived())
        {/*!
          */
            
        }
        
        /**********************************************************************/
        std::vector<Eigen::Triplet<double> >  globalTriplets() const
        {
            std::vector<Eigen::Triplet<double> > gt1=op1.globalTriplets();
            std::vector<Eigen::Triplet<double> > gt2=op2.globalTriplets();
            gt1.insert(
                        gt1.end(),
                        std::make_move_iterator(gt2.begin()),
                        std::make_move_iterator(gt2.end())
                        );
            return gt1;
        }
        
    };
    
    
    /**************************************************************************/
    template <typename T1, typename T2>
    BilinearWeakSum<T1,T2> operator+(const BilinearWeakExpression<T1>& op1,
                                     const BilinearWeakExpression<T2>& op2)
    {
        return  BilinearWeakSum<T1,T2>(op1,op2);
    }
    
    
}	// close namespace
#endif
