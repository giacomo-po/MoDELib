/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2012 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _model_SingleFieldPoint_h
#define _model_SingleFieldPoint_h

#include <FieldPoint.h>


namespace model
{
    
    /******************************************************************************/
    template<typename Field>
    struct SingleFieldPoint :
    /* inheritance */ public FieldPoint<SingleFieldPoint<Field>,Field::dim,Field>
    {
        const Eigen::Matrix<double,Field::dim,1> P;
        
        SingleFieldPoint(const Eigen::Matrix<double,Field::dim,1>& Pin, const bool& enabled) :
        /* init */ FieldPoint<SingleFieldPoint<Field>,Field::dim,Field>(enabled),
        /* init */ P(Pin)
        {
        
        }
        
        
        typename Field::MatrixType field() const
        {
            return FieldPoint<SingleFieldPoint<Field>,Field::dim,Field>::template field<Field>();
        }
    };
    
} // end namespace
#endif
