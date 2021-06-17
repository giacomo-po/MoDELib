/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_FEMnodeEvaluation_H_
#define model_FEMnodeEvaluation_H_

#include <Eigen/Dense>
#include <FEMbaseEvaluation.h>

namespace model
{
    
 
    /******************************************************************************/
    template <typename ElementType,int rows,int cols>
    struct FEMnodeEvaluation : public FEMbaseEvaluation<ElementType,rows,cols>
    {
        
        typedef Eigen::Matrix<double,ElementType::dim,1> VectorDim;
        
        const size_t pointID;
        
        /**********************************************************************/
        FEMnodeEvaluation(const size_t& _pointID,const VectorDim& _P) :
        /* init */ FEMbaseEvaluation<ElementType,rows,cols>(_P)
        /* init */,pointID(_pointID)
        {
        }
        
        /**********************************************************************/
        const FEMnodeEvaluation& operator=(const Eigen::Matrix<double,rows,cols>& val)
        {
            static_cast<Eigen::Matrix<double,rows,cols>*>(this)->operator=(val);
            return *this;
        }
        
    };

}
#endif
