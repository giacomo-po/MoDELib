/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2012 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _model_DisplacementPoint_h
#define _model_DisplacementPoint_h

#include <assert.h>
#include <FieldPoint.h>


namespace model
{
    
    /******************************************************************************/
    template<int dim>
    struct DisplacementPoint : public Eigen::Matrix<double,dim,1>
    {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        
        const size_t pointID;
        const VectorDim P;
        
        /**********************************************************************/
        DisplacementPoint(const size_t& _pointID,const VectorDim& _P) :
        /* init */ Eigen::Matrix<double,dim,1>(VectorDim::Zero()),
        /* init */ pointID(_pointID),
        /* init */ P(_P)
        {
        }
        
        const DisplacementPoint& operator=(const Eigen::Matrix<double,dim,1>& val)
        {
            static_cast<Eigen::Matrix<double,dim,1>*>(this)->operator=(val);
            return *this;
        }
        
    };
    
}
#endif
