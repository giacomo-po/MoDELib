/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2018 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2018 by Yinan Cui <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_StressStraight_H_
#define model_StressStraight_H_

#ifndef _MODEL_NON_SINGULAR_DD_
#define _MODEL_NON_SINGULAR_DD_ 1
#endif

#include <Eigen/Dense>
#include <PolycrystallineMaterialBase.h>
#include <DislocationFieldBase.h>

namespace model
{

    template <int dim,typename Scalar=double>
    class StressStraight
    {
        typedef Eigen::Matrix<Scalar,dim,dim> MatrixDim;
        typedef Eigen::Matrix<Scalar,dim,1>   VectorDim;
        MatrixDim nonSymmStress_kernel(const VectorDim& r) const;
        VectorDim displacement_kernel(const VectorDim& r) const;
        
    public:
        
        const PolycrystallineMaterialBase& material;
        const VectorDim P0;
        const VectorDim P1;
        const VectorDim b;
        const double length;
        const VectorDim t;
        const VectorDim bCt;
        
        StressStraight(const PolycrystallineMaterialBase& material_in,const VectorDim& _P0,const VectorDim& _P1, const VectorDim& _b);
        MatrixDim nonSymmStress(const VectorDim& x) const;
        MatrixDim stress(const VectorDim& x) const;
        VectorDim displacement(const VectorDim& x) const;
	};	
	
}
#endif
