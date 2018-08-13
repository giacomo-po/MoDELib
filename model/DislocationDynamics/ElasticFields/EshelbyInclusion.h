/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2018 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2018 by Yinan Cui <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_EshelbyInclusion_H_
#define model_EshelbyInclusion_H_

#ifndef _MODEL_NON_SINGULAR_DD_
#define _MODEL_NON_SINGULAR_DD_ 1
#endif

#include <Eigen/Dense>
#include <model/DislocationDynamics/Materials/Material.h>
//#include <model/DislocationDynamics/NearestNeighbor/DislocationStress.h>

namespace model
{
    
    
    template <int dim,typename Scalar=double>
    class EshelbyInclusion
    {
        typedef Eigen::Matrix<Scalar,dim,dim> MatrixDim;
        typedef Eigen::Matrix<Scalar,dim,1>   VectorDim;
        

        VectorDim P;
        double R;

    public:
        
        EshelbyInclusion(const VectorDim& _P,const double& _R) :
        /* init */ P(_P)
        /* init */,R(_R)
        {
        
        }
        
        
	};	
	
}
#endif

