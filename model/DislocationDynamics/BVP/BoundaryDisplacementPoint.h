/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2012 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _model_BoundaryDisplacementPoint_h
#define _model_BoundaryDisplacementPoint_h

#include <model/ParticleInteraction/FieldPoint.h>


namespace model
{
    
    /******************************************************************************/
    template<typename DislocationNetworkType>
    struct BoundaryDisplacementPoint :
    /* inheritance */ public FieldPoint<BoundaryDisplacementPoint<DislocationNetworkType>,DislocationNetworkType::DisplacementField::dim,typename DislocationNetworkType::DisplacementField>
    {
        typedef typename DislocationNetworkType::DisplacementField DisplacementField;
        typedef typename DislocationNetworkType::ElementType::NodeType NodeType;

        constexpr static int dim=DisplacementField::dim;
        
        const Eigen::Matrix<double,dim,1>& P;
        const Eigen::Matrix<double,dim,1> S;
        
        BoundaryDisplacementPoint(const NodeType& node) :
        /*   */ P(node.p0)
        {
        
            std::cout<<"The S vector has not been computed yet yet"<<std::endl;

            
        }
        
    };
    
} // end namespace
#endif
