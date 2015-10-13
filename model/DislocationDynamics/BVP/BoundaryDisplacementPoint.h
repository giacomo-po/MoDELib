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
        
        const size_t gID;
        const Eigen::Matrix<double,dim,1> P;
        const Eigen::Matrix<double,dim,1> S;
        
        BoundaryDisplacementPoint(const NodeType& node,
                                  const Eigen::Matrix<double,dim,1>& s_in) :
        /*   */ gID(node.gID),
        /*   */ P(node.P0),
        /*   */ S(s_in)
        {
            
        }
        
        BoundaryDisplacementPoint(const size_t gID_in,
                                  const Eigen::Matrix<double,dim,1>& P0,
                                  const Eigen::Matrix<double,dim,1>& s_in) :
        /*   */ gID(gID_in),
        /*   */ P(P0),
        /*   */ S(s_in)
        {
            
        }
        
    };
    
} // end namespace
#endif
