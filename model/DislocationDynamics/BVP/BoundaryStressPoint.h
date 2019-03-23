/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2012 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _model_BoundaryStressPoint_h
#define _model_BoundaryStressPoint_h

#include <model/DislocationDynamics/Materials/Material.h>
#include <model/ParticleInteraction/FieldPoint.h>
#include <model/FEM/Domains/BoundaryIntegrationPoint.h>


namespace model
{
    
    /******************************************************************************/
    template<typename DislocationNetworkType>
    struct BoundaryStressPoint :
    /* inheritance */ public FieldPoint<BoundaryStressPoint<DislocationNetworkType>,DislocationNetworkType::StressField::dim,typename DislocationNetworkType::StressField>,
    /* inheritance */ public BoundaryIntegrationPoint<typename DislocationNetworkType::ElementType>
    {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        
        typedef typename DislocationNetworkType::ElementType ElementType;
        typedef BoundaryIntegrationPoint<ElementType> BoundaryIntegrationPointType;
        typedef typename DislocationNetworkType::StressField StressField;
        constexpr static int dim=StressField::dim;
        typedef FieldPoint<BoundaryStressPoint<DislocationNetworkType>,dim,StressField> FieldPointType;
        
        const Eigen::Matrix<double,dim,1> P;

        /**********************************************************************/
        BoundaryStressPoint(const ElementType& _ele,
        /*               */ const Eigen::Matrix<double,dim+1,1>& _domainBary,
        /*               */ const int& _boundaryFace,
        /*               */ const double& _weight) :
        /*               */ FieldPointType(true),
        /*               */ BoundaryIntegrationPointType(_ele,_domainBary,_boundaryFace,_weight),
        /*               */ P(_ele.position(_domainBary))
        {
            
        }
        
        /**********************************************************************/
        Eigen::Matrix<double,dim,dim> operator()() const
        {
            return this->template field<StressField>();
        }
        
    };
    
} // end namespace
#endif
