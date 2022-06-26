/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationCrossSlip_h_
#define model_DislocationCrossSlip_h_

#include <memory>

#include <DislocationDynamicsModule.h>
#include <CrossSlipModels.h>

#ifndef NDEBUG
#define VerboseCrossSlip(N,x) if(verboseCrossSlip>=N){std::cout<<x;}
#else
#define VerboseCrossSlip(N,x)
#endif

namespace model
{
        
    template <typename DislocationNetworkType>
    class DislocationCrossSlip
    {
        typedef TypeTraits<DislocationNetworkType> TraitsType;
        static constexpr int dim=TraitsType::dim;
        typedef typename TraitsType::NetworkLinkType NetworkLinkType;
        typedef typename TraitsType::NetworkNodeType NetworkNodeType;
        typedef typename TraitsType::LoopNodeType LoopNodeType;
        typedef typename TraitsType::VectorDim VectorDim;
        typedef std::tuple<std::shared_ptr<NetworkNodeType>,std::shared_ptr<NetworkNodeType>,size_t,size_t> CrossSlipTupleType;
        typedef std::deque<CrossSlipTupleType> CrossSlipContainerType;
        
        DislocationNetworkType& DN; //! A reference to the DislocationNetwork
        CrossSlipContainerType crossSlipDeq;
        
        static std::shared_ptr<BaseCrossSlipModel<DislocationNetworkType>> getModel(const PolycrystallineMaterialBase& material,const int& csModel);
        
    public:
        
        const int crossSlipModel;
        const std::shared_ptr<BaseCrossSlipModel<DislocationNetworkType>> model;
        const double crossSlipDeg;
        const int verboseCrossSlip;

        DislocationCrossSlip(DislocationNetworkType& DN_in) ;
        void findCrossSlipSegments();
        void execute();

        
    };
        
}
#endif
