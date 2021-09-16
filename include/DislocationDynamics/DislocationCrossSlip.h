/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationCrossSlip_h_
#define model_DislocationCrossSlip_h_

#include <DislocationDynamicsModule.h>

//#include <utility> // for std::pair
//#include <vector>
//#include <Eigen/Dense>
//#include <DislocationNetworkTraits.h>
////#include <DislocationDynamicsModule.h>
//#include <SegmentSegmentDistance.h>
////#include <DislocationSegmentIntersection.h>
//// #include <DislocationNetworkRemesh.h>
//#include <CrossSlipModels.h>
////#include <DislocationNetwork.h>
//#include <PlanePlaneIntersection.h>
//
//
////#include <EqualIteratorRange.h>
////#include <N2IteratorRange.h>
//#include <BCClattice.h>
//#include <FCClattice.h>


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
        //typedef DislocationNetwork<dim,corder> DislocationNetworkType;
        typedef typename TraitsType::NetworkLinkType NetworkLinkType;
        typedef typename TraitsType::NetworkNodeType NetworkNodeType;
        typedef typename TraitsType::LoopNodeType LoopNodeType;
        // typedef typename DislocationNetworkType::IsNetworkEdgeType IsNetworkLinkType;
        // typedef typename DislocationNetworkType::IsNodeType IsNodeType;
        typedef typename TraitsType::VectorDim VectorDim;
        // typedef typename DislocationNetworkType::NetworkLinkContainerType NetworkLinkContainerType;
        
//        typedef std::tuple<size_t,size_t,size_t,size_t> CrossSlipTupleType;
        typedef std::tuple<std::shared_ptr<NetworkNodeType>,std::shared_ptr<NetworkNodeType>,size_t,size_t> CrossSlipTupleType;

        typedef std::deque<CrossSlipTupleType> CrossSlipContainerType;
        
        //! A reference to the DislocationNetwork
        DislocationNetworkType& DN;
        CrossSlipContainerType crossSlipDeq;
        
    public:
        
        const int verboseCrossSlip;
        const double crossSlipDeg;
        
        DislocationCrossSlip(DislocationNetworkType& DN_in) ;
        void findCrossSlipSegments();
        void execute();

        
    };
    
//    template <typename DislocationNetworkType>
//    int DislocationCrossSlip<DislocationNetworkType>::verboseCrossSlip=0;
//
//    template <typename DislocationNetworkType>
//    double DislocationCrossSlip<DislocationNetworkType>::crossSlipDeg=2.0;
    
}
#endif
