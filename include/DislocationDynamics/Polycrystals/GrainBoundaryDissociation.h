/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2016 by Nathaniel Burbery <nbur049@aucklanduni.ac.nz>
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_GrainBoundaryDissociation_H_
#define model_GrainBoundaryDissociation_H_

#include <deque>
#include <utility>
#include <math.h>
#include <TerminalColors.h>
#include <DissociateSegment.h>
 // defines mode::cout


namespace model
{
    
    template <typename DislocationNetworkType>
    class GrainBoundaryDissociation : private std::deque<DissociateSegment<typename DislocationNetworkType::LinkType> >
    {
        
        DislocationNetworkType& DN;
        
        static constexpr double chordTol=1.0;
        
    public:
        
        static bool use_GBdissociation;
        
        /**********************************************************************/
        GrainBoundaryDissociation(DislocationNetworkType& DN_in) :
        /* init */ DN(DN_in)
        {
            if(use_GBdissociation)
            {
                for (const auto& segment : DN.links() )
                {
                    
                    if ( !(segment.second.source->isBoundaryNode() && segment.second.sink->isBoundaryNode())
                        && segment.second.isGrainBoundarySegment()
                        && segment.second.chord().norm()>chordTol*DislocationNetworkRemesh<DislocationNetworkType>::Lmin)
                    {
                        this->emplace_back(segment.second);
                    }
                }
            }

        }
        
        /**********************************************************************/
        void dissociate()
        {
            if(use_GBdissociation)
            {
                std::cout<<"		grain-boundary dissociation..."<<std::flush;
                
                const auto t0= std::chrono::system_clock::now();
                size_t nDissociateted=0;
                for(const auto& gbSegment : *this)
                {
                    if(DN.link(gbSegment.sourceID,gbSegment.sinkID).first)
                    {
                        if(gbSegment.isValidDissociation)
                        {
                            DN.template disconnect<false>(gbSegment.sourceID,gbSegment.sinkID);
                            DN.connect(gbSegment.sourceID,gbSegment.sinkID,gbSegment.residualBurgers);
                            const size_t newNodeID(DN.insertVertex(gbSegment.dissociateMidpoint.cartesian(),gbSegment.grain.grainID).first->first);
                            DN.connect(gbSegment.sinkID,newNodeID,gbSegment.dissociateBurgers);
                            DN.connect(newNodeID,gbSegment.sourceID,gbSegment.dissociateBurgers);
                            nDissociateted++;
                        }
                    }
                }
                
                std::cout<<"("<<nDissociateted<<" dissociations)"<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
            }
        }
        
    };
    
    template <typename DislocationNetworkType>
    bool GrainBoundaryDissociation<DislocationNetworkType>::use_GBdissociation=false;
    
} // end namespace
#endif

