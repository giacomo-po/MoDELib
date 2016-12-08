/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2016 by Nathaniel Burbery <nbur049@aucklanduni.ac.nz>
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_GrainBoundaryTransmission_H_
#define model_GrainBoundaryTransmission_H_

#include <deque>
#include <utility>
#include <math.h>
#include <model/Utilities/TerminalColors.h>
#include <model/DislocationDynamics/Polycrystals/TransmitSegment.h>
#include <model/MPI/MPIcout.h> // defines mode::cout


namespace model
{
    
    
    
    template <typename DislocationNetworkType>
    class GrainBoundaryTransmission : private std::deque<TransmitSegment<typename DislocationNetworkType::LinkType> >
    //    :
    //    /* base */ private std::map<int,const Grain<dim>* const>,
    //    /* base */ private std::map<int,LatticePlane>
    {
        
        DislocationNetworkType& DN;
        
        static constexpr double chordTol=1.0;
        
    public:
        
        static bool use_GBtransmission;

        /**********************************************************************/
        GrainBoundaryTransmission(DislocationNetworkType& DN_in) :
        /* init */ DN(DN_in)
        {
            
            for (const auto& segment : DN.links() )
            {
                
                if ( !(segment.second.source->isBoundaryNode() && segment.second.sink->isBoundaryNode())
                    && segment.second.isGrainBoundarySegment()
                    && segment.second.chord().norm()>chordTol*DislocationNetworkRemesh<DislocationNetworkType>::Lmin)
                {
                    //Lmin=DislocationNetworkRemesh<DislocationNetworkType>::Lmin;
                    this->emplace_back(segment.second);
                }
            }
            
        }
        
        /**********************************************************************/
        size_t transmit()
        {
            model::cout<<"		grain-boundary transmission..."<<std::flush;
            
            const auto t0= std::chrono::system_clock::now();
            size_t nTransmitted=0;
            for(const auto& gbSegment : *this)
            {
                if(DN.link(gbSegment.sourceID,gbSegment.sinkID).first)
                {
                 if(gbSegment.isValidTransmission)
                 {
                     DN.template disconnect<false>(gbSegment.sourceID,gbSegment.sinkID);
                     DN.connect(gbSegment.sourceID,gbSegment.sinkID,gbSegment.residualBurgers);
                     const size_t newNodeID(DN.insertVertex(gbSegment.transmitMidpoint.cartesian(),gbSegment.grain.grainID).first->first);
                     DN.connect(gbSegment.sinkID,newNodeID,gbSegment.transmitBurgers);
                     DN.connect(newNodeID,gbSegment.sourceID,gbSegment.transmitBurgers);
                     nTransmitted++;
                 }
                }
            }
            
            
            
            model::cout<<"("<<nTransmitted<<" transmissions)"<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
            
            
            return nTransmitted;
        }
        
    };
    
    template <typename DislocationNetworkType>
    bool GrainBoundaryTransmission<DislocationNetworkType>::use_GBtransmission=false;

    
} // end namespace
#endif

