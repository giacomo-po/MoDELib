/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2016 by Nathaniel Burbery <nbur049@aucklanduni.ac.nz>
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_GrainBoundary_H_
#define model_GrainBoundary_H_

#include <deque>
#include <utility>
#include <math.h>
#include <model/Utilities/TerminalColors.h>
#include <model/DislocationDynamics/Polycrystals/TransmitSegment.h>



namespace model
{
    
    
    
    template <template DislocationNetworkType>
    class GrainBoundaryTransmission : private std::deque<TransmitSegment<typename DislocationNetworkType::LinkType> >
//    :
//    /* base */ private std::map<int,const Grain<dim>* const>,
//    /* base */ private std::map<int,LatticePlane>
    {

        DislocationNetworkType& DN;
        
    public:
        
        /**********************************************************************/
        GrainBoundaryTransmission(DislocationNetworkType& DN_in) :
        /* init */ DN(DN_in)
        {
        
            for (const auto& segment : DN.links() )
            {
                
                if (	!(segment.second.source->isBoundaryNode() && segment.second.sink->isBoundaryNode())
                    && segment.second.isGrainBoundarySegment()
                    && segment.second.chord().norm()>chordTol*DislocationNetworkRemesh<DislocationNetworkType>::Lmin)
                {
                    Lmin=DislocationNetworkRemesh<DislocationNetworkType>::Lmin;
                    this->emplace_back(segment.second);
                }
            }
            
        }
        
        /**********************************************************************/
        size_t transmit()
        {
        
        }
        
    };
    
} // end namespace
#endif

