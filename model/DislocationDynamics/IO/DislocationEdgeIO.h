/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationEdgeIO_H_
#define model_DislocationEdgeIO_H_

#include <tuple>

namespace model
{
    
    template<short unsigned int dim>
    struct DislocationEdgeIO
    {
        
        
        size_t loopID;          // sID
        size_t sourceID;          // sID
        size_t sinkID;          // sID
        int  meshLocation;    // mesh location
        
        /**********************************************************************/
        template<typename LoopLinkType>
        DislocationEdgeIO(const LoopLinkType& ll) :
        /* init */ loopID(ll.loop()->sID),
        /* init */ sourceID(ll.source()->sID),
        /* init */ sinkID(ll.sink()->sID),
        /* init */ meshLocation(ll.pLink->meshLocation())
        {
         
//            assert(0 && "FINISH HERE, THIS MUST BE COMPATIBLE WITH ID READER");
            
        }
        
        /**********************************************************************/
        DislocationEdgeIO(const size_t& loopID_in,          // sID
                          const size_t& sourceID_in,          // position
                          const size_t& sinkID_in,          // velocity
                          const int& meshLocation_in) :
        /* init */ loopID(loopID_in),
        /* init */ sourceID(sourceID_in),
        /* init */ sinkID(sinkID_in),
        /* init */ meshLocation(meshLocation_in)
        {
            
            //            assert(0 && "FINISH HERE, THIS MUST BE COMPATIBLE WITH ID READER");
            
        }
        
        /**********************************************************************/
        DislocationEdgeIO() :
        /* init */ loopID(0),
        /* init */ sourceID(0),
        /* init */ sinkID(0),
        /* init */ meshLocation(0)
        {
            
            
        }
        
        /**********************************************************************/
        DislocationEdgeIO(std::stringstream& ss) :
        /* init */ loopID(0),
        /* init */ sourceID(0),
        /* init */ sinkID(0),
        /* init */ meshLocation(0)
        {
            ss>>loopID;
            ss>>sourceID;
            ss>>sinkID;
            ss>>meshLocation;
        }

        
        /**********************************************************************/
        template <class T>
        friend T& operator << (T& os, const DislocationEdgeIO<dim>& ds)
        {
            os  << ds.loopID<<"\t"
            /**/<< ds.sourceID<<"\t"
            /**/<< ds.sinkID<<"\t"
            /**/<< ds.meshLocation;
            return os;
        }
        
	};
	
}
#endif

