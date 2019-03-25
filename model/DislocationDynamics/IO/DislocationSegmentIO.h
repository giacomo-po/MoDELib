/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationSegmentIO_H_
#define model_DislocationSegmentIO_H_

#include <tuple>
#include <iomanip>

namespace model
{
    
    template<short unsigned int dim>
    struct DislocationSegmentIO
    {
        
                typedef Eigen::Matrix<double,dim,1> VectorDim;
//        size_t loopID;          // sID
        const size_t sourceID;          // sID
        const size_t sinkID;          // sID
//        int  meshLocation;    // mesh location
        
        VectorDim b;
        VectorDim n;

        int  meshLocation;
//        int loopCounter;
//        bool isGlissile
        
//        /**********************************************************************/
//        template<typename LoopLinkType>
//        DislocationSegmentIO(const LoopLinkType& ll) :
//        /* init */ loopID(ll.loop()->sID),
//        /* init */ sourceID(ll.source()->sID),
//        /* init */ sinkID(ll.sink()->sID),
//        /* init */ meshLocation(ll.pLink->meshLocation())
//        {
//         
////            assert(0 && "FINISH HERE, THIS MUST BE COMPATIBLE WITH ID READER");
//            
//        }
        
        /**********************************************************************/
        DislocationSegmentIO(const size_t& sourceID_in,          // position
                          const size_t& sinkID_in) :
        /* init */ sourceID(sourceID_in),
        /* init */ sinkID(sinkID_in),
        /* init */ b(VectorDim::Zero()),
        /* init */ n(VectorDim::Zero()),
        /* init */ meshLocation(-1)
//        loopCounter(0),
//        isGlissile(true)
//        /* init */ meshLocation(meshLocation_in)
        {
            
            //            assert(0 && "FINISH HERE, THIS MUST BE COMPATIBLE WITH ID READER");
            
        }
        
//        /**********************************************************************/
//        DislocationSegmentIO() :
//        /* init */ loopID(0),
//        /* init */ sourceID(0),
//        /* init */ sinkID(0),
//        /* init */ meshLocation(0)
//        {
//            
//            
//        }
        
//        /**********************************************************************/
//        DislocationSegmentIO(std::stringstream& ss) :
//        /* init */ loopID(0),
//        /* init */ sourceID(0),
//        /* init */ sinkID(0),
//        /* init */ meshLocation(0)
//        {
//            ss>>loopID;
//            ss>>sourceID;
//            ss>>sinkID;
//            ss>>meshLocation;
//        }

        
        /**********************************************************************/
        template <class T>
        friend T& operator << (T& os, const DislocationSegmentIO<dim>& ds)
        {
            os  << ds.sourceID<<" "
            /**/<< ds.sinkID<<" "
            /**/<< std::setprecision(15)<<std::scientific<<ds.b.transpose()<<" "
            /**/<< std::setprecision(15)<<std::scientific<<ds.n.transpose()<<" "
            /**/<< ds.meshLocation;
            return os;
        }
        
	};
	
}
#endif

