/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PolyhedronInclusionEdgeIO_H_
#define model_PolyhedronInclusionEdgeIO_H_

#include <tuple>

namespace model
{
    
    struct PolyhedronInclusionEdgeIO
    {
        
//        typedef Eigen::Matrix<double,dim,1> VectorDim;
        
        
        size_t inclusionID;          // sID
        size_t faceID;          // sID
        size_t sourceID;          // sID
        size_t sinkID;          // sID

        
        /**********************************************************************/
        PolyhedronInclusionEdgeIO(const size_t& inclusionID_in,          // sID
                        const size_t& faceID_in,
                        const size_t& sourceID_in,
                        const size_t& sinkID_in) :
        /* init */ inclusionID(inclusionID_in),
        /* init */ faceID(faceID_in),
        /* init */ sourceID(sourceID_in),
        /* init */ sinkID(sinkID_in)
        {
            
            //            assert(0 && "FINISH HERE, THIS MUST BE COMPATIBLE WITH ID READER");
            
        }
        
        /**********************************************************************/
        PolyhedronInclusionEdgeIO() :
        /* init */ inclusionID(0),
        /* init */ faceID(0),
        /* init */ sourceID(0),
        /* init */ sinkID(0)
        {
            
            
        }
        
        /**********************************************************************/
        PolyhedronInclusionEdgeIO(std::stringstream& ss) :
        /* init */ inclusionID(0),
        /* init */ faceID(0),
        /* init */ sourceID(0),
        /* init */ sinkID(0)
        {
            ss>>inclusionID;
            ss>>faceID;
            ss>>sourceID;
            ss>>sinkID;
        }

        
        /**********************************************************************/
        template <class T>
        friend T& operator << (T& os, const PolyhedronInclusionEdgeIO& ds)
        {
            os  << ds.inclusionID<<"\t"
            /**/<< ds.faceID<<"\t"
            /**/<< ds.sourceID<<"\t"
            /**/<< ds.sinkID;
            return os;
        }
        
	};
	
}
#endif

