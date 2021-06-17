/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PeriodicLoopLinkIO_H_
#define model_PeriodicLoopLinkIO_H_

#include <tuple>
#include <iomanip>

namespace model
{
    
    template<short unsigned int dim>
    struct PeriodicLoopLinkIO
    {
        
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef Eigen::Matrix<double,dim-1,1> VectorLowerDim;

        size_t periodicLoopID;          // sID
        size_t sourceID;          // sID
        size_t sinkID;          // sID
//        bool isBoundary;
        
        /**********************************************************************/
        template<typename PeriodicLoopLinkType>
        PeriodicLoopLinkIO(const PeriodicLoopLinkType& link) :
        /* init */ periodicLoopID(link.periodicLoop.sID)
        /* init */,sourceID(link.source->sID)
        /* init */,sinkID(link.sink->sID)
//        /* init */,isBoundary(!link.twin)
        {
        }

        /**********************************************************************/
        PeriodicLoopLinkIO(const size_t& pID_in,          // position
                           const size_t& sourceID_in,
                           const size_t& sinkID_in,
                           const bool& isBnd) :
        /* init */ periodicLoopID(pID_in)
        /* init */,sourceID(sourceID_in)
        /* init */,sinkID(sinkID_in)
//        /* init */,isBoundary(isBnd)
        {
            
            //            assert(0 && "FINISH HERE, THIS MUST BE COMPATIBLE WITH ID READER");
            
        }
        
        /**********************************************************************/
        PeriodicLoopLinkIO() :
        /* init */ periodicLoopID(0)
        /* init */,sourceID(0)
        /* init */,sinkID(0)
//        /* init */,isBoundary(false)
        {
            
            
        }
        
        /**********************************************************************/
        PeriodicLoopLinkIO(std::stringstream& ss) :
        /* init */ periodicLoopID(0)
        /* init */,sourceID(0)
        /* init */,sinkID(0)
//        /* init */,isBoundary(false)
        {
            ss>>periodicLoopID;
            ss>>sourceID;
            ss>>sinkID;
//            ss>>isBoundary;
        }
        
        
        /**********************************************************************/
        template <class T>
        friend T& operator << (T& os, const PeriodicLoopLinkIO<dim>& ds)
        {
            os  << ds.periodicLoopID<<" "
            /**/<< ds.sourceID<<" "
            /**/<< ds.sinkID;
//            /**/<< ds.isBoundary;
            return os;
        }
        
	};
	
}
#endif

