/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationLoopNodeIO_H_
#define model_DislocationLoopNodeIO_H_

#include <tuple>
#include <iomanip>
#include <limits>
#include <Eigen/Dense>


namespace model
{
    
    template<short unsigned int dim>
    struct DislocationLoopNodeIO
    {
        
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        
        size_t sID;          // sID
        size_t loopID;
        VectorDim P;          // position
        size_t networkNodeID;          // sID

        
        /**********************************************************************/
        template<typename DislocationLoopNodeType>
        DislocationLoopNodeIO(const DislocationLoopNodeType& dn) :
        /* init */ sID(dn.sID),
        /* init */ loopID(dn.loop()? dn.loop()->sID : std::numeric_limits<size_t>::max),
        /* init */ P(dn.get_P()),
        /* init */ networkNodeID(dn.networkNode? dn.networkNode->sID : std::numeric_limits<size_t>::max)
        {
         
            
        }
        
        /**********************************************************************/
        DislocationLoopNodeIO(const size_t& sID_in,          // sID
                                const size_t& loopID_in,
                          const VectorDim& P_in,          // position
                          const size_t& nnID) :
        /* init */ sID(sID_in),
        /* init */ loopID(loopID_in),
        /* init */ P(P_in),
        /* init */ networkNodeID(nnID)
        {// Constructor for MicrostructureGenerator
        }
        
        /**********************************************************************/
        DislocationLoopNodeIO() :
        /* init */ sID(0),
        /* init */ loopID(0),
        /* init */ P(VectorDim::Zero()),
        /* init */ networkNodeID(0)
        {
        }

        /**********************************************************************/
        DislocationLoopNodeIO(std::stringstream& ss) :
        /* init */ sID(0),
        /* init */ loopID(0),
        /* init */ P(VectorDim::Zero()),
        /* init */ networkNodeID(0)
        {
            ss>>sID;
            ss>>loopID;
            for(int d=0;d<dim;++d)
            {
                ss>>P(d);
            }
            ss>>networkNodeID;
        }

        
        /**********************************************************************/
        template <class T>
        friend T& operator << (T& os, const DislocationLoopNodeIO<dim>& ds)
        {
            os  << ds.sID<<"\t"
            /**/<< ds.loopID<<"\t"
            /**/<< std::setprecision(15)<<std::scientific<<ds.P.transpose()<<"\t"
            /**/<< ds.networkNodeID;
            return os;
        }
        
	};
	
}
#endif

