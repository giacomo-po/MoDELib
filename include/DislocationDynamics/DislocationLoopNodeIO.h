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

        size_t loopID;
        size_t sID;          // sID
        VectorDim P;          // position
        size_t networkNodeID;          // sID
        VectorDim periodicShift;
        short int edgeID;          // sID

        
        /**********************************************************************/
        template<typename DislocationLoopNodeType>
        DislocationLoopNodeIO(const DislocationLoopNodeType& dn) :
        /* init */ loopID(dn.loop()->sID)
        /* init */,sID(dn.sID)
        /* init */,P(dn.get_P())
        /* init */,networkNodeID(dn.networkNode->sID)
        /* init */,periodicShift(dn.periodicPlanePatch()? dn.periodicPlanePatch()->shift : VectorDim::Zero())
        /* init */,edgeID(dn.periodicPlaneEdge? dn.periodicPlaneEdge->edgeID : -1)
        {
         
            
        }
        
        /**********************************************************************/
        DislocationLoopNodeIO(const size_t& sID_in,          // sID
                                const size_t& loopID_in,
                          const VectorDim& P_in,          // position
                          const size_t& nnID,
                            const VectorDim& shift,
                              const short int& edgeID_in) :
        /* init */ loopID(loopID_in),
        /* init */ sID(sID_in),
        /* init */ P(P_in),
        /* init */ networkNodeID(nnID),
        /* init */ periodicShift(shift),
        /* init */ edgeID(edgeID_in)
        {// Constructor for MicrostructureGenerator
        }
        
        /**********************************************************************/
        DislocationLoopNodeIO() :
        /* init */ loopID(0),
        /* init */ sID(0),
        /* init */ P(VectorDim::Zero()),
        /* init */ networkNodeID(0),
        /* init */ periodicShift(VectorDim::Zero()),
        /* init */ edgeID(-1)
        {
        }

        /**********************************************************************/
        DislocationLoopNodeIO(std::stringstream& ss) :
        /* init */ loopID(0),
        /* init */ sID(0),
        /* init */ P(VectorDim::Zero()),
        /* init */ networkNodeID(0),
        /* init */ periodicShift(VectorDim::Zero()),
        /* init */ edgeID(-1)
        {
            ss>>loopID;
            ss>>sID;
            for(int d=0;d<dim;++d)
            {
                ss>>P(d);
            }
            ss>>networkNodeID;
            for(int d=0;d<dim;++d)
            {
                ss>>periodicShift(d);
            }
            ss>>edgeID;
        }

        
        /**********************************************************************/
        template <class T>
        friend T& operator << (T& os, const DislocationLoopNodeIO<dim>& ds)
        {
            os  << ds.loopID<<"\t"
            /**/<< ds.sID<<"\t"
            /**/<< std::setprecision(15)<<std::scientific<<ds.P.transpose()<<"\t"
            /**/<< ds.networkNodeID<<"\t"
            /**/<< std::setprecision(15)<<std::scientific<<ds.periodicShift.transpose()<<"\t"
            /**/<< ds.edgeID;
            return os;
        }
        
	};
	
}
#endif

