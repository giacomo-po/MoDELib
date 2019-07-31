/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_GlidePlaneBoundaryIO_H_
#define model_GlidePlaneBoundaryIO_H_

#include <tuple>
#include <iomanip>
#include <Eigen/Dense>
#include <MeshBoundarySegment.h>

namespace model
{
    
    template<short unsigned int dim>
    struct GlidePlaneBoundaryIO
    {
        
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        
        size_t glidePlaneID;          // sID
        VectorDim P0;          // position
        VectorDim P1;          // velocity
        VectorDim outNormal;          // velocity
        
        /**********************************************************************/
        GlidePlaneBoundaryIO(const size_t& glidePlaneID_in,
                             const MeshBoundarySegment<dim>& seg) :
        /* init */ glidePlaneID(glidePlaneID_in),
        /* init */ P0(seg.P0),
        /* init */ P1(seg.P1),
        /* init */ outNormal(seg.outNormal)
        {
         
//            assert(0 && "FINISH HERE, THIS MUST BE COMPATIBLE WITH ID READER");
            
        }
        
        /**********************************************************************/
        GlidePlaneBoundaryIO(const size_t& glidePlaneID_in,          // sID
                          const VectorDim& P0_in,          // position
                          const VectorDim& P1_in,
                          const VectorDim& n_in) :
        /* init */ glidePlaneID(glidePlaneID_in),
        /* init */ P0(P0_in),
        /* init */ P1(P1_in),
        /* init */ outNormal(n_in)
        {// Constructor for MicrostructureGenerator
        }
        
        /**********************************************************************/
        GlidePlaneBoundaryIO() :
        /* init */ glidePlaneID(0),
        /* init */ P0(VectorDim::Zero()),
        /* init */ P1(VectorDim::Zero()),
        /* init */ outNormal(VectorDim::Zero())
        {
        }

        
        /**********************************************************************/
        GlidePlaneBoundaryIO(std::stringstream& ss) :
        /* init */ glidePlaneID(0),
        /* init */ P0(VectorDim::Zero()),
        /* init */ P1(VectorDim::Zero()),
        /* init */ outNormal(VectorDim::Zero())
        {
            
            ss>>glidePlaneID;
            for(int d=0;d<dim;++d)
            {
                ss>>P0(d);
            }
            for(int d=0;d<dim;++d)
            {
                ss>>P1(d);
            }
            for(int d=0;d<dim;++d)
            {
                ss>>outNormal(d);
            }
        }

        
        /**********************************************************************/
        template <class T>
        friend T& operator << (T& os, const GlidePlaneBoundaryIO<dim>& ds)
        {
            os  << ds.glidePlaneID<<"\t"
            /**/<< std::setprecision(15)<<std::scientific<<ds.P0.transpose()<<"\t"
            /**/<< std::setprecision(15)<<std::scientific<<ds.P1.transpose()<<"\t"
            /**/<< std::setprecision(15)<<std::scientific<<ds.outNormal.transpose();
            return os;
        }
        
	};
	
}
#endif

