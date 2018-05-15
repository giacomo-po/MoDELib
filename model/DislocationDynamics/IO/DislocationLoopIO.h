/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationLoopIO_H_
#define model_DislocationLoopIO_H_

#include <tuple>

namespace model
{
    
    template<short unsigned int dim>
    struct DislocationLoopIO
    {
        
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        
         size_t sID;          // sID
         VectorDim B;          // position
         VectorDim N;          // velocity
         VectorDim P;          // velocity
         size_t grainID;          // component ID
         std::tuple<double,double,double> loopLength;
        
        /**********************************************************************/
        template<typename DislocationLoopType>
        DislocationLoopIO(const DislocationLoopType& dL) :
        /* init */ sID(dL.sID),
        /* init */ B(dL.flow().cartesian()),
        /* init */ N(dL.glidePlane.unitNormal),
        /* init */ P(dL.glidePlane.P),
        /* init */ grainID(dL.grain.grainID),
        /* init */ loopLength(dL.network().outputLoopLength? dL.loopLength() : std::make_tuple(0.0,0.0,0.0))
        {
            
        }
        
        /**********************************************************************/
        DislocationLoopIO(const size_t& sID_in,         // sID
                          const VectorDim& B_in,          // position
                          const VectorDim& N_in,          // velocity
                          const VectorDim& P_in,          // velocity
                          const size_t& grainID) :
        /* init */ sID(sID_in),
        /* init */ B(B_in),
        /* init */ N(N_in),
        /* init */ P(P_in),
        /* init */ grainID(grainID),
        /* init */ loopLength(std::make_tuple(0.0,0.0,0.0))
        {
            
        }
        
        /**********************************************************************/
        DislocationLoopIO() :
        /* init */ sID(0),
        /* init */ B(VectorDim::Zero()),
        /* init */ N(VectorDim::Zero()),
        /* init */ P(VectorDim::Zero()),
        /* init */ grainID(0),
        /* init */ loopLength(std::make_tuple(0.0,0.0,0.0))
        {
            
        }
        
        /**********************************************************************/
        DislocationLoopIO(std::stringstream& ss) :
        /* init */ sID(0),
        /* init */ B(VectorDim::Zero()),
        /* init */ N(VectorDim::Zero()),
        /* init */ P(VectorDim::Zero()),
        /* init */ grainID(0),
        /* init */ loopLength(std::make_tuple(0.0,0.0,0.0))
        {
            ss>>sID;
            for(int d=0;d<dim;++d)
            {
                ss>>B(d);
            }
            for(int d=0;d<dim;++d)
            {
                ss>>N(d);
            }
            for(int d=0;d<dim;++d)
            {
                ss>>P(d);
            }
            ss>>grainID;
            double l1,l2,l3;
            ss>>l1>>l2>>l3;
            loopLength=std::make_tuple(l1,l2,l3);
        }
        
        /**********************************************************************/
        template <class T>
        friend T& operator << (T& os, const DislocationLoopIO<dim>& ds)
        {
            os  << ds.sID<<"\t"
            /**/<< std::setprecision(15)<<std::scientific<<ds.B.transpose()<<"\t"
            /**/<< std::setprecision(15)<<std::scientific<<ds.N.transpose()<<"\t"
            /**/<< std::setprecision(15)<<std::scientific<<ds.P.transpose()<<"\t"
            /**/<< ds.grainID<<"\t"
            /**/<< std::get<0>(ds.loopLength)<<"\t"<< std::get<1>(ds.loopLength)<<"\t"<< std::get<2>(ds.loopLength);
            return os;
        }
        
	};
	
}
#endif

