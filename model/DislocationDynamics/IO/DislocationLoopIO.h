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
        
        const size_t loopID;          // loopID
        const VectorDim B;          // position
        const VectorDim N;          // velocity
        const VectorDim P;          // velocity
        const size_t grainID;          // component ID
        
        /**********************************************************************/
        template<typename DislocationLoopType>
        DislocationLoopIO(const DislocationLoopType& dL) :
        /* init */ loopID(dL.sID),
        /* init */ B(dL.flow().cartesian()),
        /* init */ N(dL.glidePlane.n.cartesian()),
        /* init */ P(dL.glidePlane.P.cartesian()),
        /* init */ grainID(dL.grain.grainID)
        {
            
        }
        
        /**********************************************************************/
        DislocationLoopIO(const size_t& loopID_in,         // loopID
                          const VectorDim& B_in,          // position
                          const VectorDim& N_in,          // velocity
                          const VectorDim& P_in,          // velocity
                          const size_t& grainID) :
        /* init */ loopID(loopID_in),
        /* init */ B(B_in),
        /* init */ N(N_in),
        /* init */ P(P_in),
        /* init */ grainID(grainID)
        {
            
        }

        
        /**********************************************************************/
        template <class T>
        friend T& operator << (T& os, const DislocationLoopIO<dim>& ds)
        {
            os  << ds.loopID<<"\t"
            /**/<< std::setprecision(15)<<std::scientific<<ds.B.transpose()<<"\t"
            /**/<< std::setprecision(15)<<std::scientific<<ds.N.transpose()<<"\t"
            /**/<< std::setprecision(15)<<std::scientific<<ds.P.transpose()<<"\t"
            /**/<< ds.grainID;

            return os;
        }
        
	};
	
}
#endif

