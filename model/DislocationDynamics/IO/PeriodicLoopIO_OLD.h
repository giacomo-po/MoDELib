/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PeriodicLoopIO_H_
#define model_PeriodicLoopIO_H_

#include <tuple>

namespace model
{
    
    template<short unsigned int dim>
    struct PeriodicLoopIO
    {
        
//        enum DislocationLoopType{GLISSILELOOP,SESSILELOOP,VIRTUALLOOP};

        
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        
         size_t sID;          // sID
//         VectorDim B;          // position
//         VectorDim N;          // velocity
//         VectorDim P;          // velocity
//         size_t grainID;          // component ID
//         int loopType;
//         long int periodicLoopID;          // sID
//         std::tuple<double,double,double> loopLength;
        
        /**********************************************************************/
        template<typename PeriodicLoopType>
        PeriodicLoopIO(const PeriodicLoopType& pL) :
        /* init */ sID(pL.sID)
//        /* init */,B(dL.flow().cartesian())
//        /* init */,N(dL.glidePlane? (dL.slippedArea()>FLT_EPSILON? dL.rightHandedUnitNormal() : dL.glidePlane->unitNormal) : VectorDim::Zero())
//        /* init */,P(dL.glidePlane? dL.glidePlane->P : dL.links().begin()->second->source()->get_P() )
//        /* init */,grainID(dL.grain.grainID)
//        /* init */,loopType(dL.loopType)
//        /* init */,periodicLoopID(dL.periodicLoop? dL.periodicLoop->sID : -1)
//        /* init */,loopLength(dL.network().outputLoopLength? dL.loopLength() : std::make_tuple(0.0,0.0,0.0))
        {
            
        }
        
        /**********************************************************************/
        PeriodicLoopIO(const size_t& sID_in         // sID
//                          const VectorDim& B_in,          // position
//                          const VectorDim& N_in,          // velocity
//                          const VectorDim& P_in,          // velocity
//                          const size_t& grainID_in,
//                          const int& loopType_in,
//                          const long int& periodicLoopID_in
                       ) :
        /* init */ sID(sID_in)
//        /* init */,B(B_in)
//        /* init */,N(N_in)
//        /* init */,P(P_in)
//        /* init */,grainID(grainID_in)
//        /* init */,loopType(loopType_in)
//        /* init */,periodicLoopID(periodicLoopID_in)
//        /* init */,loopLength(std::make_tuple(0.0,0.0,0.0))
        {
            
        }
        
        /**********************************************************************/
        PeriodicLoopIO() :
        /* init */ sID(0)
//        /* init */,B(VectorDim::Zero())
//        /* init */,N(VectorDim::Zero())
//        /* init */,P(VectorDim::Zero())
//        /* init */,grainID(0)
//        /* init */,loopType(0)
//        /* init */,periodicLoopID(-1)
//        /* init */,loopLength(std::make_tuple(0.0,0.0,0.0))
        {
            
        }
        
        /**********************************************************************/
        PeriodicLoopIO(std::stringstream& ss) :
        /* init */ sID(0)
//        /* init */,B(VectorDim::Zero())
//        /* init */,N(VectorDim::Zero())
//        /* init */,P(VectorDim::Zero())
//        /* init */,grainID(0)
//        /* init */,loopType(0)
//        /* init */,periodicLoopID(-1)
//        /* init */,loopLength(std::make_tuple(0.0,0.0,0.0))
        {
            ss>>sID;
//            for(int d=0;d<dim;++d)
//            {
//                ss>>B(d);
//            }
//            for(int d=0;d<dim;++d)
//            {
//                ss>>N(d);
//            }
//            for(int d=0;d<dim;++d)
//            {
//                ss>>P(d);
//            }
//            ss>>grainID;
//            ss>>loopType;
//            ss>>periodicLoopID;
//            double l1,l2,l3;
//            ss>>l1>>l2>>l3;
//            loopLength=std::make_tuple(l1,l2,l3);
        }
        
        /**********************************************************************/
        template <class T>
        friend T& operator << (T& os, const PeriodicLoopIO<dim>& ds)
        {
            os  << ds.sID<<"\t";
//            /**/<< std::setprecision(15)<<std::scientific<<ds.B.transpose()<<"\t"
//            /**/<< std::setprecision(15)<<std::scientific<<ds.N.transpose()<<"\t"
//            /**/<< std::setprecision(15)<<std::scientific<<ds.P.transpose()<<"\t"
//            /**/<< ds.grainID<<"\t"
//            /**/<< ds.loopType<<"\t"
//            /**/<< ds.periodicLoopID<<"\t"
//            /**/<< std::get<0>(ds.loopLength)<<"\t"<< std::get<1>(ds.loopLength)<<"\t"<< std::get<2>(ds.loopLength);
            return os;
        }
        
	};
	
}
#endif

