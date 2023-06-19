/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationNodeIO_H_
#define model_DislocationNodeIO_H_

#include <tuple>
#include <iomanip>
#include <Eigen/Dense>


namespace model
{
    
    template<short unsigned int dim>
    struct DislocationNodeIO
    {
        
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        
        size_t sID;          // sID
        VectorDim P;          // position
        VectorDim V;          // velocity
        double velocityReduction;             // velocity reduction factor
        int  meshLocation;    // mesh location
//        size_t masterID;          // sID

        
        /**********************************************************************/
        template<typename DislocationNodeType>
        DislocationNodeIO(const DislocationNodeType& dn) :
        /* init */ sID(dn.sID),
        /* init */ P(dn.get_P()),
        /* init */ V(dn.get_V()),
        /* init */ velocityReduction(dn.velocityReduction()),
//        /* init */ snID(dn.pSN()->sID),
        /* init */ meshLocation(dn.meshLocation())
//        /* init */ masterID(dn.masterNode? dn.masterNode->sID : sID)
        {
         
            
        }
        
        /**********************************************************************/
        DislocationNodeIO(const size_t& sID_in,          // sID
                          const VectorDim& P_in,          // position
                          const VectorDim& V_in,          // velocity
                          const double& velocityReduction_in,             // velocity reduction factor
//                          const size_t& snID_in,          // component ID
                          const int&  meshLocation_in) :
        /* init */ sID(sID_in),
        /* init */ P(P_in),
        /* init */ V(V_in),
        /* init */ velocityReduction(velocityReduction_in),
//        /* init */ snID(snID_in),
        /* init */ meshLocation(meshLocation_in)
//        /* init */ masterID(sID)
        {// Constructor for MicrostructureGenerator
        }
        
        /**********************************************************************/
        DislocationNodeIO() :
        /* init */ sID(0),
        /* init */ P(VectorDim::Zero()),
        /* init */ V(VectorDim::Zero()),
        /* init */ velocityReduction(1.0),
//        /* init */ snID(0),
        /* init */ meshLocation(0)
//        /* init */ masterID(0)
        {
        }

        /**********************************************************************/
        DislocationNodeIO(std::stringstream& ss) :
        /* init */ sID(0),
        /* init */ P(VectorDim::Zero()),
        /* init */ V(VectorDim::Zero()),
        /* init */ velocityReduction(1.0),
//        /* init */ snID(0),
        /* init */ meshLocation(0)
//        /* init */ masterID(0)
        {
            ss>>sID;
            for(int d=0;d<dim;++d)
            {
                ss>>P(d);
            }
            for(int d=0;d<dim;++d)
            {
                ss>>V(d);
            }
            ss>>velocityReduction;
//            ss>>snID;
            ss>>meshLocation;
//            ss>>masterID;
        }

        
        /**********************************************************************/
        template <class T>
        friend T& operator << (T& os, const DislocationNodeIO<dim>& ds)
        {
            os  << ds.sID<<"\t"
            /**/<< std::setprecision(15)<<std::scientific<<ds.P.transpose()<<"\t"
            /**/<< ds.V.transpose()<<"\t"
            /**/<< ds.velocityReduction<<"\t"
//            /**/<< ds.snID<<"\t"
            /**/<< ds.meshLocation;
//            /**/<< ds.masterID;
            return os;
        }
        
	};
	
}
#endif

