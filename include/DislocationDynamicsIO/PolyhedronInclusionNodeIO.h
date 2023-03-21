/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PolyhedronInclusionNodeIO_H_
#define model_PolyhedronInclusionNodeIO_H_

#include <tuple>

namespace model
{
    
    template<short unsigned int dim>
    struct PolyhedronInclusionNodeIO
    {
        
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        
//        size_t inclusionID;
        size_t nodeID;          // sID
        VectorDim P;
        
        
        /**********************************************************************/
        PolyhedronInclusionNodeIO(
                                  //const size_t& inclusionID_in,
                                  const size_t& nodeID_in,          // sID
                   const VectorDim& P_in) :
//        /* init */ inclusionID(inclusionID_in),
        /* init */ nodeID(nodeID_in),
        /* init */ P(P_in)
        {
            
            //            assert(0 && "FINISH HERE, THIS MUST BE COMPATIBLE WITH ID READER");
            
        }
        
        /**********************************************************************/
        PolyhedronInclusionNodeIO() :
//        /* init */ inclusionID(0),
        /* init */ nodeID(0),
        /* init */ P(VectorDim::Zero())
        {
            
            
        }
        
        /**********************************************************************/
        PolyhedronInclusionNodeIO(std::stringstream& ss) :
//        /* init */ inclusionID(0),
        /* init */ nodeID(0),
        /* init */ P(VectorDim::Zero())
        {
//            ss>>inclusionID;
            ss>>nodeID;
            for(int d=0;d<dim;++d)
            {
                ss>>P(d);
            }
        }

        
        /**********************************************************************/
        template <class T>
        friend T& operator << (T& os, const PolyhedronInclusionNodeIO<dim>& ds)
        {
            os  //<< ds.inclusionID<<"\t"
            <<ds.nodeID<<"\t"
            /**/<< ds.P.transpose();
            return os;
        }
        
	};
	
}
#endif

