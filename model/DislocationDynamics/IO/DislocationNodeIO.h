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

namespace model
{
    
    template<short unsigned int dim>
    struct DislocationNodeIO
    {
        
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        
        const size_t sID;          // sID
        const VectorDim P;          // position
        const VectorDim V;          // velocity
        const double f;             // velocity reduction factor
        const size_t snID;          // component ID
        const bool  meshLocaiton;   // mesh location
        
        /**********************************************************************/
        template<typename DislocationNodeType>
        DislocationNodeIO(const DislocationNodeType& dn) :
        /* init */ sID(dn.sID),
        /* init */ P(dn.get_P()),
        /* init */ V(dn.get_V()),
        /* init */ f(dn.velocityReduction()),
        /* init */ snID(dn.pSN()->sID),
        /* init */ meshLocaiton(dn.meshLocation()==onMeshBoundary)
        {
         
            assert(0 && "FINISH HERE, THIS MUST BE COMPATIBLE WITH ID READER");
            
        }

        
	};
	
}
#endif

