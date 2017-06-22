/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2011 by Benjamin Ramirez<ramirezbrf@gmail.com>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SPLINENODE_H_
#define model_SPLINENODE_H_

#include <model/LoopNetwork/LoopNode.h>
#include <model/Geometry/Splines/SplineNodeBase.h>


namespace model
{
	
	/**************************************************************************/
	/* class template SplineNodeBase: general case                            */
	/**************************************************************************/
	template <typename Derived,	short unsigned int dim, short unsigned int corder>
	class SplineNode
    {
		
	public:
//		SplineNode()
//        {
//			std::cout<<"This Node type has not been implemented"<<std::endl;
//			assert(0);
//		};
	};
    
    /**************************************************************************/
    /* class template SplineNodeBase: C1-continuous nodes                     */
    /**************************************************************************/
    template <typename Derived,	short unsigned int dim>
    class SplineNode<Derived,dim,1> :
    /* inherits */  public LoopNode<Derived>,
    /* inherits */  public SplineNodeBase<Derived,dim,1>
    
    {
        
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        
        constexpr static int corder=1;
        typedef SplineNodeBase<Derived,dim,1> SplineNodeBaseType;
        
    public:
        
        /*************************************************/
        SplineNode(const VectorDim& P_in) :
        /* init list */ SplineNodeBaseType(P_in)
        //        /* init list */ P(roundP(P_in))
        {
            
        }
        
        //		SplineNode()
        //        {
        //			std::cout<<"This Node type has not been implemented"<<std::endl;
        //			assert(0);
        //		};
    };
		
}

//#include "model/Geometry/Splines/SplineNode_Hermite.h"
//#include "model/Geometry/Splines/SplineNodeCatmullRom.h"

#endif
