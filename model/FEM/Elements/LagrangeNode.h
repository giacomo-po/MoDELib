/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LagrangeNode_H_
#define model_LagrangeNode_H_

#include <set>
#include <Eigen/Dense>

namespace model
{
    
	
    /**************************************************************************/
	/**************************************************************************/
//	template<int _dim>
	template<typename ElementType>
	struct LagrangeNode : public std::set<const ElementType*>
//	struct LagrangeNode : public std::vector<const ElementType*>
    {
        static constexpr int dim=ElementType::dim;
//        static constexpr int dim=_dim;

        typedef Eigen::Matrix<double,dim,1> PositionType;
        
        const PositionType P0;
        const size_t gID; // global ID
        
        /**********************************************************************/
        LagrangeNode(const PositionType& p,
                     const size_t& gid) :
//                     const ElementType* const pEle) :
        /* init list */ P0(p),
        /* init list */ gID(gid)
        {
//            std::cout<<"Lagrange Node constructor..."<<this<<" ";
////            auto success=this->insert(pEle);
//            const bool success=this->emplace(pEle).second;
//            assert(success && "COULD NOT INSERT ELEMENT IN LAGRANGE NODE.");
//            std::cout<<this->size()<<std::endl;
//            show();
        }
        
//        void show() const
//        {
//            for (auto ele : *this)
//            {
//                std::cout<<"node "<<gID<<" ("<<this<< ") of element "<<ele<<std::endl;
//            }
//        }
        
    };
    
    
}	// close namespace
#endif