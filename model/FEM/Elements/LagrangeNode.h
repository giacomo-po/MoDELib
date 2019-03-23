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
	template<typename ElementType>
	struct LagrangeNode : public std::set<const ElementType*>
    {
        static constexpr int dim=ElementType::dim;

        typedef Eigen::Matrix<double,dim,1> PositionType;
        
        const PositionType P0;
        const size_t gID; // global ID
        
        /**********************************************************************/
        LagrangeNode(const PositionType& p,
                     const size_t& gid) :
        /* init list */ P0(p),
        /* init list */ gID(gid)
        {

        }
        
        /**********************************************************************/
        Eigen::Matrix<double,dim,1> outNormal() const
        {
            Eigen::Matrix<double,dim,1> temp(Eigen::Matrix<double,dim,1>::Zero());
            for(auto ele : *this)
            {
                const Eigen::Matrix<double,dim+1,1> bary(ele->simplex.pos2bary(P0));
                for(int k=0;k<dim+1;++k)
                {
                    if (std::fabs(bary(k))<FLT_EPSILON && ele->simplex.child(k).isBoundarySimplex())
                    {
                        temp += ele->simplex.nda.col(k).normalized();
                    }
                }
            }
            const double tempNorm(temp.norm());
            return tempNorm>FLT_EPSILON? (temp/tempNorm).eval() : Eigen::Matrix<double,dim,1>::Zero();
        }
        
    };
    
    
}	// close namespace
#endif