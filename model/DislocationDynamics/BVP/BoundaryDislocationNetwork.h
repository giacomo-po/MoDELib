/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_BoundaryDislocationNetwork_H_
#define model_BoundaryDislocationNetwork_H_

#include <model/DislocationDynamics/StressStraight.h>

//#include <memory> // unique_ptr
//#include <model/FEM/FiniteElement.h>
//#include <model/DislocationDynamics/Materials/Material.h>
//#include <model/DislocationDynamics/BVP/DislocationNegativeFields.h>
//#include <model/Mesh/SimplicialMesh.h>
//#include <model/FEM/WeakForms/JGNselector.h>

namespace model
{
	
	template <int dim>
	class BoundaryDislocationNetwork : public std::deque<StressStraight<dim>>
    {
        
        typedef Eigen::Matrix<double,dim,1>     VectorDim;
        typedef Eigen::Matrix<double,dim,dim>   MatrixDim;
        
    public:
        
        /**********************************************************************/
        template <typename DislocationNetworkType>
        void update(const DislocationNetworkType& DN)
        {
            const double L=1000000;
            
            this->clear();
            
            for (typename DislocationNetworkType::NetworkLinkContainerType::const_iterator linkIter=DN.linkBegin();linkIter!=DN.linkEnd();++linkIter)
            {
                if(linkIter->second->is_boundarySegment())
                {
                    VectorDim n = linkIter->second->source->bndNormal()-linkIter->second->source->bndNormal().dot(linkIter->second->glidePlaneNormal)*linkIter->second->glidePlaneNormal;
                    VectorDim P0 = linkIter->second->source->get_P();
                    VectorDim P1 = P0 + n*L;
                    this->emplace_back(P0,P1,linkIter->second->Burgers);
                    
                    n = linkIter->second->sink->bndNormal()-linkIter->second->sink->bndNormal().dot(linkIter->second->glidePlaneNormal)*linkIter->second->glidePlaneNormal;
                    P0 = linkIter->second->sink->get_P() - n*L;
                    P1 = linkIter->second->sink->get_P();
                    this->emplace_back(P0,P1,linkIter->second->Burgers);
                }
            }
        }
        
        /**********************************************************************/
        MatrixDim stress(const VectorDim& x) const
        {
            MatrixDim temp(MatrixDim::Zero());
            for(int k=0;k<this->size();++k)
            {
                temp += this->operator[](k).stress(x);
            }
            return temp;
        }
        
	};
	
    
} // namespace model
#endif


