/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2018 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2018 by Yinan Cui <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationLoopPatches_cpp_
#define model_DislocationLoopPatches_cpp_


#include <DislocationLoopPatches.h>

namespace model
{


template <int dim>
DislocationLoopPatches<dim>::DislocationLoopPatches(const std::shared_ptr<PeriodicGlidePlaneType>& periodicGlidePlane_in) :
/* */ periodicGlidePlane(periodicGlidePlane_in)
{
    
}

template <int dim>
void DislocationLoopPatches<dim>::update(const std::vector<VectorDim>& linkShifts, const std::vector<Eigen::Matrix<double,dim,1>>& nodePos)
{
    if(periodicGlidePlane)
    {
        LocalPatchPositionType::clear();
        GlobalPatchPositionType::clear();
        std::vector<Eigen::Matrix<double,dim-1,1>> localNodePos;
        for(const auto& globalPos : nodePos)
        {
            localNodePos.push_back(periodicGlidePlane->referencePlane->localPosition(globalPos));
        }
        
        const auto loopPatches(periodicGlidePlane->filledPatches(linkShifts));
        for(const auto& patch : loopPatches)
        {
            // Collect local patch  nodes pos
            std::vector<Eigen::Matrix<double,dim-1,1>> localPatchPos;
            for(const auto& edge : patch->edges())
            {
                localPatchPos.push_back(*edge->source);
            }
            
            const auto localLoopPosOnPeriodicPlane(SutherlandHodgman::clip(localNodePos,localPatchPos));
            std::vector<Eigen::Matrix<double,dim,1>> globalLoopPos;
            std::vector<Eigen::Matrix<double,dim-1,1>> localLoopPos;
            
            for(const auto& localPos : localLoopPosOnPeriodicPlane)
            {
                //                        const VectorDim globalPos(periodicGlidePlane->referencePlane->globalPosition(localPos)+patch->shift);
                
                globalLoopPos.push_back(periodicGlidePlane->referencePlane->globalPosition(localPos)+patch->shift);
                localLoopPos.push_back(patch->glidePlane->localPosition(globalLoopPos.back()));
                //
                //
                //                        localLoopPosOnPatchPlane.push_back(patch->glidePlane->localPosition(globalPos));
            }
            GlobalPatchPositionType::emplace(patch,globalLoopPos);
            LocalPatchPositionType::emplace(patch,localLoopPos);
            
        }
    }
}
//        const auto loopPatches(periodicGlidePlane->filledPatches(linkShifts));

template <int dim>
const typename DislocationLoopPatches<dim>::LocalPatchPositionType& DislocationLoopPatches<dim>::localPatches() const
{
    return *this;
}

template <int dim>
const typename DislocationLoopPatches<dim>::GlobalPatchPositionType& DislocationLoopPatches<dim>::globalPatches() const
{
    return *this;
}

template <int dim>
const typename DislocationLoopPatches<dim>::VectorDim DislocationLoopPatches<dim>::orientedArea() const
{
    VectorDim nA(VectorDim::Zero());
    for(const auto& patch : globalPatches())
    {
        
        const VectorDim& Ps(patch.second[0]);
        for(size_t k0=0;k0<patch.second.size();++k0)
        {
            const size_t k1(k0+1<patch.second.size()? k0+1 : 0);
            const VectorDim& P0(patch.second[k0]);
            const VectorDim& P1(patch.second[k1]);
            nA+= 0.5*(P0-Ps).cross(P1-P0);
        }
    }
    
    return nA;
}


template struct DislocationLoopPatches<3>;


}
#endif

