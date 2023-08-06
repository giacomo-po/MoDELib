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
        localPatches().clear();
        globalPatches().clear();
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
                globalLoopPos.push_back(periodicGlidePlane->referencePlane->globalPosition(localPos)+patch->shift);
                localLoopPos.push_back(patch->glidePlane->localPosition(globalLoopPos.back()));
            }
            globalPatches().emplace(patch,globalLoopPos);
            localPatches().emplace(patch,localLoopPos);
            
        }
    }
}
//        const auto loopPatches(periodicGlidePlane->filledPatches(linkShifts));

template <int dim>
const typename DislocationLoopPatches<dim>::LocalPatchPositionType& DislocationLoopPatches<dim>::localPatches() const
{
//    return *this;
    return _localPatches;

}

template <int dim>
const typename DislocationLoopPatches<dim>::GlobalPatchPositionType& DislocationLoopPatches<dim>::globalPatches() const
{
//    return *this;
    return _globalPatches;
}

template <int dim>
typename DislocationLoopPatches<dim>::LocalPatchPositionType& DislocationLoopPatches<dim>::localPatches()
{
//    return *this;
    return _localPatches;

}

template <int dim>
typename DislocationLoopPatches<dim>::GlobalPatchPositionType& DislocationLoopPatches<dim>::globalPatches()
{
//    return *this;
    return _globalPatches;
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

template <int dim>
double DislocationLoopPatches<dim>::solidAngle(const VectorDim& x) const
{
    double temp(0.0);
    const VectorDim rhA(orientedArea());
    const double rhAnorm(rhA.norm());
    if(rhAnorm>FLT_EPSILON)
    {
        const VectorDim rhN(rhA/rhAnorm);
        for(const auto& patch : globalPatches())
        {
            if(patch.second.size())
            {
                const VectorDim& planePoint(patch.second[0]);
                const double posNorm((x-planePoint).norm());
                const double dotProd((x-planePoint).dot(rhN));
                if(std::fabs(dotProd)>FLT_EPSILON*posNorm)
                {// x is outside the plane of the loop
                    const VectorDim s(sgn(dotProd)*rhN); // s points along +n for points above, and along -n for points below
                    for(size_t k=0;k<patch.second.size();++k)
                    {
                        const size_t k1(k<patch.second.size()-1? k+1 : 0);
                        
                        VectorDim e1(patch.second[k]-x);
                        const double e1Norm(e1.norm());
                        if(e1Norm>FLT_EPSILON)
                        {
                            e1/=e1Norm;
                            VectorDim Y1(patch.second[k1]-x);
                            const double Y1norm(Y1.norm());
                            if(Y1norm>FLT_EPSILON)
                            {
                                Y1/=Y1norm;
                                VectorDim e3(e1.cross(Y1));
                                const double e3Norm(e3.norm());
                                if(e3Norm>FLT_EPSILON)
                                {// e1 and Y1 are not align. If they are the projection on the unit sphere is a point and therefore there is no contribution to solid angle
                                    e3/=e3Norm; // normalize e3
                                    const VectorDim e2(e3.cross(e1));
                                    const double ydy(e1.dot(Y1));
                                    const double w=sqrt((1.0-ydy)/(1.0+ydy));
                                    const double oneA2=sqrt(1.0+DislocationFieldBase<dim>::a2);
                                    const double s3(s.dot(e3));
                                    const double s3A2=sqrt(std::pow(s3,2)+DislocationFieldBase<dim>::a2);
                                    temp+=2.0*s3/oneA2/s3A2*atan(s3A2*w/(oneA2-s.dot(e1)-s.dot(e2)*w));
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return temp;
}

template <int dim>
template <typename T>
int DislocationLoopPatches<dim>::sgn(const T& val)
{
    return (val > T(0)) - (val < T(0));
}

template struct DislocationLoopPatches<3>;


}
#endif

