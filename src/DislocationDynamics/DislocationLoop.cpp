/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationLoop_cpp_
#define model_DislocationLoop_cpp_


#include <DislocationLoop.h>
#include <DislocationFieldBase.h>
#include <SutherlandHodgman.h>
#include <DislocationLoopPatches.h>

namespace model
{
template <int dim, short unsigned int corder>
DislocationLoop<dim,corder>::DislocationLoop(LoopNetworkType* const net,
                                             const VectorDim& B,
                                             const std::shared_ptr<GlidePlaneType>& glidePlane_in) :
/* init */ Loop<DislocationLoop>(net,glidePlane_in->grain.singleCrystal->rationalLatticeDirection(B))
/* init */,glidePlane(glidePlane_in)
/* init */,periodicGlidePlane(this->network().periodicGlidePlaneFactory? this->network().periodicGlidePlaneFactory->getFromKey(glidePlane->key) : nullptr)
/* init */,_patches(periodicGlidePlane)
/* init */,grain(glidePlane->grain)
/* init */,loopType(this->flow().dot(glidePlane->n)==0? DislocationLoopIO<dim>::GLISSILELOOP : DislocationLoopIO<dim>::SESSILELOOP)
/* init */,nA(VectorDim::Zero())
/* init */,nAR(VectorDim::Zero())
/* init */,_slippedArea(0.0)
/* init */,_slippedAreaRate(0.0)
/* init */,_rightHandedUnitNormal(VectorDim::Zero())
/* init */,_rightHandedNormal(*grain.singleCrystal)
/* init */,_slipSystem(nullptr)

{
    VerboseDislocationLoop(1,"Constructing DislocationLoop "<<this->tag()<<std::endl;);
    //        assert(this->flow().dot(glidePlane->n)==0);
}

template <int dim, short unsigned int corder>
void DislocationLoop<dim,corder>::computeStackingFaultForces()
{
    if(this->slipSystem() && this->glidePlane)
    {
        const double eps=1.0e-2;

        for(const auto& loopLink : this->loopLinks())
        {
            if(loopLink->networkLink())
            {
                
                VectorDim outDir((loopLink->sink->get_P() - loopLink->source->get_P()).cross(this->rightHandedUnitNormal()));
                const double outDirNorm(outDir.norm());
                if(outDirNorm>FLT_EPSILON)
                {
                    outDir/=outDirNorm;
                    std::vector<std::pair<VectorDim,VectorDim>> qPointSlip(loopLink->networkLink()->quadraturePoints().size(),std::make_pair(VectorDim::Zero(),VectorDim::Zero())); // accumulated slip vectors (outside,inside) for each qPoint
                    
                    for(const auto& weakSourceLoop : this->network().loops())
                    {
                        const auto sourceLoop(weakSourceLoop.second.lock());
                        if(sourceLoop->slipSystem())
                        {
                            if(this->slipSystem()->n==sourceLoop->slipSystem()->n)
                            {// same glide plane family
                                
                                for(size_t q=0;q<loopLink->networkLink()->quadraturePoints().size();++q)
                                {
                                    const auto& qPoint(loopLink->networkLink()->quadraturePoints()[q]);
                                    qPointSlip[q].first -=sourceLoop->windingNumber(qPoint.r + eps*outDir)*sourceLoop->burgers(); // slip vector is negative the burgers vector
                                    qPointSlip[q].second-=sourceLoop->windingNumber(qPoint.r - eps*outDir)*sourceLoop->burgers(); // slip vector is negative the burgers vector

                                }
                            }
                        }
                    }
                    
                    for(size_t q=0;q<loopLink->networkLink()->quadraturePoints().size();++q)
                    {
                        if((qPointSlip[q].first-qPointSlip[q].second).squaredNorm()>FLT_EPSILON)
                        {
                            auto& qPoint(loopLink->networkLink()->quadraturePoints()[q]);
                            if(qPoint.inclusionID<0)
                            {// qPoint is not inside an inclusion, we use the matrix gamma-surface

                                const double gamma1(this->slipSystem()->misfitEnergy(qPointSlip[q].first));  // outer point
                                const double gamma2(this->slipSystem()->misfitEnergy(qPointSlip[q].second)); // inner point

                                double gammaNoise(0.0);
                                if(this->slipSystem()->planeNoise)
                                {
                                    if(loopLink->networkLink()->glidePlanes().size()==1)
                                    {
                                        const auto& glidePlane(**loopLink->networkLink()->glidePlanes().begin());
                                        gammaNoise=std::get<2>(this->slipSystem()->gridInterp(qPoint.r-glidePlane.P));
                                    }
                                }
                                qPoint.stackingFaultForce+= -(gamma2-gamma1+gammaNoise)*outDir;
                            }
                            else
                            {// qPoint is inside an inclusion, we use the inclusion gamma-surface

                                const auto& secondPhase(this->network().eshelbyInclusions().at(qPoint.inclusionID)->secondPhase);
                                if(secondPhase)
                                {
                                        const double gamma1(secondPhase->misfitEnergy(qPointSlip[q].first ,this->slipSystem()));  // outer point
                                        const double gamma2(secondPhase->misfitEnergy(qPointSlip[q].second,this->slipSystem())); // inner point
                                    qPoint.stackingFaultForce+= -(gamma2-gamma1)*outDir;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}


template <int dim, short unsigned int corder>
void DislocationLoop<dim,corder>::crossSlipBranches(std::deque<std::pair<std::deque<std::shared_ptr<LoopNodeType>>,int>>& csNodes) const
{
    
    std::deque<std::deque<std::pair<const LoopLinkType*,int>>> csBranches;
    //std::cout<<"Loop "<<this->tag()<<std::endl;
    if(this->network().crossSlipModel)
    {
        std::deque<std::pair<const LoopLinkType*,int>> currentBranch; // pair<link,cross-slip slipSystem ID>
        for(const auto& link : this->linkSequence())
        {
            //std::cout<<"Link "<<link->tag()<<std::endl;
            
            
            if(link->hasNetworkLink())
            {
                const auto isCSLink(this->network().crossSlipModel->isCrossSlipLink(*link->networkLink()));
                //std::cout<<isCSLink.first<<", "<<isCSLink.second.first<<", "<<isCSLink.second.second<<std::endl;
                
                if(isCSLink.first)
                {// a cross-slip segment
                    if(currentBranch.empty())
                    {// start of new branch
                        currentBranch.emplace_back(link,isCSLink.second.second);
                        //std::cout<<"A"<<std::endl;
                        
                    }
                    else
                    {// existing branch
                        if(currentBranch.back().second==isCSLink.second.second)
                        {
                            currentBranch.emplace_back(link,isCSLink.second.second);
                            //std::cout<<"B"<<std::endl;
                            
                        }
                        else
                        {// close, store, and push to currentBranch,
                            csBranches.push_back(currentBranch);
                            currentBranch.clear();
                            currentBranch.emplace_back(link,isCSLink.second.second);
                            //std::cout<<"C"<<std::endl;
                            
                        }
                    }
                }
                else
                {// not a cross-slip segment
                    if(!currentBranch.empty())
                    {// close and store branch if not empty
                        csBranches.push_back(currentBranch);
                        currentBranch.clear();
                        //std::cout<<"D"<<std::endl;
                    }
                }
            }
            else
            {
                if(!currentBranch.empty())
                {// close and store branch if not empty
                    currentBranch.emplace_back(link,currentBranch.back().second);
                    //std::cout<<"E"<<std::endl;
                }
            }
            
        }
        if(!currentBranch.empty())
        {// close and store branch if not empty
            csBranches.push_back(currentBranch);
            currentBranch.clear();
            //std::cout<<"F"<<std::endl;
        }
        
        //std::cout<<"Loop "<<this->tag()<<", csBranches.size()="<<csBranches.size()<<std::endl;
        if(csBranches.size()>1)
        {// Inserted two or more branches. Merge last and first branch if possible
            if(   csBranches[0].front().first->prev==csBranches.back().back().first
               && csBranches[0].front().second==csBranches.back().back().second)
            {
                for(typename std::deque<std::pair<const LoopLinkType*,int>>::reverse_iterator rIter = csBranches.back().rbegin();
                    rIter != csBranches.back().rend(); ++rIter)
                {
                    csBranches[0].push_front(*rIter);
                }
                
                //                    for(const auto& pair : csBranches.front())
                //                    {
                //                        csBranches.back().push_back(pair);
                //                    }
                csBranches.pop_back();
            }
        }
        
        //std::cout<<"Adding csNodes"<<std::endl;
        for(const auto& branch : csBranches)
        {
            if(branch.size())
            {
                //std::cout<<"A"<<std::endl;
                
                csNodes.emplace_back(std::deque<std::shared_ptr<LoopNodeType>>(),branch.back().second);
                for(const auto& pair : branch)
                {
                    //std::cout<<"B"<<std::endl;
                    
                    if(!pair.first->source->periodicPlaneEdge.first && !pair.first->source->periodicPlaneEdge.second)
                    {// not a boundary node
                        //std::cout<<"C"<<std::endl;
                        
                        csNodes.back().first.emplace_back(pair.first->source);
                        //std::cout<<"D"<<std::endl;
                        
                    }
                }
                //                    //std::cout<<"E"<<std::endl;
                //                    auto e1(branch.back().first->sink);
                //                    auto e2(csNodes.back().first.back());
                //                    bool aa(branch.back().first->sink!=csNodes.back().first.back());
                //                    //std::cout<<"E1 "<<aa<<std::endl;
                //                    bool bb(branch.back().first->sink->periodicPlaneEdge.first);
                //                    //std::cout<<"E2 "<<bb<<std::endl;
                //                    bool cc(branch.back().first->sink->periodicPlaneEdge.second);
                //                    //std::cout<<"E3 "<<cc<<std::endl;
                
                if(csNodes.back().first.size())
                {
                    if(    branch.back().first->sink!=csNodes.back().first.back()
                       && !branch.back().first->sink->periodicPlaneEdge.first && !branch.back().first->sink->periodicPlaneEdge.second)
                    {
                        //std::cout<<"F"<<std::endl;
                        
                        csNodes.back().first.emplace_back(branch.back().first->sink);
                        //std::cout<<"G"<<std::endl;
                        
                    }
                }
            }
        }
        
        //            //std::cout<<"Loop "<<this->tag()<<", csBranches.size()="<<csBranches.size()<<std::endl;
        
        
        //            //std::cout<<"Loop "<<this->tag()<<": csBranches.size()="<<csBranches.size()<<std::endl;
        //            for(const auto& brach : csBranches)
        //            {
        //                for(const auto& pair : brach)
        //                {
        //                    //std::cout<<pair.first->tag()<<std::endl;
        //                }
        //            }
        
    }
    
    
    
    
    
}

template <int dim, short unsigned int corder>
int DislocationLoop<dim,corder>::windingNumber(const VectorDim& pt)
{/*!\param[in] pts,vector of positions about which to compute the winding number of this loop
  
  */
    int wn(0);
    for(const auto& pair : _patches.localPatches())
    {
        if(pair.first->glidePlane)
        {
            if(pair.first->glidePlane->contains(pt))
            {
                const auto localPt(pair.first->glidePlane->localPosition(pt));
                wn+=Polygon2D::windingNumber(localPt,pair.second);
            }
        }
    }
    return wn;
}

template <int dim, short unsigned int corder>
int DislocationLoop<dim,corder>::windingNumber(const Eigen::Matrix<double,dim-1,1>& localPt,const std::shared_ptr<GlidePlane<dim>>& ptPlane)
{/*!\param[in] pts,vector of positions about which to compute the winding number of this loop
  
  */
    int wn(0);
    for(const auto& pair : _patches.localPatches())
    {
        if(pair.first->glidePlane)
        {
            if(pair.first->glidePlane==ptPlane)
            {
                //const auto localPt(pair.first->glidePlane->localPosition(pt));
                wn+=Polygon2D::windingNumber(localPt,pair.second);
            }
        }
    }
    return wn;
}



template <int dim, short unsigned int corder>
DislocationLoop<dim,corder>::~DislocationLoop()
{
    VerboseDislocationLoop(1,"Destroying DislocationLoop "<<this->sID<<std::endl;);
    
}

template <int dim, short unsigned int corder>
const DislocationLoopPatches<dim>& DislocationLoop<dim,corder>::patches() const
{
    return _patches;
}

template <int dim, short unsigned int corder>
std::shared_ptr<typename TypeTraits<DislocationLoop<dim,corder>>::LoopType> DislocationLoop<dim,corder>::clone() const
{
    return std::shared_ptr<typename TypeTraits<DislocationLoop<dim,corder>>::LoopType>(new DislocationLoop(this->p_network(),burgers(),glidePlane));
}

template <int dim, short unsigned int corder>
void DislocationLoop<dim,corder>::initFromFile(const std::string& fileName)
{
    verboseDislocationLoop=TextFileParser(fileName).readScalar<int>("verboseDislocationLoop",true);
}

template <int dim, short unsigned int corder>
const double& DislocationLoop<dim,corder>::slippedArea() const
{
    return _slippedArea;
}

template <int dim, short unsigned int corder>
const double& DislocationLoop<dim,corder>::slippedAreaRate() const
{
    return _slippedAreaRate;
}

template <int dim, short unsigned int corder>
const typename DislocationLoop<dim,corder>::VectorDim& DislocationLoop<dim,corder>::rightHandedUnitNormal() const
{
    return _rightHandedUnitNormal;
}

template <int dim, short unsigned int corder>
const typename DislocationLoop<dim,corder>::ReciprocalLatticeDirectionType& DislocationLoop<dim,corder>::rightHandedNormal() const
{
    return _rightHandedNormal;
}

template <int dim, short unsigned int corder>
bool DislocationLoop<dim,corder>::isVirtualBoundaryLoop() const
{
    return loopType==DislocationLoopIO<dim>::VIRTUALLOOP;
}

template <int dim, short unsigned int corder>
double DislocationLoop<dim,corder>::planarSolidAngle(const VectorDim& x,
                                                     const VectorDim& planePoint,
                                                     const VectorDim& rhN,
                                                     const std::vector<std::pair<VectorDim,VectorDim>>& polygonSegments)
{
    double temp(0.0);
    const double posNorm((x-planePoint).norm());
    const double dotProd((x-planePoint).dot(rhN));
    if(std::fabs(dotProd)>FLT_EPSILON*posNorm)
    {// x is outside the plane of the loop
        const VectorDim s(sgn(dotProd)*rhN); // s points along +n for points above, and along -n for points below
        for(const auto& pair : polygonSegments)
        {
            VectorDim e1(pair.first-x);
            const double e1Norm(e1.norm());
            if(e1Norm>FLT_EPSILON)
            {
                e1/=e1Norm;
                VectorDim Y1(pair.second-x);
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
    return temp;
}

template <int dim, short unsigned int corder>
template <typename T>
int DislocationLoop<dim,corder>::sgn(const T& val)
{
    return (val > T(0)) - (val < T(0));
}

template <int dim, short unsigned int corder>
double DislocationLoop<dim,corder>::solidAngle(const VectorDim& x) const
{
    double temp(0.0);
    if(glidePlane)
    {
        if(_slippedArea>FLT_EPSILON)
        {// a right-handed normal for the loop can be determined
            temp+=_patches.solidAngle(x);
//            for(const auto& pair : _patches.globalPatches())
//            {
////                const auto patchGlidePlane(pair.first->patchBoundary->referencePlane);
////                std::vector<std::pair<VectorDim,VectorDim>> segments;
////                for(size_t k=0;k<pair.second.size();++k)
////                {
////                    const size_t k1(k+1==pair.second.size()? 0 : k+1);
//////                    segments.emplace_back(patchGlidePlane->globalPosition(pair.second[k]),patchGlidePlane->globalPosition(pair.second[k1]));
////                    segments.emplace_back(pair.second[k],pair.second[k1]);
////                }
////                temp+=planarSolidAngle(x,patchGlidePlane->P,rightHandedUnitNormal(),segments);
//                temp+=solidAngle(x);
//            }
            
            
            //                std::vector<std::pair<VectorDim,VectorDim>> segments;
            //                for(const auto& loopLink : this->loopLinks())
            //                {
            //                    segments.emplace_back(loopLink->source->get_P(),loopLink->sink->get_P());
            //                }
            //                temp+=planarSolidAngle(x,glidePlane->P,rightHandedUnitNormal(),segments);
        }
    }
    else
    {
        const auto linkSeq(this->linkSequence());
        assert(linkSeq.size()==4);
        
        std::vector<std::pair<VectorDim,VectorDim>> triangle0;
        triangle0.emplace_back(linkSeq[0]->source->get_P(),linkSeq[0]->sink->get_P());
        triangle0.emplace_back(linkSeq[1]->source->get_P(),linkSeq[1]->sink->get_P());
        triangle0.emplace_back(linkSeq[1]->sink->get_P(),linkSeq[0]->source->get_P());
        const VectorDim planePoint0(triangle0[0].first);
        VectorDim rhN0(VectorDim::Zero());
        for(const auto& pair : triangle0)
        {
            rhN0+= 0.5*(pair.first-planePoint0).cross(pair.second-pair.first);
        }
        const double rhN0norm(rhN0.norm());
        if(rhN0norm>FLT_EPSILON)
        {
            temp+=planarSolidAngle(x,planePoint0,rhN0/rhN0norm,triangle0);
        }
        
        std::vector<std::pair<VectorDim,VectorDim>> triangle1;
        triangle1.emplace_back(linkSeq[2]->source->get_P(),linkSeq[2]->sink->get_P());
        triangle1.emplace_back(linkSeq[3]->source->get_P(),linkSeq[3]->sink->get_P());
        triangle1.emplace_back(linkSeq[3]->sink->get_P(),linkSeq[2]->source->get_P());
        const VectorDim planePoint1(triangle1[0].first);
        VectorDim rhN1(VectorDim::Zero());
        for(const auto& pair : triangle1)
        {
            rhN1+= 0.5*(pair.first-planePoint1).cross(pair.second-pair.first);
        }
        const double rhN1norm(rhN1.norm());
        if(rhN1norm>FLT_EPSILON)
        {
            temp+=planarSolidAngle(x,planePoint1,rhN1/rhN1norm,triangle1);
        }
    }
    return temp;
}

template <int dim, short unsigned int corder>
std::tuple<double,double,double,double> DislocationLoop<dim,corder>::loopLength() const
{
    double freeLength=0.0;
    double boundaryLength=0.0;
    double junctionLength=0.0;
    for(const auto& link : this->loopLinks())
    {
        if(link->networkLink())
        {
            if(link->networkLink()->isBoundarySegment())
            {
                boundaryLength+=(link->sink->get_P()-link->source->get_P()).norm();
            }
            else
            {
                if(link->networkLink()->loopLinks().size()==1)
                {
                    freeLength+=(link->sink->get_P()-link->source->get_P()).norm();
                }
                else
                {
                    junctionLength+=(link->sink->get_P()-link->source->get_P()).norm();
                }
            }
        }
    }
    return std::make_tuple(freeLength,junctionLength,boundaryLength,slippedArea());
}

template <int dim, short unsigned int corder>
std::shared_ptr<SlipSystem> DislocationLoop<dim,corder>::searchSlipSystem() const
{
    for(const auto& ss : grain.singleCrystal->slipSystems())
    {
        if(ss->isSameAs(this->flow(),_rightHandedNormal))
        {
            return ss;
        }
    }
    return std::shared_ptr<SlipSystem>(nullptr);
}

template <int dim, short unsigned int corder>
typename DislocationLoop<dim,corder>::MatrixDim DislocationLoop<dim,corder>::plasticDistortion() const
{
    return -burgers()*nA.transpose()/this->network().mesh.volume();
}

//Added by Yash

template <int dim, short unsigned int corder>
typename DislocationLoop<dim,corder>::MatrixDim DislocationLoop<dim,corder>::plasticDistortionRate() const
{
    return -burgers()*nAR.transpose()/this->network().mesh.volume();
}

template <int dim, short unsigned int corder>
typename DislocationLoop<dim,corder>::VectorDim DislocationLoop<dim,corder>::burgers() const
{
    return this->flow().cartesian();
}

template <int dim, short unsigned int corder>
void DislocationLoop<dim,corder>::updateSlipSystem()
{
    VerboseDislocationLoop(3,"DislocationLoop "<<this->tag()<<std::endl;);
    
    const SlipSystem* const oldSlipSystem(_slipSystem.get());
    if(glidePlane)
    {// a glide plane exists
        if(_slipSystem)
        {// a current slip system exists
            if(_slipSystem->isSameAs(this->flow(),_rightHandedNormal))
            {// currenst slip system still valid, don't do anything
            }
            else
            {// currenst slip system not valid
                _slipSystem=searchSlipSystem();
            }
        }
        else
        {// a current slip system does not exist
            _slipSystem=searchSlipSystem();
        }
    }
    else
    {// no glide plane
        _slipSystem=nullptr;
    }
    
    if(_slipSystem.get()!=oldSlipSystem)
    {
        for(const auto& link : this->loopLinks())
        {
            if(link->networkLink())
            {
                link->networkLink()->updateSlipSystem();
            }
        }
    }
    if(_slipSystem)
    {
        VerboseDislocationLoop(3,"_slipSystem.s= "<<_slipSystem->s.cartesian().transpose()<<std::endl;);
        VerboseDislocationLoop(3,"_slipSystem.unitNormal= "<<_slipSystem->unitNormal.transpose()<<std::endl;);
    }
    else
    {
        VerboseDislocationLoop(3,"slipSystem NOT FOUND"<<std::endl;);
    }
}

template <int dim, short unsigned int corder>
const std::shared_ptr<SlipSystem>&  DislocationLoop<dim,corder>::slipSystem() const
{
    return _slipSystem;
}

template <int dim, short unsigned int corder>
void DislocationLoop<dim,corder>::updateGeometry()
{
    VerboseDislocationLoop(2,"DislocationLoop "<<this->sID<<" updateGeometry"<<std::endl;);
    nA.setZero();
//    nAR.setZero();
    if(glidePlane)
    {// remove numerical errors by projecting along glidePlane->unitNormal
        if(this->loopLinks().size())
        {
            const VectorDim P0((*this->loopLinks().begin())->source->get_P());
            for(const auto& loopLink : this->loopLinks())
            {
//                const VectorDim linkChord(loopLink->sink->get_P()-loopLink->source->get_P());
                nA+= 0.5*(loopLink->source->get_P()-P0).cross(loopLink->sink->get_P()-loopLink->source->get_P());
//                nAR+=0.5*(loopLink->source->networkNode->get_V()+loopLink->sink->networkNode->get_V()).cross(linkChord);
            }
        }
        _slippedArea=nA.norm();
//        _slippedAreaRate=nA.dot(nAR)/_slippedArea;
        _rightHandedUnitNormal= _slippedArea>FLT_EPSILON? (nA/_slippedArea).eval() : VectorDim::Zero();
        const double nnDot(_rightHandedUnitNormal.dot(glidePlane->unitNormal));
        _rightHandedNormal= nnDot>=0.0? glidePlane->n : ReciprocalLatticeDirection<dim>(glidePlane->n*(-1));
        _rightHandedUnitNormal=_rightHandedNormal.cartesian().normalized();
        VerboseDislocationLoop(3,"_rightHandedUnitNormal= "<<_rightHandedUnitNormal.transpose()<<std::endl;);
        //            updateSlipSystem();
    }
    else
    {
        if(this->loopLinks().size())
        {
            const VectorDim P0((*this->loopLinks().begin())->source->get_P());
            for(const auto& loopLink : this->loopLinks())
            {
                nA+= 0.5*(loopLink->source->get_P()-P0).cross(loopLink->sink->get_P()-loopLink->source->get_P());
            }
            Plane<dim> tempPlane(P0,nA);
            bool isPlanar(true);
            for(const auto& loopLink : this->loopLinks())
            {
                isPlanar= (isPlanar && tempPlane.contains(loopLink->source->get_P()));
            }
            if(!isPlanar)
            {
                nA.setZero();
            }
            _slippedArea=nA.norm();
//            _slippedAreaRate=0.0; //Slipped area rate will be zero because it is sessile (Discuss this with Dr. Po)
            _rightHandedUnitNormal= _slippedArea>FLT_EPSILON? (nA/_slippedArea).eval() : VectorDim::Zero();
            _rightHandedNormal= grain.singleCrystal->reciprocalLatticeDirection(_rightHandedUnitNormal);
            VerboseDislocationLoop(3,"non-glide _rightHandedUnitNormal= "<<_rightHandedUnitNormal.transpose()<<std::endl;);
        }
    }
    updateSlipSystem();
  
    if(periodicGlidePlane)
    {
        std::vector<VectorDim> linkShifts;
        std::vector<VectorDim> globalNodePos;

        const auto linkSeq(this->linkSequence());
        for(const auto& link : linkSeq)
        {// Collect patches of the loop
            if(link->periodicPlanePatch())
            {
                linkShifts.push_back(link->periodicPlanePatch()->shift);
            }
            if(   !link->source->periodicPlaneEdge.first
               && !link->source->periodicPlaneEdge.second)
            {// source is not a bnd node
                globalNodePos.push_back(link->source->get_P());
            }
        }
        _patches.update(linkShifts,globalNodePos);
    }
}

template <int dim, short unsigned int corder>
void DislocationLoop<dim,corder>::updateRates()
{
    VerboseDislocationLoop(2,"DislocationLoop "<<this->sID<<" updateRates"<<std::endl;);
    nAR.setZero();
    if(glidePlane)
    {// remove numerical errors by projecting along glidePlane->unitNormal
        if(this->loopLinks().size())
        {
            for(const auto& loopLink : this->loopLinks())
            {
                const VectorDim linkChord(loopLink->sink->get_P()-loopLink->source->get_P());
                nAR+=0.5*(loopLink->source->networkNode->get_V()+loopLink->sink->networkNode->get_V()).cross(linkChord);
            }
        }
        _slippedAreaRate=nA.dot(nAR)/_slippedArea;
    }
}

//template <int dim, short unsigned int corder>
//const std::vector<typename DislocationLoop<dim,corder>::VectorDim>& DislocationLoop<dim,corder>::barycentricNodes() const
//{
//    return _barycentricNodes;
//}

template <int dim, short unsigned int corder>
int DislocationLoop<dim,corder>::verboseDislocationLoop=0;

template class DislocationLoop<3,0>;

}
#endif
