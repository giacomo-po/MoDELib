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

namespace model
{
    template <int dim, short unsigned int corder>
    DislocationLoop<dim,corder>::DislocationLoop(LoopNetworkType* const net,
                                     const VectorDim& B,
                                     const std::shared_ptr<GlidePlaneType>& glidePlane_in) :
    /* init */ Loop<DislocationLoop>(net,glidePlane_in->grain.rationalLatticeDirection(B))
    /* init */,glidePlane(glidePlane_in)
    /* init */,periodicGlidePlane(this->network().periodicGlidePlaneFactory? this->network().periodicGlidePlaneFactory->getFromKey(glidePlane->key) : nullptr)
    /* init */,grain(glidePlane->grain)
    /* init */,loopType(DislocationLoopIO<dim>::GLISSILELOOP)
    /* init */,nA(VectorDim::Zero())
    /* init */,nAR(VectorDim::Zero())
    /* init */,_slippedArea(0.0)
    /* init */,_slippedAreaRate(0.0)
    /* init */,_rightHandedUnitNormal(VectorDim::Zero())
    /* init */,_rightHandedNormal(grain)
    /* init */,_slipSystem(nullptr)

    {
        VerboseDislocationLoop(1,"Constructing DislocationLoop "<<this->tag()<<std::endl;);
        assert(this->flow().dot(glidePlane->n)==0);
    }

        /**********************************************************************/
    template <int dim, short unsigned int corder>
    DislocationLoop<dim,corder>::DislocationLoop(LoopNetworkType* const net,
                          const VectorDim& B,
                          const int& grainID,
                          const int& _loopType) :
    /* base init */ Loop<DislocationLoop>(net,net->poly.grain(grainID).rationalLatticeDirection(B))
    /*      init */,glidePlane(nullptr)
    /*      init */,periodicGlidePlane(nullptr)
    /*      init */,grain(net->poly.grain(grainID))
    /*      init */,loopType(_loopType)
    /*      init */,nA(VectorDim::Zero())
    /*      init */,nAR(VectorDim::Zero())
    /*      init */,_slippedArea(0.0)
    /*      init */,_slippedAreaRate(0.0)
    /*      init */,_rightHandedUnitNormal(VectorDim::Zero())
    /*      init */,_rightHandedNormal(grain)
    /*      init */,_slipSystem(nullptr)
    {
        VerboseDislocationLoop(1,"Constructing DislocationLoop "<<this->sID<<" without plane."<<std::endl;);
    }

    template <int dim, short unsigned int corder>
    DislocationLoop<dim,corder>::~DislocationLoop()
    {
        VerboseDislocationLoop(1,"Destroying DislocationLoop "<<this->sID<<std::endl;);

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
                            const double oneA2=sqrt(1.0+DislocationStress<dim>::a2);
                            const double s3(s.dot(e3));
                            const double s3A2=sqrt(std::pow(s3,2)+DislocationStress<dim>::a2);
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
                std::vector<std::pair<VectorDim,VectorDim>> segments;
                for(const auto& loopLink : this->loopLinks())
                {
                    segments.emplace_back(loopLink->source->get_P(),loopLink->sink->get_P());
                }
                temp+=planarSolidAngle(x,glidePlane->P,rightHandedUnitNormal(),segments);
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
    std::tuple<double,double,double> DislocationLoop<dim,corder>::loopLength() const
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
        return std::make_tuple(freeLength,junctionLength,boundaryLength);
    }

    template <int dim, short unsigned int corder>
    std::shared_ptr<SlipSystem> DislocationLoop<dim,corder>::searchSlipSystem() const
    {
        for(const auto& ss : grain.slipSystems())
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
        nAR.setZero();
        if(glidePlane)
        {// remove numerical errors by projecting along glidePlane->unitNormal
            if(this->loopLinks().size())
            {
                const VectorDim P0((*this->loopLinks().begin())->source->get_P());
                for(const auto& loopLink : this->loopLinks())
                {
                    const VectorDim linkChord(loopLink->sink->get_P()-loopLink->source->get_P());
                    nA+= 0.5*(loopLink->source->get_P()-P0).cross(loopLink->sink->get_P()-loopLink->source->get_P());
                    // nAR+=0.5*V0.cross(linkChord)+0.25*(loopLink->sink->networkNode->get_V()-loopLink->source->networkNode->get_V()).cross(linkChord);
                    nAR+=0.5*(loopLink->source->networkNode->get_V()+loopLink->sink->networkNode->get_V()).cross(linkChord);
                }
            }
            _slippedArea=nA.norm();
            _slippedAreaRate=nA.dot(nAR)/_slippedArea;
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
                    isPlanar*=tempPlane.contains(loopLink->source->get_P());
                }
                if(!isPlanar)
                {
                    nA.setZero();
                }
                _slippedArea=nA.norm();
                _slippedAreaRate=0.0; //Slipped area rate will be zero because it is sessile (Discuss this with Dr. Po)
                _rightHandedUnitNormal= _slippedArea>FLT_EPSILON? (nA/_slippedArea).eval() : VectorDim::Zero();
                _rightHandedNormal= grain.reciprocalLatticeDirection(_rightHandedUnitNormal);
                VerboseDislocationLoop(3,"non-glide _rightHandedUnitNormal= "<<_rightHandedUnitNormal.transpose()<<std::endl;);
            }
        }
        updateSlipSystem();
    }

    template <int dim, short unsigned int corder>
    int DislocationLoop<dim,corder>::verboseDislocationLoop=0;
    
    template class DislocationLoop<3,0>;

}
#endif
