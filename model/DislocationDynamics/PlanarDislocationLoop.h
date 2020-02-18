/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PlanarDislocationLoop_H_
#define model_PlanarDislocationLoop_H_

#include <memory>

#include <Loop.h>
#include <LatticeVector.h>
#include <LatticePlaneBase.h>
#include <LatticePlane.h>
#include <Grain.h>
#include <GlidePlane.h>
#include <SlipSystem.h>
//#include <PeriodicDislocationLoopPair.h>
//#include <BoundaryLoopLinkSequence.h>
#include <PeriodicDislocationLoop.h>
//#include <PlanarDislocationLoopIO.h>
//#include <PlanarPolygon.h>

#ifndef NDEBUG
#define VerbosePlanarDislocationLoop(N,x) if(verbosePlanarDislocationLoop>=N){model::cout<<x;}
#else
#define VerbosePlanarDislocationLoop(N,x)
#endif


namespace model
{
    template <typename Derived>
    class PlanarDislocationLoop : public Loop<Derived>
    {
        
    public:
        
        
        constexpr static int dim=TypeTraits<Derived>::dim;

        typedef PlanarDislocationLoop<Derived> PlanarDislocationLoopType;
        typedef Loop<Derived> BaseLoopType;
        typedef typename BaseLoopType::LoopLinkType LoopLinkType;
        typedef typename TypeTraits<Derived>::LoopNetworkType LoopNetworkType;
        typedef typename TypeTraits<Derived>::LoopType LoopType;
        typedef typename TypeTraits<Derived>::NodeType NodeType;
        
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        typedef GlidePlane<dim> GlidePlaneType;
        typedef Eigen::Matrix<long int,dim+1,1> GlidePlaneKeyType;
        typedef PeriodicDislocationLoop<LoopNetworkType> PeriodicDislocationLoopType;

        const std::shared_ptr<GlidePlaneType> glidePlane;
        const Grain<dim>& grain;
        const int loopType;
        const VectorDim periodicShift;
        const std::shared_ptr<PeriodicDislocationLoopType> periodicLoop;

        static int verbosePlanarDislocationLoop;

        
    private:
        
        Eigen::Matrix<double,dim,1> nA;
        double _slippedArea;
        Eigen::Matrix<double,dim,1> _rightHandedUnitNormal;
        ReciprocalLatticeDirection<dim> _rightHandedNormal;
        std::shared_ptr<SlipSystem> _slipSystem;

        
        template <typename T>
        static int sgn(const T& val)
        {
            return (val > T(0)) - (val < T(0));
        }
        
        std::shared_ptr<PeriodicDislocationLoopType> getPeriodicLoop() const
        {
            if(this->network().simulationParameters.isPeriodicSimulation())
            {
                if(glidePlane)
                {
                    return this->network().periodicDislocationLoopFactory->get(*glidePlane);
                }
                else
                {
                    return nullptr;
                }
            }
            else
            {
                return nullptr;
            }
        }
        
    public:
        
        
        std::shared_ptr<SlipSystem> searchSlipSystem() const
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

        /******************************************************************/
        void updateSlipSystem()
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
                for(const auto& link : this->links())
                {
                    link.second->pLink->updateSlipSystem();
                }
            }
            
            if(_slipSystem)
            {
                    VerbosePlanarDislocationLoop(3,"_slipSystem= "<<_slipSystem->s.cartesian().transpose()<<std::endl;);
                    VerbosePlanarDislocationLoop(3,"_slipSystem= "<<_slipSystem->unitNormal.transpose()<<std::endl;);
            }


        }
        
        /******************************************************************/
        static void initFromFile(const std::string& fileName)
        {
            verbosePlanarDislocationLoop=TextFileParser(fileName).readScalar<int>("verbosePlanarDislocationLoop",true);
        }
        
        
        /**********************************************************************/
        template<typename FLowType>
        PlanarDislocationLoop(LoopNetworkType* const dn,
                              const FLowType& B,
                              const std::shared_ptr<GlidePlaneType>& glidePlane_in) :
        /* base init */ BaseLoopType(dn,B)
        /*      init */,glidePlane(glidePlane_in)
        /*      init */,grain(glidePlane->grain)
        /*      init */,loopType(DislocationLoopIO<dim>::GLISSILELOOP)
        /*      init */,periodicShift(VectorDim::Zero())
        /*      init */,periodicLoop(getPeriodicLoop())
        /*      init */,nA(VectorDim::Zero())
        /*      init */,_slippedArea(0.0)
        /*      init */,_rightHandedUnitNormal(VectorDim::Zero())
        /*      init */,_rightHandedNormal(grain)
        /*      init */,_slipSystem(nullptr)
        {
            VerbosePlanarDislocationLoop(1,"Constructing PlanarDislocationLoop "<<this->sID<<std::endl;);
            
            assert(this->flow().dot(glidePlane->n)==0);
        }
        
        /**********************************************************************/
        template<typename FLowType>
        PlanarDislocationLoop(LoopNetworkType* const dn,
                              const FLowType& B,
                              const std::shared_ptr<GlidePlaneType>& glidePlane_in,
//                              const std::shared_ptr<PeriodicDislocationLoopType>& pLoop_in,
                              const VectorDim& shift_in
                              ) :
        /* base init */ BaseLoopType(dn,B)
        /*      init */,glidePlane(glidePlane_in)
        /*      init */,grain(glidePlane->grain)
        /*      init */,loopType(DislocationLoopIO<dim>::GLISSILELOOP)
        /*      init */,periodicShift(shift_in)
        /*      init */,periodicLoop(getPeriodicLoop())
        /*      init */,nA(VectorDim::Zero())
        /*      init */,_slippedArea(0.0)
        /*      init */,_rightHandedUnitNormal(VectorDim::Zero())
        /*      init */,_rightHandedNormal(grain)
        /*      init */,_slipSystem(nullptr)
        {
            VerbosePlanarDislocationLoop(1,"Constructing PlanarDislocationLoop "<<this->sID<<" in PeriodicLoop "<<periodicLoop->sID<<std::endl;);

            assert(this->flow().dot(glidePlane->n)==0);

//            periodicLoop->addLoop(this->p_derived());
        }
        
        /**********************************************************************/
        template<typename FLowType>
        PlanarDislocationLoop(LoopNetworkType* const dn,
                              const FLowType& B,
                              const int& grainID,
                              const int& _loopType) :
        /* base init */ BaseLoopType(dn,B)
        /*      init */,glidePlane(nullptr)
        /*      init */,grain(dn->poly.grain(grainID))
        /*      init */,loopType(_loopType)
        /*      init */,periodicShift(VectorDim::Zero())
        /*      init */,periodicLoop(nullptr)
        /*      init */,nA(VectorDim::Zero())
        /*      init */,_slippedArea(0.0)
        /*      init */,_rightHandedUnitNormal(VectorDim::Zero())
        /*      init */,_rightHandedNormal(grain)
        /*      init */,_slipSystem(nullptr)
        {
            VerbosePlanarDislocationLoop(1,"Constructing PlanarDislocationLoop "<<this->sID<<" without plane."<<std::endl;);
        }
        
        /**********************************************************************/
        PlanarDislocationLoop(const PlanarDislocationLoop& other) :
        /* base init */ BaseLoopType(other)
        /*      init */,glidePlane(other.glidePlane)
        /*      init */,grain(other.grain)
        /*      init */,loopType(other.loopType)
        /*      init */,periodicShift(other.periodicShift)
        /*      init */,periodicLoop(getPeriodicLoop())
        /*      init */,nA(other.nA)
        /*      init */,_slippedArea(0.0)
        /*      init */,_rightHandedUnitNormal(VectorDim::Zero())
        /*      init */,_rightHandedNormal(grain)
        /*      init */,_slipSystem(nullptr)
        {
            VerbosePlanarDislocationLoop(1,"Copy-constructing PlanarDislocationLoop "<<this->sID<<std::endl;);
        }
        
        /**********************************************************************/
        ~PlanarDislocationLoop()
        {
            VerbosePlanarDislocationLoop(1,"Destroying PlanarDislocationLoop "<<this->sID<<std::endl;);
            
//            if(periodicLoop)
//            {
//                periodicLoop->removeLoop(this->p_derived());
//            }

        }
        
        /**********************************************************************/
        void addLoopLink(LoopLinkType* const pL)
        {
            Loop<Derived>::addLoopLink(pL); // forward to base class
            if(periodicLoop)
            {
                periodicLoop->addLoopLink(pL);
            }
            //            updateBoundaryDecomposition();
        }
        
        /**********************************************************************/
        void removeLoopLink(LoopLinkType* const pL)
        {
            Loop<Derived>::removeLoopLink(pL); // forward to base class
            if(periodicLoop)
            {
                periodicLoop->removeLoopLink(pL);
            }
            //            updateBoundaryDecomposition();
        }
        
        /**********************************************************************/
        static double planarSolidAngle(const VectorDim& x,
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
        
        /**********************************************************************/
        double solidAngle(const VectorDim& x) const
        {
            double temp(0.0);
            if(glidePlane)
            {
                if(_slippedArea>FLT_EPSILON)
                {// a right-handed normal for the loop can be determined
                    
                    std::vector<std::pair<VectorDim,VectorDim>> segments;
                    for(const auto& loopLink : this->links())
                    {
                        segments.emplace_back(loopLink.second->source()->get_P(),loopLink.second->sink()->get_P());
                    }
                    temp+=planarSolidAngle(x,glidePlane->P,rightHandedUnitNormal(),segments);
                    
                    //                    const double posNorm((x-glidePlane->P).norm());
                    //                    const double dotProd((x-glidePlane->P).dot(rightHandedUnitNormal()));
                    //                    if(std::fabs(dotProd)>FLT_EPSILON*posNorm)
                    //                    {// x is outside the plane of the loop
                    //                        const VectorDim s(sgn(dotProd)*rightHandedUnitNormal()); // s points along +n for points above, and along -n for points below
                    //                        for(const auto& loopLink : this->links())
                    //                        {
                    //                            VectorDim e1(loopLink.second->source()->get_P()-x);
                    //                            const double e1Norm(e1.norm());
                    //                            if(e1Norm>FLT_EPSILON)
                    //                            {
                    //                                e1/=e1Norm;
                    //                                VectorDim Y1(loopLink.second->sink()->get_P()-x);
                    //                                const double Y1norm(Y1.norm());
                    //                                if(Y1norm>FLT_EPSILON)
                    //                                {
                    //                                    Y1/=Y1norm;
                    //                                    VectorDim e3(e1.cross(Y1));
                    //                                    const double e3Norm(e3.norm());
                    //                                    if(e3Norm>FLT_EPSILON)
                    //                                    {// e1 and Y1 are not align. If they are the projection on the unit sphere is a point and therefore there is no contribution to solid angle
                    //                                        e3/=e3Norm; // normalize e3
                    //                                        const VectorDim e2(e3.cross(e1));
                    //                                        const double ydy(e1.dot(Y1));
                    //                                        const double w=sqrt((1.0-ydy)/(1.0+ydy));
                    //                                        const double oneA2=sqrt(1.0+DislocationStress<dim>::a2);
                    //                                        const double s3(s.dot(e3));
                    //                                        const double s3A2=sqrt(std::pow(s3,2)+DislocationStress<dim>::a2);
                    //                                        temp+=2.0*s3/oneA2/s3A2*atan(s3A2*w/(oneA2-s.dot(e1)-s.dot(e2)*w));
                    //                                    }
                    //                                }
                    //                            }
                    //                        }
                    //                    }
                }
            }
            else
            {
                const auto linkSeq(this->linkSequence());
                assert(linkSeq.size()==4);
                
                std::vector<std::pair<VectorDim,VectorDim>> triangle0;
                triangle0.emplace_back(linkSeq[0]->source()->get_P(),linkSeq[0]->sink()->get_P());
                triangle0.emplace_back(linkSeq[1]->source()->get_P(),linkSeq[1]->sink()->get_P());
                triangle0.emplace_back(linkSeq[1]->sink()->get_P(),linkSeq[0]->source()->get_P());
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
                triangle1.emplace_back(linkSeq[2]->source()->get_P(),linkSeq[2]->sink()->get_P());
                triangle1.emplace_back(linkSeq[3]->source()->get_P(),linkSeq[3]->sink()->get_P());
                triangle1.emplace_back(linkSeq[3]->sink()->get_P(),linkSeq[2]->source()->get_P());
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
                
        /**********************************************************************/
        bool isVirtualBoundaryLoop() const
        {
            return loopType==DislocationLoopIO<dim>::VIRTUALLOOP;
        }
        
        /**********************************************************************/
        bool isPureVirtualBoundaryLoop() const
        {
            bool temp(isVirtualBoundaryLoop());
            if(temp)
            {
                for(const auto& link : this->links())
                {
                    temp*=link.second->pLink->isVirtualBoundarySegment();
                    if(!temp)
                    {
                        break;
                    }
                }
            }
            return temp;
        }
        
        /**********************************************************************/
        VectorDim burgers() const
        {
            return this->flow().cartesian();
        }

        /**********************************************************************/
        const std::shared_ptr<SlipSystem>&  slipSystem() const
        {
            return _slipSystem;
        }

        
        /**********************************************************************/
        void updateGeometry()
        {
            nA.setZero();
            VerbosePlanarDislocationLoop(3,"PlanarDislocationLoop "<<this->sID<<std::endl;);

            if(glidePlane)
            {// remove numerical errors by projecting along glidePlane->unitNormal
                
                if(this->links().size())
                {
                    const VectorDim P0(this->links().begin()->second->source()->get_P());
                    for(const auto& loopLink : this->links())
                    {
                        nA+= 0.5*(loopLink.second->source()->get_P()-P0).cross(loopLink.second->sink()->get_P()-loopLink.second->source()->get_P());
                    }
                }
                
                
                //                nA=nnDot*glidePlane->unitNormal;
                //                nA=nA.dot(glidePlane->n.cartesian().normalized())*glidePlane->n.cartesian().normalized();
                
                _slippedArea=nA.norm();
                _rightHandedUnitNormal= _slippedArea>FLT_EPSILON? (nA/_slippedArea).eval() : VectorDim::Zero();
                
                const double nnDot(_rightHandedUnitNormal.dot(glidePlane->unitNormal));
                _rightHandedNormal= nnDot>=0.0? glidePlane->n : ReciprocalLatticeDirection<dim>(glidePlane->n*(-1));
                _rightHandedUnitNormal=_rightHandedNormal.cartesian().normalized();
                VerbosePlanarDislocationLoop(3,"_rightHandedUnitNormal= "<<_rightHandedUnitNormal.transpose()<<std::endl;);

                updateSlipSystem();
            }
        }
        
        /**********************************************************************/
        const double& slippedArea() const
        {
            return _slippedArea;
        }
        
        /**********************************************************************/
        const VectorDim& rightHandedUnitNormal() const
        {
            return _rightHandedUnitNormal;
        }
        
        /**********************************************************************/
        const ReciprocalLatticeDirection<dim>& rightHandedNormal() const
        {
            return _rightHandedNormal;
        }
        
        /**********************************************************************/
        MatrixDim plasticDistortion() const
        {
            return -burgers()*nA.transpose()/this->network().mesh.volume();
        }
    };
    
    template <typename Derived>
    int PlanarDislocationLoop<Derived>::verbosePlanarDislocationLoop=0;
    
}
#endif
