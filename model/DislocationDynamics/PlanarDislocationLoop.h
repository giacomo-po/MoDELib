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
#include <GlidePlaneObserver.h>
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
    //                            public PlanarPolygon
    {
        
        constexpr static int dim=TypeTraits<Derived>::dim;
        
        
        Eigen::Matrix<double,dim,1> nA;
        double _slippedArea;
        Eigen::Matrix<double,dim,1> _rightHandedNormal;
        
    public:
        
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        
        
        //        constexpr static int dim=_dim;
        typedef PlanarDislocationLoop<Derived> PlanarDislocationLoopType;
        typedef Loop<Derived> BaseLoopType;
        typedef typename BaseLoopType::LoopLinkType LoopLinkType;
        typedef typename TypeTraits<Derived>::LoopNetworkType LoopNetworkType;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        typedef GlidePlane<dim> GlidePlaneType;
        typedef GlidePlaneObserver<dim> GlidePlaneObserverType;
        typedef Eigen::Matrix<long int,dim+1,1> GlidePlaneKeyType;
        
        
        const Grain<dim>& grain;
        const std::shared_ptr<GlidePlaneType> glidePlane;
        static int verbosePlanarDislocationLoop;
        
        //        const GlidePlaneType& glidePlane;
        //        const bool isGlissile;
        
        
        //        /**********************************************************************/
        //        static bool allowedSlipSystem(const LatticeVector<dim>& b,
        //                                      const ReciprocalLatticeDirection<dim>& n,
        //                                      const Grain<dim>& gr)
        //        {
        //
        //            std::cout<<"PlanarDislocationLoop::FINISH HERE. ANOTHER CONDITION FOR isGlissile SHOULD BE THAT N IS ONE OF THE ALLOWED SLIP PLANES"<<std::endl;
        //
        //            return true;
        //        }
        
        /******************************************************************/
        static void initFromFile(const std::string& fileName)
        {
            verbosePlanarDislocationLoop=TextFileParser("inputFiles/DD.txt").readScalar<int>("verbosePlanarDislocationLoop",true);
        }
        
        
        /**********************************************************************/
        template<typename FLowType>
        PlanarDislocationLoop(LoopNetworkType* const dn,
                              const FLowType& B,
                              const VectorDim& N,
                              const VectorDim& P,
                              const int& grainID) :
        /* base init */ BaseLoopType(dn,B),
        //        /*      init */ PlanarPolygon(fabs(B.dot(N))<FLT_EPSILON? B : N.cross(VectorDim::Random()),N),
        /*      init */ nA(VectorDim::Zero()),
        /*      init */ _slippedArea(0.0),
        /*      init */ _rightHandedNormal(VectorDim::Zero()),
        /*      init */ grain(dn->poly.grain(grainID)),
        /*      init */ glidePlane(dn->sharedGlidePlane(dn->mesh,dn->poly.grain(grainID),P,N))
        //        /*      init */ glidePlane(*_glidePlane->get()),
        //        /*      init */,isGlissile(_isGlissile)
        {
            VerbosePlanarDislocationLoop(1,"Constructing PlanarDislocationLoop "<<this->sID<<std::endl;);
            
            //            glidePlane->addLoop(this);
            //            glidePlane->addParentSharedPtr(&glidePlane,isGlissile,this->sID);
            glidePlane->addParentSharedPtr(&glidePlane);
            
        }
        
        /**********************************************************************/
        template<typename FLowType>
        PlanarDislocationLoop(LoopNetworkType* const dn,
                              const FLowType& B,
                              const int& grainID) :
        /* base init */ BaseLoopType(dn,B),
        //        /*      init */ PlanarPolygon(fabs(B.dot(N))<FLT_EPSILON? B : N.cross(VectorDim::Random()),N),
        /*      init */ nA(VectorDim::Zero()),
        /*      init */ _slippedArea(0.0),
        /*      init */ _rightHandedNormal(VectorDim::Zero()),
        /*      init */ grain(dn->poly.grain(grainID)),
        /*      init */ glidePlane(nullptr)
        //        /*      init */ glidePlane(*_glidePlane->get()),
        //        /*      init */,isGlissile(true)
        {
            VerbosePlanarDislocationLoop(1,"Constructing PlanarDislocationLoop "<<this->sID<<" without plane."<<std::endl;);
            
            //            glidePlane->addLoop(this);
            if(!isVirtualBoundaryLoop())
            {
                glidePlane->addParentSharedPtr(&glidePlane);
            }
            //
            
        }
        
        /**********************************************************************/
        PlanarDislocationLoop(const PlanarDislocationLoop& other) :
        /* base init */ BaseLoopType(other),
        //        /*      init */ PlanarPolygon(other),
        /* init */ nA(other.nA),
        /* init */ _slippedArea(0.0),
        /* init */ _rightHandedNormal(VectorDim::Zero()),
        /* init */ grain(other.grain),
        /* init */ glidePlane(other.glidePlane)
        //        /* init */ glidePlane(*_glidePlane->get()),
        //        /* init */,isGlissile(other.isGlissile)
        {
            VerbosePlanarDislocationLoop(1,"Copy-constructing PlanarDislocationLoop "<<this->sID<<std::endl;);
            
            //            glidePlane->addLoop(this);
            
            if(!isVirtualBoundaryLoop())
            {
                glidePlane->addParentSharedPtr(&glidePlane);
            }
        }
        
        /**********************************************************************/
        ~PlanarDislocationLoop()
        {
            VerbosePlanarDislocationLoop(1,"Destroying PlanarDislocationLoop "<<this->sID<<std::endl;);
            
            //            glidePlane->removeLoop(this);
            if(!isVirtualBoundaryLoop())
            {
                glidePlane->removeParentSharedPtr(&glidePlane);
            }
        }
        
        /**********************************************************************/
        VectorDim diracStringDirection(const VectorDim& x) const
        {
            VectorDim temp(VectorDim::Zero());
            if(glidePlane && _slippedArea>FLT_EPSILON)
            {// a right-handed normal for the loop can be determined
                if((x-glidePlane->P).dot(rightHandedNormal())>=0.0)
                {
                    temp= rightHandedNormal();
                }
                else
                {
                    temp=-rightHandedNormal();
                }
            }
            return temp;
        }
        
        /**********************************************************************/
        double solidAngle(const VectorDim& x) const
        {
            double temp(0.0);
            if(isVirtualBoundaryLoop())
            {
                assert(0 && "FINISH HERE. NEED TO BREAK LOOP IN TWO TRIAGLES");
            }
            else
            {
                const VectorDim s(diracStringDirection(x));
                if(s.squaredNorm())
                {
                    for(const auto& loopLink : this->links())
                    {
                        const VectorDim e1((loopLink.second->source()->get_P()-x).normalized());
                        const VectorDim Y1((loopLink.second->  sink()->get_P()-x).normalized());
                        const VectorDim e3(e1.cross(Y1).normalized());
                        const VectorDim e2(e3.cross(e1));
                        const double ydy(e1.dot(Y1));
                        const double w=sqrt((1.0-ydy)/(1.0+ydy));
                        temp+=2.0*atan(s.dot(e3)*w/(1.0-s.dot(e1)-s.dot(e2)*w));
                    }
                }
                else
                {
                    assert(0 && "FINISH HERE. NEED TO HANDLE SPECIAL CASES OF x on the loop plane. Inside or outside");
                }
            }
            return temp;
        }
        
        /**********************************************************************/
        bool isVirtualBoundaryLoop() const
        {
            return glidePlane.get()==nullptr;
        }
        
        //        /**********************************************************************/
        //        bool isPureVirtualBoundaryLoop() const
        //        {
        //            bool temp(isVirtualBoundaryLoop());
        //            if(temp)
        //            {
        //                for(const auto& loopLink : this->links())
        //                {
        //                    temp*=loopLink->plink->isVirtualBoundarySegment();
        //                    if(!temp)
        //                    {
        //                        break;
        //                    }
        //                }
        //            }
        //
        //            return temp;
        //        }
        
        /**********************************************************************/
        VectorDim burgers() const
        {
            return this->flow().cartesian();
        }
        
        //        /**********************************************************************/
        //        std::tuple<double,double,double> loopLength() const
        //        {
        //            double freeLength=0.0;
        //            double boundaryLength=0.0;
        //            double junctionLength=0.0;
        //
        //            for(const auto& link : this->links())
        //            {
        //                if(link.second->pLink->isBoundarySegment())
        //                {
        //                    boundaryLength+=(link.second->sink()->get_P()-link.second->source()->get_P()).norm();
        //                }
        //                else
        //                {
        //                    if(link.second->pLink->loopLinks().size()==1)
        //                    {
        //                        freeLength+=(link.second->sink()->get_P()-link.second->source()->get_P()).norm();
        //                    }
        //                    else
        //                    {
        //                        junctionLength+=(link.second->sink()->get_P()-link.second->source()->get_P()).norm();
        //                    }
        //                }
        //            }
        //
        //            return std::make_tuple(freeLength,junctionLength,boundaryLength);
        //        }
        
        /**********************************************************************/
        void updateGeometry()
        {
            nA.setZero();
            if(this->links().size())
            {
                const VectorDim P0(this->links().begin()->second->source()->get_P());
                for(const auto& loopLink : this->links())
                {
                    nA+= 0.5*(loopLink.second->source()->get_P()-P0).cross(loopLink.second->sink()->get_P()-loopLink.second->source()->get_P());
                }
            }
            
            _slippedArea=nA.norm();
            _rightHandedNormal= _slippedArea>FLT_EPSILON? (nA/_slippedArea).eval() : VectorDim::Zero();
        }
        
        /**********************************************************************/
        const double& slippedArea() const
        {
            return _slippedArea;
        }
        
        /**********************************************************************/
        const VectorDim& rightHandedNormal() const
        {
            return _rightHandedNormal;
        }
        
        /**********************************************************************/
        MatrixDim plasticDistortion() const
        {
            return -burgers()*nA.transpose()/this->network().mesh.volume();
        }
        
        //        /**********************************************************************/
        //        template <class T>
        //        friend T& operator << (T& os, const PlanarDislocationLoopType& dL)
        //        {
        //            os<< PlanarDislocationLoopIO<dim>(dL);
        //            return os;
        //        }
        
    };
    
    template <typename Derived>
    int PlanarDislocationLoop<Derived>::verbosePlanarDislocationLoop=0;
    
}
#endif

//        /**********************************************************************/
//        MatrixDim plasticDistortion() const
//        {
//            return -this->flow().cartesian()*nA.transpose();
//        }


//        static chordInPlane(LoopLinkType* const pL,)


//        /**********************************************************************/
//        void addLink(LoopLinkType* const pL)
//        {
//            BaseLoopType::addLink(pL);
//
////            const VectorDim chord(pL->source()->get_P()-pL->sink()->get_P());
////            const double chornNorm(chord.norm());
////            if(chornNorm>FLT_EPSILON)
////            {
////            assert(std::fabs((chord/chornNorm).dot(glidePlane->unitNormal))<FLT_EPSILON && "Chord does not belong to plane");
////            }
//
//        }

//        /**********************************************************************/
//        void removeLink(LoopLinkType* const pL)
//        {
//            LoopType::removeLink(pL);
//
//        }

//    template <int _dim, short unsigned int corder, typename InterpolationType>
//    bool PlanarDislocationLoop<_dim,corder,InterpolationType>::outputLoopLength=false;

///**********************************************************************/
//void update()
//{
//    // Collect node positions along loop
//    //            std::deque<VectorDim,Eigen::aligned_allocator<VectorDim>> nodePos;
//    //            std::deque<size_t> nodeIDs;
//    //            for(const auto& nodePair : this->nodeSequence)
//    //            {
//    //                nodePos.push_baack(nodePair.first->get_P());
//    //                nodeIDs.push_back(nodePair.first->sID);
//    //            }
//    //
//    //            // Call PlanarPolygon
//    //            this->assignPoints(nodePos);
//    //            std::deque<std::array<size_t, 3>> triangles=this->triangulate();
//    //
//    //            nA.setZero();
//    //            for(const auto& triIDs : triangles)
//    //            {
//    //                const VectorDim ndA=(nodePos[std::get<0>(triIDs)]-nodePos[std::get<2>(triIDs)]).cross(nodePos[std::get<1>(triIDs)]-nodePos[std::get<0>(triIDs)]);
//    //                assert(nA.dot(ndA)>=0.0 && "Tringles must all be right-handed");
//    //                nA+=ndA;
//    //            }
//    
//
//    
//}