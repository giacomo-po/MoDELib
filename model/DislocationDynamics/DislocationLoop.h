/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po         <gpo@ucla.edu>.
 * Copyright (C) 2019 by Yash Pachaury      <ypachaur@purdue.edu>
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationLoop_h_
#define model_DislocationLoop_h_

#include <memory>

#include <TypeTraits.h>
#include <Loop.h>
//#include <LatticeVector.h>
//#include <LatticePlaneBase.h>
//#include <LatticePlane.h>
//#include <Grain.h>
//#include <GlidePlane.h>
////#include <DislocationLoopIO.h>
////#include <PlanarDislocationLoop.h>
//#include <SlipSystem.h>

//#include <PlanarPolygon.h>
//#include <DislocationNetwork.h>

#ifndef NDEBUG
#define VerboseDislocationLoop(N,x) if(verboseDislocationLoop>=N){model::cout<<x;}
#else
#define VerboseDislocationLoop(N,x)
#endif

namespace model
{
    template <int _dim, short unsigned int corder, typename InterpolationType>
    class DislocationLoop : public Loop<DislocationLoop<_dim,corder,InterpolationType>>
    {

    public:
        
        typedef TypeTraits<DislocationLoop<_dim,corder,InterpolationType>> TraitsType;
        typedef typename TraitsType::LoopNetworkType LoopNetworkType;
        typedef typename TraitsType::LoopType LoopType;
        typedef typename TraitsType::LoopNodeType LoopNodeType;
        typedef typename TraitsType::LoopLinkType LoopLinkType;
        typedef typename TraitsType::NetworkNodeType NetworkNodeType;
        typedef typename TraitsType::NetworkLinkType NetworkLinkType;
        typedef typename TraitsType::FlowType FlowType;
        
        static int verboseDislocationLoop;

        
        DislocationLoop(LoopNetworkType* const,
                        const FlowType&);
        
        std::shared_ptr<LoopType> clone() const;
        
        static void initFromFile(const std::string&);
        
//        typedef DislocationLoop<_dim,corder,InterpolationType> LoopType;

//
//        constexpr static int dim=_dim;
//
////        std::shared_ptr<SlipSystem> slipSystem;
////
////        /**********************************************************************/
////        void updateSlipSystem()
////        {
////            slipSystem=nullptr;
//////            const ReciprocalLatticeDirection<dim> rightHandedDir(this->grain.reciprocalLatticeDirection(this->rightHandedNormal()));
////            for(const auto& ss : this->grain.slipSystems())
////            {
////                if(  ((this->flow()-ss->s).squaredNorm()==0 && (this->rightHandedNormal()-ss->n).squaredNorm()==0)
////                   ||((this->flow()+ss->s).squaredNorm()==0 && (this->rightHandedNormal()+ss->n).squaredNorm()==0))
////                {
////                    slipSystem=ss;
////                }
////            }
////        }
//
//    public:
//
//        typedef DislocationLoop<dim,corder,InterpolationType> DislocationLoopType;
//        typedef PlanarDislocationLoop<DislocationLoopType> BaseLoopType;
//        typedef typename TypeTraits<DislocationLoopType>::LoopNetworkType LoopNetworkType;
//        typedef PeriodicDislocationLoop<LoopNetworkType> PeriodicDislocationLoopType;
//        typedef Eigen::Matrix<double,dim,1> VectorDim;
//        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
//
//
//
//        /**********************************************************************/
//        DislocationLoop(LoopNetworkType* const dn,
//                        const VectorDim& B,
//                        const std::shared_ptr<GlidePlane<dim>>& glidePlane) :
////        /* init */ BaseLoopType(dn,glidePlane->grain.latticeVector(B),glidePlane)
//        /* init */ BaseLoopType(dn,glidePlane->grain.rationalLatticeDirection(B),glidePlane)
//        //        /* init */,slipSystem(nullptr)
//        //        /* init */,isGlissile(this->flow().dot(this->glidePlane->n)==0)
//        {
//        }
//
//        /**********************************************************************/
//        DislocationLoop(LoopNetworkType* const dn,
//                        const VectorDim& B,
//                        const std::shared_ptr<GlidePlane<dim>>& glidePlane,
////                        const std::shared_ptr<PeriodicDislocationLoopType>& pLoop,
//                        const VectorDim& shift) :
////        /* init */ BaseLoopType(dn,glidePlane->grain.latticeVector(B),glidePlane,pLoop,shift)
//        /* init */ BaseLoopType(dn,glidePlane->grain.rationalLatticeDirection(B),glidePlane,shift)
//        //        /* init */,slipSystem(nullptr)
//        //        /* init */,isGlissile(this->flow().dot(this->glidePlane->n)==0)
//        {
//        }
//
//        /**********************************************************************/
//        DislocationLoop(LoopNetworkType* const dn,
//                        const VectorDim& B,
//                        const int& grainID,
//                        const int& _loopType) :
//        /* init */ BaseLoopType(dn,dn->poly.grain(grainID).rationalLatticeDirection(B),grainID,_loopType)
////        /* init */,slipSystem(nullptr)
////        /* init */,isGlissile(false)
//        {// Virtual dislocation loop
//        }
//
////        /**********************************************************************/
////        DislocationLoop(const DislocationLoop& other) :
////        /* base init */ BaseLoopType(other)
////        /* init */,slipSystem(other.slipSystem)
//////        /* init */,isGlissile(other.isGlissile)
////        {
////        }
//
////        /**********************************************************************/
////        void updateGeometry()
////        {
////            BaseLoopType::updateGeometry();
////            updateSlipSystem();
////        }
//
//        /**********************************************************************/
//        std::tuple<double,double,double> loopLength() const
//        {
//            double freeLength=0.0;
//            double boundaryLength=0.0;
//            double junctionLength=0.0;
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
//            return std::make_tuple(freeLength,junctionLength,boundaryLength);
//        }
//
//
//        /**********************************************************************/
//        template <class T>
//        friend T& operator << (T& os, const DislocationLoopType& dL)
//        {
//            os<< DislocationLoopIO<dim>(dL);
//            return os;
//        }
//
   };
    
}
#endif
