/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationLoop_H_
#define model_DislocationLoop_H_

#include <memory>

#include <model/LoopNetwork/Loop.h>
#include <model/LatticeMath/LatticeVector.h>
#include <model/LatticeMath/LatticePlaneBase.h>
#include <model/LatticeMath/LatticePlane.h>
#include <model/DislocationDynamics/Polycrystals/Grain.h>
#include <model/DislocationDynamics/GlidePlanes/GlidePlane.h>
#include <model/DislocationDynamics/GlidePlanes/GlidePlaneObserver.h>
#include <model/DislocationDynamics/IO/DislocationLoopIO.h>
#include <model/Geometry/PlanarPolygon.h>


namespace model
{
    template <int _dim, short unsigned int corder, typename InterpolationType>
    class DislocationLoop : public Loop<DislocationLoop<_dim,corder,InterpolationType>>,
                            public PlanarPolygon
    {
    
        Eigen::Matrix<double,_dim,1> nA;
        double _slippedArea;
        Eigen::Matrix<double,_dim,1> _rightHandedNormal;
        
    public:
        
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        
        constexpr static int dim=_dim;
        typedef DislocationLoop<dim,corder,InterpolationType> DislocationLoopType;
        typedef Loop<DislocationLoopType> BaseLoopType;
        typedef typename BaseLoopType::LoopLinkType LoopLinkType;
        typedef typename TypeTraits<DislocationLoopType>::LoopNetworkType LoopNetworkType;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        typedef GlidePlane<dim> GlidePlaneType;
        typedef GlidePlaneObserver<dim> GlidePlaneObserverType;
        typedef Eigen::Matrix<long int,dim+1,1> GlidePlaneKeyType;
        
        
        const Grain<dim>& grain;
        const std::shared_ptr<GlidePlaneType> _glidePlane;
        const GlidePlaneType& glidePlane;
        const bool isGlissile;
        
        
        /**********************************************************************/
        static bool allowedSlipSystem(const LatticeVector<dim>& b,
                                      const ReciprocalLatticeDirection<dim>& n,
                                      const Grain<dim>& gr)
        {
        
            std::cout<<"DislocationLoop::FINISH HERE. ANOTHER CONDITION FOR isGlissile SHOULD BE THAT N IS ONE OF THE ALLOWED SLIP PLANES"<<std::endl;

            return true;
        }

        
        /**********************************************************************/
        DislocationLoop(LoopNetworkType* const dn,
                        const VectorDim& B,
                        const VectorDim& N,
                        const VectorDim& P,
                        const int& grainID) :
        /* base init */ BaseLoopType(dn,dn->poly.grain(grainID).latticeVector(B)),
        /*      init */ PlanarPolygon(fabs(B.dot(N))<FLT_EPSILON? B : N.cross(VectorDim::Random()),N),
        /*      init */ nA(VectorDim::Zero()),
        /*      init */ _slippedArea(0.0),
        /*      init */ _rightHandedNormal(VectorDim::Zero()),
        /*      init */ grain(dn->poly.grain(grainID)),
//        /*      init */ _glidePlane(dn->sharedGlidePlane(dn->mesh,grain.lattice(),grain.grainID,grain.grainID,P,N)),
        /*      init */ _glidePlane(dn->sharedGlidePlane(dn->mesh,dn->poly.grain(grainID),P,N)),
        /*      init */ glidePlane(*_glidePlane.get()),
        /*      init */ isGlissile(this->flow().dot(glidePlane.n)==0 && allowedSlipSystem(this->flow(),glidePlane.n,grain))
        {
            model::cout<<"Creating DislocationLoop "<<this->sID<<", glissile="<<isGlissile<<std::endl;
            
//            _glidePlane->addLoop(this);
            _glidePlane->addParentSharedPtr(&_glidePlane,isGlissile,this->sID);

        }
        
        /**********************************************************************/
        DislocationLoop(const DislocationLoop& other) :
        /* base init */ BaseLoopType(other),
        /*      init */ PlanarPolygon(other),
        /* init */ nA(other.nA),
        /* init */ _slippedArea(0.0),
        /* init */ _rightHandedNormal(VectorDim::Zero()),
        /* init */ grain(other.grain),
        /* init */ _glidePlane(other._glidePlane),
        /* init */ glidePlane(*_glidePlane.get()),
        /* init */ isGlissile(other.isGlissile)
        {
//            model::cout<<"Copying DislocationLoop "<<this->sID<<std::endl;
            
//            _glidePlane->addLoop(this);
            _glidePlane->addParentSharedPtr(&_glidePlane,isGlissile,this->sID);

        }
        
        /**********************************************************************/
        ~DislocationLoop()
        {
//            model::cout<<"Destroying DislocationLoop "<<this->sID<<std::endl;

            
//            _glidePlane->removeLoop(this);
            _glidePlane->removeParentSharedPtr(&_glidePlane,isGlissile,this->sID);
        }
        
        /**********************************************************************/
        VectorDim burgers() const
        {
            return this->flow().cartesian();
        }
        
        /**********************************************************************/
        std::tuple<double,double,double> loopLength() const
        {
            double freeLength=0.0;
            double boundaryLength=0.0;
            double junctionLength=0.0;

            for(const auto& link : this->links())
            {
                if(link.second->pLink->isBoundarySegment())
                {
                    boundaryLength+=(link.second->sink()->get_P()-link.second->source()->get_P()).norm();
                }
                else
                {
                    if(link.second->pLink->loopLinks().size()==1)
                    {
                        freeLength+=(link.second->sink()->get_P()-link.second->source()->get_P()).norm();
                    }
                    else
                    {
                        junctionLength+=(link.second->sink()->get_P()-link.second->source()->get_P()).norm();
                    }
                }
            }
            
            return std::make_tuple(freeLength,junctionLength,boundaryLength);
        }

        /**********************************************************************/
        void update()
        {
            nA.setZero();
            for(const auto& loopLink : this->links())
            {
                nA+= 0.5*(loopLink.second->source()->get_P()-glidePlane.P).cross(loopLink.second->sink()->get_P()-loopLink.second->source()->get_P());
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
        
        /**********************************************************************/
        template <class T>
        friend T& operator << (T& os, const DislocationLoopType& dL)
        {
            os<< DislocationLoopIO<dim>(dL);
            return os;
        }
        
    };
    
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
////            assert(std::fabs((chord/chornNorm).dot(glidePlane.unitNormal))<FLT_EPSILON && "Chord does not belong to plane");
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
//    bool DislocationLoop<_dim,corder,InterpolationType>::outputLoopLength=false;

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
