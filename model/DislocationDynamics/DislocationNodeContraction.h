/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2018 by Yinan Cui  <cuiyinan@g.ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationNodeContraction_H_
#define model_DislocationNodeContraction_H_

#include <memory>
#include <tuple>
#include <Eigen/Dense>
#include <TypeTraits.h>
#include <GlidePlane.h>
//#include <BoundingMeshSegments.h>
#include <ConfinedDislocationObject.h>
#include <Grain.h>
#include <FiniteLineSegment.h>

#ifndef NDEBUG
#define VerboseNodeContraction(N,x) if(verboseNodeContraction>=N){model::cout<<x;}
#else
#define VerboseNodeContraction(N,x)
#endif

namespace model
{
    template <typename DislocationNetworkType>
    class DislocationNodeContraction
    {
        
        static constexpr int dim=DislocationNetworkType::dim;
        typedef typename DislocationNetworkType::NodeType NodeType;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        
        DislocationNetworkType& DN;

    public:
        
        const int verboseNodeContraction;

        /**********************************************************************/
        DislocationNodeContraction(DislocationNetworkType& DN_in) :
        /* init */ DN(DN_in)
        /* init */,verboseNodeContraction(TextFileParser("inputFiles/DD.txt").readScalar<int>("verboseNodeContraction",true))
        {
            
        }
        
        /**********************************************************************/
        bool contractSecondAndVirtual(std::shared_ptr<NodeType> nA,
                                      std::shared_ptr<NodeType> nB)
        {
            switch (DN.simulationParameters.simulationType)
            {
                case DefectiveCrystalParameters::FINITE_FEM:
                {
                    if(nB->virtualBoundaryNode())
                    {
                        assert(nA->virtualBoundaryNode());
                        DN.contractSecond(nA->virtualBoundaryNode(),nB->virtualBoundaryNode());
                    }
                    break;
                }
                    
//                case DefectiveCrystalParameters::PERIODIC:
//                {
//                 // FINISH HERE
//                    break;
//                }
            }
            return DN.contractSecond(nA,nB);
        }
        
        /**********************************************************************/
        bool contractYoungest(std::shared_ptr<NodeType> nA,
                              std::shared_ptr<NodeType> nB)
        {
            return nA->sID<nB->sID? contractSecondAndVirtual(nA,nB) : contractSecondAndVirtual(nB,nA);
            
        }
        
        /**********************************************************************/
        bool contractToPosition(std::shared_ptr<NodeType> nA,
                                std::shared_ptr<NodeType> nB,
                                const VectorDim& X,
                                const double& maxRange)
        {
            
            bool movedA=false;
            bool movedB=false;
            
            if(   nA->isMovableTo(X)
               && nB->isMovableTo(X)
               && (nA->get_P()-X).norm()+(nB->get_P()-X).norm()<maxRange)
            {
                movedA=nA->set_P(X);
                movedB=nB->set_P(X);
                VerboseNodeContraction(2,"contractToPosition"<<std::endl;);
                VerboseNodeContraction(2,"movedA="<<movedA<<std::endl;);
                VerboseNodeContraction(2,"movedB="<<movedB<<std::endl;);
                assert(movedA && movedB && "COULD NOT MOVE NODES");
            }
            
            return (movedA && movedB)? contractYoungest(nA,nB) : false;
        }
        
        /**********************************************************************/
        bool contract(std::shared_ptr<NodeType> nA,
                      std::shared_ptr<NodeType> nB)
        {
            
            VerboseNodeContraction(1,"DislocationNodeContraction::contract "<<nA->sID<<" "<<nB->sID<<std::endl;);
            
            const bool nAisMovable=nA->isMovableTo(nB->get_P());
            const bool nBisMovable=nB->isMovableTo(nA->get_P());
            
            if(nAisMovable && nBisMovable)
            {
                VerboseNodeContraction(1,"DislocationNodeContraction case 1a"<<std::endl;);
                return contractYoungest(nA,nB);
            }
            else if(nAisMovable && !nBisMovable)
            {
                VerboseNodeContraction(1,"DislocationNodeContraction case 1b"<<std::endl;);
                return contractSecondAndVirtual(nB,nA);
            }
            else if(!nAisMovable && nBisMovable)
            {
                VerboseNodeContraction(1,"DislocationNodeContraction case 1c"<<std::endl;);
                return contractSecondAndVirtual(nA,nB);
            }
            else
            {// nA and nB cannot be moved to each other. The calculation of a third point is necessary
                VerboseNodeContraction(1,"DislocationNodeContraction case 1d"<<std::endl;);

                const double maxRange=4.0*(nA->get_P()-nB->get_P()).norm();
                
                
                if(nA->isOnBoundary() || nB->isOnBoundary())
                {// either one of the nodes is a boundary node. Therefore the contraction point must be a boundary node
                    
                    //BoundingMeshSegments<dim> temp(nA->boundingBoxSegments(),nB->boundingBoxSegments());
                    const ConfinedDislocationObject<dim> cdo(*nA,*nB);
                    const BoundingMeshSegments<dim>& temp(cdo.boundingBoxSegments());

                    VerboseNodeContraction(1,"temp.size="<<temp.size()<<std::endl;);
                    
                    switch (temp.size())
                    {
                            
                        case 0:
                        {// no intersection of the bounding boxes
                            VerboseNodeContraction(1,"DislocationNodeContraction case 1f"<<std::endl;);
                            return false;
                            break;
                        }
                            
                        case 1:
                        {
                            const std::shared_ptr<MeshBoundarySegment<dim>>& mbs(temp.front());
                            if((mbs->P0-mbs->P1).norm()<FLT_EPSILON)
                            {// a unique intersection point of the bounding boxes exist
                                VerboseNodeContraction(1,"DislocationNodeContraction case 1d"<<std::endl;);
                                const VectorDim X=0.5*(mbs->P0+mbs->P1);
                                return contractToPosition(nA,nB,X,maxRange);
                            }
                            else
                            {// two possible intersection points of the bounding boxes exist
                                const bool firstIsCloser=(nA->get_P()-mbs->P0).norm()+(nB->get_P()-mbs->P0).norm()<(nA->get_P()-mbs->P1).norm()+(nB->get_P()-mbs->P1).norm();
                                const VectorDim X= firstIsCloser? mbs->P0 : mbs->P1;
                                const VectorDim Y= firstIsCloser? mbs->P1 : mbs->P0;
                                
                                VerboseNodeContraction(1,"DislocationNodeContraction case 1dX"<<std::endl;);
                                const bool Xcontracted=contractToPosition(nA,nB,X,maxRange);
                                if(Xcontracted)
                                {
                                    return true;
                                }
                                else
                                {
                                    VerboseNodeContraction(1,"DislocationNodeContraction case 1dY"<<std::endl;);
                                    const bool Ycontracted=contractToPosition(nA,nB,Y,maxRange);
                                    if(Ycontracted)
                                    {
                                        return true;
                                    }
                                    else
                                    {
                                        return false;
                                    }
                                }
                            }
                            break;
                        }
                            
                        default:
                        {// bounding boxes intersect in more than one line
                            
                            std::map<double,VectorDim> vertexMap;
                            for(const auto& seg : temp)
                            {
                                const double firstRange=(nA->get_P()-seg->P0).norm()+(nB->get_P()-seg->P0).norm();
                                if(firstRange<maxRange)
                                {
                                    vertexMap.insert(std::make_pair(firstRange,seg->P0));
                                }
                                
                                const double secondRange=(nA->get_P()-seg->P1).norm()+(nB->get_P()-seg->P1).norm();
                                if(secondRange<maxRange)
                                {
                                    vertexMap.insert(std::make_pair(secondRange,seg->P1));
                                }
                            }
                            
                            bool success=false;
                            for(const auto& pair : vertexMap)
                            {
                                success=contractToPosition(nA,nB,pair.second,maxRange);
                                if(success)
                                {
                                    break;
                                }
                            }
                            return success;
                            
                            break;
                        }
                            
                    }
                    
                }
                else
                {// neither nA nor nB are on bounding box
                    VerboseNodeContraction(1,"DislocationNodeContraction case 5"<<std::endl;);
                    if(nA->glidePlaneIntersections() && nB->glidePlaneIntersections())
                    {// both nodes confined by more then one plane
                        
                        SegmentSegmentDistance<dim> ssd(nA->glidePlaneIntersections()->P0,nA->glidePlaneIntersections()->P1,
                                                        nB->glidePlaneIntersections()->P0,nB->glidePlaneIntersections()->P1);
                        
                        const auto iSeg=ssd.intersectionSegment();
                        if(iSeg.size()==1)
                        {// incident intersection
                            VerboseNodeContraction(1,"DislocationNodeContraction case 5a"<<std::endl;);
                            return contractToPosition(nA,nB,std::get<0>(iSeg[0]),maxRange);
                        }
                        else if(iSeg.size()==2)
                        {// coincident intersection
                            VerboseNodeContraction(1,"DislocationNodeContraction case 5b"<<std::endl;);
                            FiniteLineSegment<dim> ls(std::get<0>(iSeg[0]),std::get<0>(iSeg[1]));
                            return contractToPosition(nA,nB,ls.snap(0.5*(nA->get_P()+nB->get_P())),maxRange);
                        }
                        else
                        {// parallel or skew intersection
                            VerboseNodeContraction(1,"DislocationNodeContraction case 5c"<<std::endl;);
                            return false;
                        }
                    }
                    else if(nA->glidePlaneIntersections() && !nB->glidePlaneIntersections())
                    {// nA confined by more then one plane, nB confined by only one plane
                        assert(nB->glidePlanes().size()==1);
                        PlaneLineIntersection<dim> pli((*nB->glidePlanes().begin())->P,
                                                       (*nB->glidePlanes().begin())->unitNormal,
                                                       nA->glidePlaneIntersections()->P0, // origin of line
                                                       nA->glidePlaneIntersections()->P1-nA->glidePlaneIntersections()->P0 // line direction
                                                       );
                        
                        // THERE SHOULD BE A PlaneSegmentIntersection class, which intersects the plane with a finite segment
                        
                        if(pli.type==PlaneLineIntersection<dim>::COINCIDENT)
                        {// nothing to do, _glidePlaneIntersections remains unchanged
                            VerboseNodeContraction(1,"DislocationNodeContraction case 6a"<<std::endl;);
                            return contractToPosition(nA,nB,nA->glidePlaneIntersections()->snap(0.5*(nA->get_P()+nB->get_P())),maxRange);
                        }
                        else if(pli.type==PlaneLineIntersection<dim>::INCIDENT)
                        {// _glidePlaneIntersections becomes a singular point
                            VerboseNodeContraction(1,"DislocationNodeContraction case 6b"<<std::endl;);
                            FiniteLineSegment<dim> cutLine(nA->glidePlaneIntersections()->P0,nA->glidePlaneIntersections()->P1);
                            if((pli.P-cutLine.snap(pli.P)).squaredNorm()<FLT_EPSILON)
                            {// intersection point is inside mesh
                                VerboseNodeContraction(1,"DislocationNodeContraction case 6b1"<<std::endl;);
                                return contractToPosition(nA,nB,pli.P,maxRange);
                            }
                            else
                            {
                                return false;
                            }
                        }
                        else
                        {// parallel planes, cannot contract
                            VerboseNodeContraction(1,"DislocationNodeContraction case 6c"<<std::endl;);
                            return false;
                        }
                    }
                    else if(!nA->glidePlaneIntersections() && nB->glidePlaneIntersections())
                    {// same as previous case, so switch nA and nB
                        VerboseNodeContraction(1,"DislocationNodeContraction case 7a"<<std::endl;);
                        return contract(nB,nA);
                    }
                    else
                    {// both nodes confined by only one plane
                        
                        
                        if(nA->glidePlanes().size()==1 && nB->glidePlanes().size()==1)
                        {
                            const PlanePlaneIntersection<dim>& ppi(DN.glidePlaneFactory.glidePlaneIntersection(*nA->glidePlanes().begin(),*nB->glidePlanes().begin()));
                            
                            if(ppi.type==PlanePlaneIntersection<dim>::COINCIDENT)
                            {// the contraction point can be the averago of nA and nB, which should be internal for convex domains
                                VerboseNodeContraction(1,"DislocationNodeContraction case 8a"<<std::endl;);
                                return contractToPosition(nA,nB,nA->glidePlaneIntersections()->snap(0.5*(nA->get_P()+nB->get_P())),maxRange);
                            }
                            else if(ppi.type==PlanePlaneIntersection<dim>::INCIDENT)
                            {
                                VerboseNodeContraction(1,"DislocationNodeContraction case 8b"<<std::endl;);
                                const ConfinedDislocationObject<dim> cdo(*nA,*nB);
                                const BoundingMeshSegments<dim>& temp(cdo.boundingBoxSegments());
                                switch (temp.size())
                                {
                                    case 2:
                                    {
                                        const std::shared_ptr<MeshBoundarySegment<dim>>& seg0(temp.front());
                                        const std::shared_ptr<MeshBoundarySegment<dim>>& seg1(temp.back());
                                        
                                        FiniteLineSegment<dim> cutLine(0.5*(seg0->P0+seg0->P1),0.5*(seg1->P0+seg1->P1));
                                        return contractToPosition(nA,nB,cutLine.snap(0.5*(nA->get_P()+nB->get_P())),maxRange);
                                        break;
                                    }
                                    case 1:
                                    {
                                        const std::shared_ptr<MeshBoundarySegment<dim>>& seg0(temp.front());
                                        
                                        FiniteLineSegment<dim> cutLine(seg0->P0,seg0->P1);
                                        return contractToPosition(nA,nB,cutLine.snap(0.5*(nA->get_P()+nB->get_P())),maxRange);
                                        break;
                                    }
                                    default:
                                    {// intersection line outside mesh
                                        return false;
                                    }
                                }
                            }
                            else
                            {
                                VerboseNodeContraction(1,"DislocationNodeContraction case 8c"<<std::endl;);
                                return false;
                            }
                        }
                        else if(nA->glidePlanes().size()==0 && nB->glidePlanes().size()==1)
                        {// nA has no GlidePlane
                            const VectorDim dir(nA->invariantDirectionOfMotion());
                            if(dir.norm()>FLT_EPSILON)
                            {// nB can be moved along dir
                                PlaneLineIntersection<dim> pli((*nB->glidePlanes().begin())->P,
                                                               (*nB->glidePlanes().begin())->unitNormal,
                                                               nA->get_P(), // origin of line
                                                               dir // line direction
                                                               );
                                
                                if(pli.type==PlaneLineIntersection<dim>::COINCIDENT)
                                {// nothing to do, _glidePlaneIntersections remains unchanged
                                    VerboseNodeContraction(1,"DislocationNodeContraction case 9a"<<std::endl;);
                                    return contractToPosition(nA,nB,nA->get_P(),maxRange);
                                }
                                else if(pli.type==PlaneLineIntersection<dim>::INCIDENT)
                                {// _glidePlaneIntersections becomes a singular point
                                    VerboseNodeContraction(1,"DislocationNodeContraction case 9b"<<std::endl;);
                                    FiniteLineSegment<dim> cutLine(nA->glidePlaneIntersections()->P0,nA->glidePlaneIntersections()->P1);
                                    return contractToPosition(nA,nB,pli.P,maxRange);
//
//                                    if((pli.P-cutLine.snap(pli.P)).squaredNorm()<FLT_EPSILON)
//                                    {// intersection point is inside mesh
//                                        VerboseNodeContraction(1,"DislocationNodeContraction case96b1"<<std::endl;);
//                                    }
//                                    else
//                                    {
//                                        return false;
//                                    }
                                }
                                else
                                {// parallel planes, cannot contract
                                    VerboseNodeContraction(1,"DislocationNodeContraction case 9c"<<std::endl;);
                                    return false;
                                }
                            }
                            else
                            {
                                return false;
                            }
                        }
                        else if(nA->glidePlanes().size()==0 && nB->glidePlanes().size()==0)
                        {// nB has no GlidePlane
                            return false;
                        }
                        else
                        {// neither nA nor nB have a GlidePlane
                            assert(nA->glidePlanes().size()<=1);
                            assert(nB->glidePlanes().size()<=1);
                            return false;
                        }
                        

                    }
                }
            }
            
        }
    };
    
    // Static data
//    template <typename DislocationNetworkType>
//    int DislocationNodeContraction<DislocationNetworkType>::verboseNodeContraction=0;
    
}
#endif


