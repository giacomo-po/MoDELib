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
#include <model/Utilities/TypeTraits.h>
#include <model/DislocationDynamics/GlidePlanes/GlidePlane.h>
#include <model/DislocationDynamics/BoundingLineSegments.h>
#include <model/DislocationDynamics/Polycrystals/Grain.h>
#include <model/Geometry/LineSegment.h>

#ifndef NDEBUG
#define VerboseNodeContraction(N,x) if(verboseNodeContraction>=N){model::cout<<x;}
#else
#define VerboseNodeContraction(N,x)
#endif

namespace model
{
    template <typename DislocationNetworkType>
    struct DislocationNodeContraction
    {
        
        static constexpr int dim=DislocationNetworkType::dim;
        typedef typename DislocationNetworkType::NodeType NodeType;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        
        static int verboseNodeContraction;
        
        
        DislocationNetworkType& DN;
        
        /**********************************************************************/
        DislocationNodeContraction(DislocationNetworkType& DN_in) : DN(DN_in)
        {
            
        }
        
        /**********************************************************************/
//old version of the contractYoungest
//        bool contractYoungest(std::shared_ptr<NodeType> nA,
//                              std::shared_ptr<NodeType> nB)
//        {
//            return nA->sID<nB->sID? DN.contractSecond(nA,nB) : DN.contractSecond(nB,nA);
//
//        }
// Chnaged with respect to the new implementation.
        bool contractYoungest(std::shared_ptr<NodeType> nA,
        		std::shared_ptr<NodeType> nB)
        {
        if ((nA->imageBoundaryNodeContainer().size()>0 || nB->imageBoundaryNodeContainer().size()>0) || (nA->isImageBoundaryNode || nB->isImageBoundaryNode))
        	return contractSecondWithImages(nA,nB);
        else
        	return nA->sID<nB->sID? DN.contractSecond(nA,nB) : DN.contractSecond(nB,nA);

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
        
        bool contractSecondWithImages(std::shared_ptr<NodeType> nA,
                             std::shared_ptr<NodeType> nB)
        {	if ((!nA->isImageBoundaryNode && !nB->isImageBoundaryNode)
        		&& (nA->imageBoundaryNodeContainer().size()>0 || nB->imageBoundaryNodeContainer().size()>0))
        	{//First contract the images and then contract the real nodes
        	size_t temp=0;
        	if (nA->imageBoundaryNodeContainer().size()>0 && nB->imageBoundaryNodeContainer().size()>0)
        	{ size_t count=0;

        		for (const auto& iternA:nA->imageBoundaryNodeContainer())
        		{

        			for (const auto& iternB:nB->imageBoundaryNodeContainer())
        			{
        				if (iternA.first==iternB.first)
        				{
        					if(DN.contractSecond(iternA.second,iternB.second));
        					std::cout<<"Contracting image node"<<iternB.second->sID<<std::endl;
        					count++;
        				}
        			}
        		}
        		if (count==std::min(nA->imageBoundaryNodeContainer().size(),nB->imageBoundaryNodeContainer().size())) //TODO after contraction, container size is samwe
        			temp=1;
        		if (temp)
        			{

        				std::cout<<"Contracting "<<nB->sID<<std::endl;
        				return DN.contractSecond(nA,nB);
        			}
        		else
        			assert(0 && "All the image nodes could not be contracted");
        	}
        	else if (nA->imageBoundaryNodeContainer().size()>0 && nB->imageBoundaryNodeContainer().size()==0)
        	{
        		model::cout<<greenColor<<"executing else statement in contract second with images"<<std::endl;
//        		temp=(nA->imageBoundaryNodeContainer().size()>0)?DN.contractSecond(nB,nA):DN.contractSecond(nA,nB);
        		temp=DN.contractSecond(nA,nB);
        		return temp;
        	}
        	else if (nA->imageBoundaryNodeContainer().size()==0 && nB->imageBoundaryNodeContainer().size()>0)
        	{
        		model::cout<<greenColor<<"executing else statement in contract second with images"<<std::endl;
        		//        		temp=(nA->imageBoundaryNodeContainer().size()>0)?DN.contractSecond(nB,nA):DN.contractSecond(nA,nB);
        		temp=DN.contractSecond(nB,nA);		//TODO check this if this can change the contract sequence
        		return temp;
        	}
        	else
        	{
        		assert(0 && "Case for the real boundary node within contract second not taken into consideration");
        	}

        	}
       //TODO: Check the condition for the contraction of the image node and real node for one more time
        else if ((nA->isImageBoundaryNode && !nB->isImageBoundaryNode))
        {
        	return DN.contractSecond(nA,nB);
        }

        else if ((!nA->isImageBoundaryNode && nB->isImageBoundaryNode))
        {
//        	for (const auto& nA1:nA->imageBoundaryNodeContainer())
//        	{
//        		for (const auto& nB1:nB->masterNode->imageBoundaryNodeContainer())
//        		{
//        			if ( (nA1.first==nB1.first) && ( nB1.second->sID==nB->sID))
//        			{
//        				//condition satisfied for the movement of the real node
//        				if (nB->masterNode->isMovableTo(nA1.second->get_P()))
//        				{
//        					auto nB2=DN.sharedNode(nB->masterNode->sID);
//        					//contract the master node
//        					contractSecondWithImages(nA1.second,nB2.second);
////        					DN.contract(nA1.second->sID,nB->masterNode->sID);
//        				}
//        			}
//        		}
//        	}
        	//check if nA is movable to nB
        	if (nA->isMovableTo(nB->get_P()))
        	{
        		return DN.contractSecond(nB,nA);
        	}
        	else
        	{
        		assert (0 && "All conditions for movable not taken into consideration in Dislocation Node contraction with images");
        	}
        }
        else if ((nA->isImageBoundaryNode && nB->isImageBoundaryNode))
        {
        	//go to master node and check if nA and nB directions are the same
        	for (const auto& nA1:nA->masterNode->imageBoundaryNodeContainer())
        	{
        		for (const auto& nB1:nB->masterNode->imageBoundaryNodeContainer())
        		{
        			if ((nA1.first==nB1.first) && (nA1.second->sID==nA->sID && nB1.second->sID==nB->sID))
        			{
        				//condition satisfied for the movement of the real node
        				if (nB->masterNode->isMovableTo(nA->masterNode->get_P()))
        				{
        					//contract the master node
        					auto nA2=DN.sharedNode(nA->masterNode->sID);
        					auto nB2=DN.sharedNode(nB->masterNode->sID);

        					contractSecondWithImages(nA2.second,nB2.second);
        				}
        			}
        		}
        	}
        }
        else if (nA->imageBoundaryNodeContainer().size()==0 && nB->imageBoundaryNodeContainer().size()==0)
        	{
        		//just contract the real nodes
        		return DN.contractSecond(nA,nB);

        	}
        else
        {
        	assert(0 && "Condition for Node contraction with images not taken into consideration.");
        }
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
                return contractSecondWithImages(nB,nA);
              // return contractSecond(nB,nA);

            }
            else if(!nAisMovable && nBisMovable)
            {
                VerboseNodeContraction(1,"DislocationNodeContraction case 1c"<<std::endl;);
                return contractSecondWithImages(nA,nB);
//                return contractSecond(nA,nB);
            }
            else
            {// nA and nB cannot be moved to each other. The calculation of a third point is necessary
                
                const double maxRange=4.0*(nA->get_P()-nB->get_P()).norm();
                
                
                if(nA->isOnBoundingBox() || nB->isOnBoundingBox())
                {// either one of the nodes is a boundary node. Therefore the contraction point must be a boundary node
                    
                    BoundingLineSegments<dim> temp(nA->boundingBoxSegments(),nB->boundingBoxSegments());
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
                            if((temp[0].first-temp[0].second).norm()<FLT_EPSILON)
                            {// a unique intersection point of the bounding boxes exist
                                VerboseNodeContraction(1,"DislocationNodeContraction case 1d"<<std::endl;);
                                const VectorDim X=0.5*(temp[0].first+temp[0].second);
                                return contractToPosition(nA,nB,X,maxRange);
                            }
                            else
                            {// two possible intersection points of the bounding boxes exist
                                const bool firstIsCloser=(nA->get_P()-temp[0].first).norm()+(nB->get_P()-temp[0].first).norm()<(nA->get_P()-temp[0].second).norm()+(nB->get_P()-temp[0].second).norm();
                                const VectorDim X= firstIsCloser? temp[0].first : temp[0].second;
                                const VectorDim Y= firstIsCloser? temp[0].second : temp[0].first;
                                
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
                                const double firstRange=(nA->get_P()-seg.first).norm()+(nB->get_P()-seg.first).norm();
                                if(firstRange<maxRange)
                                {
                                    vertexMap.insert(std::make_pair(firstRange,seg.first));
                                }
                                
                                const double secondRange=(nA->get_P()-seg.second).norm()+(nB->get_P()-seg.second).norm();
                                if(secondRange<maxRange)
                                {
                                    vertexMap.insert(std::make_pair(secondRange,seg.second));
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
                {// neither a nor b are on bounding box
                    if(nA->glidePlaneIntersections().size() && nB->glidePlaneIntersections().size())
                    {// both nodes confined by more then one plane
                        
                        SegmentSegmentDistance<dim> ssd(nA->glidePlaneIntersections()[0].first,nA->glidePlaneIntersections()[0].second,
                                                        nB->glidePlaneIntersections()[0].first,nB->glidePlaneIntersections()[0].second);
                        
                        const auto iSeg=ssd.intersectionSegment();
                        if(iSeg.size()==1)
                        {// incident intersection
                            VerboseNodeContraction(1,"DislocationNodeContraction case 5a"<<std::endl;);
                            return contractToPosition(nA,nB,std::get<0>(iSeg[0]),maxRange);
                        }
                        else if(iSeg.size()==2)
                        {// coincident intersection
                            VerboseNodeContraction(1,"DislocationNodeContraction case 5b"<<std::endl;);
                            LineSegment<dim> ls(std::get<0>(iSeg[0]),std::get<0>(iSeg[1]));
                            return contractToPosition(nA,nB,ls.snap(0.5*(nA->get_P()+nB->get_P())),maxRange);
                        }
                        else
                        {// parallel or skew intersection
                            VerboseNodeContraction(1,"DislocationNodeContraction case 5c"<<std::endl;);
                            return false;
                        }
                    }
                    else if(nA->glidePlaneIntersections().size() && !nB->glidePlaneIntersections().size())
                    {// nA confined by more then one plane, nB confined by only one plane
                        PlaneLineIntersection<dim> pli(nB->meshPlane(0).P,
                                                       nB->meshPlane(0).unitNormal,
                                                       nA->glidePlaneIntersections()[0].first, // origin of line
                                                       nA->glidePlaneIntersections()[0].second-nA->glidePlaneIntersections()[0].first // line direction
                                                       );
                        
                        if(pli.type==PlaneLineIntersection<dim>::COINCIDENT)
                        {// nothing to do, _glidePlaneIntersections remains unchanged
                            VerboseNodeContraction(1,"DislocationNodeContraction case 6a"<<std::endl;);
                            return contractToPosition(nA,nB,nA->snapToMeshPlaneIntersection(0.5*(nA->get_P()+nB->get_P())),maxRange);
                        }
                        else if(pli.type==PlaneLineIntersection<dim>::INCIDENT)
                        {// _glidePlaneIntersections becomes a singular point
                            VerboseNodeContraction(1,"DislocationNodeContraction case 6b"<<std::endl;);
                            return contractToPosition(nA,nB,pli.P,maxRange);
                        }
                        else
                        {// parallel planes, cannot contract
                            VerboseNodeContraction(1,"DislocationNodeContraction case 6c"<<std::endl;);
                            return false;
                        }
                    }
                    else if(!nA->glidePlaneIntersections().size() && nB->glidePlaneIntersections().size())
                    {
                        VerboseNodeContraction(1,"DislocationNodeContraction case 7a"<<std::endl;);
                        return contract(nB,nA);
                    }
                    else
                    {// both nodes confined by only one plane
                        
                        assert(nA->meshPlanes().size()==1);
                        assert(nB->meshPlanes().size()==1);
                        
                        const PlanePlaneIntersection<dim>& ppi(DN.glidePlaneIntersection(&nA->meshPlane(0),&nB->meshPlane(0)));
                        
                        if(ppi.type==PlanePlaneIntersection<dim>::COINCIDENT)
                        {
                            VerboseNodeContraction(1,"DislocationNodeContraction case 8a"<<std::endl;);
                            return contractToPosition(nA,nB,nA->snapToMeshPlaneIntersection(0.5*(nA->get_P()+nB->get_P())),maxRange);
                        }
                        else if(ppi.type==PlanePlaneIntersection<dim>::INCIDENT)
                        {
                            VerboseNodeContraction(1,"DislocationNodeContraction case 8b"<<std::endl;);
                            const double u=(0.5*(nA->get_P()+nB->get_P())-ppi.P).dot(ppi.d);
                            return contractToPosition(nA,nB,ppi.P+u*ppi.d,maxRange);
                        }
                        else
                        {
                            VerboseNodeContraction(1,"DislocationNodeContraction case 8c"<<std::endl;);
                            return false;
                        }
                    }
                }
            }
            
        }
    };
    
    // Static data
    template <typename DislocationNetworkType>
    int DislocationNodeContraction<DislocationNetworkType>::verboseNodeContraction=0;
    

}
#endif


