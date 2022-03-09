/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po         <gpo@ucla.edu>.
 * Copyright (C) 2011 by Benjamin Ramirez   <ramirezbrf@gmail.com>.
 * Copyright (C) 2011 by Tamer Crsoby       <tamercrosby@gmail.com>,
 * Copyright (C) 2011 by Can Erel           <canerel55@gmail.com>,
 * Copyright (C) 2011 by Mamdouh Mohamed    <msm07d@fsu.edu>
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationNode_cpp_
#define model_DislocationNode_cpp_

#include <DislocationNode.h>


namespace model
{
    
    template <int dim, short unsigned int corder>
    DislocationNode<dim,corder>::DislocationNode(LoopNetworkType* const net,
                                                                   const VectorDim& P,
                                                                   const VectorDim& V,
                                                                   const double& vrc) :
    /* init */ NetworkNode<DislocationNode>(net)
    /* init */,SplineNodeType(P)
    // /* init */,ConfinedDislocationObjectType(this->get_P())
    /* init */,p_Simplex(get_includingSimplex(this->get_P(),(const Simplex<dim,dim>*) NULL))
    /* init */,velocity(V)
    /* init */,vOld(velocity)
    /* init */,velocityReductionCoeff(vrc)
    /* init */,virtualNode(nullptr)
    /* init */,masterNode(nullptr)
    {
        VerboseDislocationNode(1, "  Creating Network Node " << this->tag() <<" @ "<<this->get_P().transpose() << std::endl;);
    }

    template <int dim, short unsigned int corder>
    DislocationNode<dim,corder>::~DislocationNode()
    {
        VerboseDislocationNode(1, "  Destroying Network Node " << this->tag() << std::endl;);

    }
    
    template <int dim, short unsigned int corder>
    std::shared_ptr<DislocationNode<dim,corder>> DislocationNode<dim,corder>::clone() const
    {
        return std::shared_ptr<DislocationNode<dim,corder>>(new DislocationNode(this->p_network(),this->get_P(),get_V(),velocityReduction()));
    }
    
    template <int dim, short unsigned int corder>
    const Simplex<dim,dim>* DislocationNode<dim,corder>::get_includingSimplex(const VectorDim& X,const Simplex<dim,dim>* const guess) const
    {
        std::pair<bool,const Simplex<dim,dim>*> temp(false,NULL);
        if (guess==NULL)
        {
            temp=this->network().mesh.search(X);
        }
        else
        {
            const auto grains(this->grains());
            if(grains.size()==1)
            {// node only in one region
                if((*grains.begin())->grainID!=guess->region->regionID)
                {
                    temp=this->network().mesh.searchRegion((*grains.begin())->grainID,X);
                }
                else
                {
                    temp=this->network().mesh.searchRegionWithGuess(X,guess);
                }
            }
            else
            {
                temp=this->network().mesh.searchWithGuess(X,guess);
            }
        }
        if(!temp.first) // PlanarDislocationNode not found inside mesh
        {// Detect if the PlanarDislocationNode is sligtly outside the boundary
            int faceID;
            const double baryMin(temp.second->pos2bary(X).minCoeff(&faceID));
            const bool isApproxOnBoundary(std::fabs(baryMin)<1.0e3*FLT_EPSILON && temp.second->child(faceID).isBoundarySimplex());
            if(!isApproxOnBoundary)
            {
                std::cout<<"PlanarDislocationNode "<<this->sID<<" @ "<<X.transpose()<<std::endl;
                std::cout<<"Simplex "<<temp.second->xID<<std::endl;
                std::cout<<"bary "<<temp.second->pos2bary(X)<<std::endl;
                std::cout<<"face of barymin is "<<temp.second->child(faceID).xID<<std::endl;
                std::cout<<"face of barymin is boundary Simplex? "<<temp.second->child(faceID).isBoundarySimplex()<<std::endl;
                std::cout<<"face of barymin is region-boundary Simplex? "<<temp.second->child(faceID).isRegionBoundarySimplex()<<std::endl;
                assert(0 && "DISLOCATION NODE OUTSIDE MESH.");
            }
        }
        return temp.second;
    }

    // template <int dim, short unsigned int corder>
    // void DislocationNode<dim,corder>::addLoopNode(LoopNodeType* const pN)
    // {
    //     NetworkNode<DislocationNode>::addLoopNode(pN);

    //     //Check if the node is connected
    //     if (pN->next.first && pN->prev.first)
    //     {
    //         if (pN->periodicPlanePatch())
    //         {
    //             this->addGlidePlane(pN->periodicPlanePatch()->glidePlane.get());
    //         }
    //     }

    // }

    // template <int dim, short unsigned int corder>
    // void DislocationNode<dim,corder>::removeLoopNode(LoopNodeType* const pN)
    // {
    //     NetworkNode<DislocationNode>::removeLoopNode(pN);
    //     this->glidePlanes().clear();
    //     for (const auto &loopN : this->loopNodes())
    //     {
    //         if (loopN->next.first && loopN->prev.first)
    //         {
    //             if (loopN->periodicPlanePatch())
    //             {
    //                 this->addGlidePlane(loopN->periodicPlanePatch()->glidePlane.get());
    //             }
    //         }
    //     }
    // }
    
    
    template <int dim, short unsigned int corder>
    void DislocationNode<dim,corder>::projectVelocity()
    {
        
        VectorOfNormalsType temp;
        
        for(const auto& loopNode : this->loopNodes())
        {
//            for(const auto& loop : loopNode->loops())
//            {
            const auto& loop(loopNode->loop());
            if(loop)
            {
                if(loop->loopType==DislocationLoopIO<dim>::GLISSILELOOP)
                {
                    if(loop->glidePlane)
                    {
                        temp.push_back(loop->glidePlane->unitNormal);
                    }
                }
                else if(loop->loopType==DislocationLoopIO<dim>::SESSILELOOP)
                {
                    velocity.setZero();
                    break;
                }
            }
            else
            {
                assert(false && "no loop in loop");
            }

//            }

        }
        
        if(this->glidePlanes().size()>=dim)
        {
            velocity.setZero();
        }
        
        if(velocity.squaredNorm()>FLT_EPSILON)
        {
            
            
            if(!this->network().simulationParameters.isPeriodicSimulation())
            {// Use boundary planes to confine velocity in case of non-periodic simulation
                
                for(const auto& face : this->meshFaces())
                {
                    temp.push_back(face->asPlane().unitNormal);
                }
            }

            // std::cout<<this->sID<<" Temp size for gramSchmidt before"<<temp.size()<<std::endl;
            
            GramSchmidt::orthoNormalize(temp);
            
            // std::cout<<this->sID<<" Temp size for gramSchmidt after"<<temp.size()<<std::endl;

            for(const auto& vec : temp)
            {
                velocity-=velocity.dot(vec)*vec;
            }
            
        }
    }
    
    template <int dim, short unsigned int corder>
    void DislocationNode<dim,corder>::set_V(const VectorDim& vNew)
    {
        vOld=velocity; // store current value of velocity before updating
        
        //            velocity=this->prjM*vNew; // kill numerical errors from the iterative solver
        velocity=vNew;
        projectVelocity();
        
        if(this->network().use_velocityFilter && !isBoundaryNode())
        {
            const double filterThreshold=0.05*velocity.norm()*vOld.norm()+FLT_EPSILON;
            
            if(velocity.dot(vOld)<-filterThreshold)
            {
                velocityReductionCoeff*=this->network().velocityReductionFactor;
            }
            else if(velocity.dot(vOld)>filterThreshold)
            {
                velocityReductionCoeff/=this->network().velocityReductionFactor;
            }
            else
            {
                // don't change velocityReductionCoeff
            }
            if(velocityReductionCoeff>1.0)
            {
                velocityReductionCoeff=1.0;
            }
            if(velocityReductionCoeff<0.005)
            {
                velocityReductionCoeff=0.005;
            }
            velocity*=velocityReductionCoeff;
        }
        // std::cout<<"Time integration methdo is "<<this->network().simulationParameters.timeIntegrationMethod<<std::endl;
        // if (this->network().simulationParameters.timeIntegrationMethod==0)
        // {
        //     const double v_crit(0.1); //TODO : Read from the text file
        //     if (velocity.norm()>v_crit)
        //     {
        //         const double alpha(0.1/velocity.norm());
        //         velocity*=alpha;
        //     }
        // }
    }
    


    // template <int dim, short unsigned int corder>
    // void DislocationNode<dim,corder>::initFromFile(const std::string& fileName)
    // {
    //     this->network().use_velocityFilter=TextFileParser(fileName).readScalar<double>("this->network().use_velocityFilter",true);
    //     this->network().velocityReductionFactor=TextFileParser(fileName).readScalar<double>("this->network().velocityReductionFactor",true);
    //     assert(this->network().velocityReductionFactor>0.0 && this->network().velocityReductionFactor<=1.0);
    //     verboseDislocationNode=TextFileParser(fileName).readScalar<int>("verboseDislocationNode",true);
    // }

    template <int dim, short unsigned int corder>
    void DislocationNode<dim,corder>::updateGeometry()
    {
        for(const auto& neighboor : this->neighbors())
        {
            std::get<1>(neighboor.second)->updateGeometry();
        }
    }

    template <int dim, short unsigned int corder>
    bool DislocationNode<dim,corder>::trySet_P(const typename DislocationNode<dim,corder>::VectorDim& newP)
    {
        VerboseDislocationNode(1, " Try  Setting P for Network Node " << this->tag() << std::endl;);
        const VectorDim snappedPosition(this->snapToGlidePlanesinPeriodic(newP));
        if((this->get_P()-snappedPosition).norm()>FLT_EPSILON)
        {
            for (auto &loopNode : this->loopNodes())
            {
                loopNode->set_P(loopNode->periodicPlanePatch() ? snappedPosition - loopNode->periodicPlanePatch()->shift : snappedPosition);
            }
        }

        return true;
    }
    
    // template <int dim, short unsigned int corder>
    template <int dim, short unsigned int corder>
    bool DislocationNode<dim,corder>::set_P(const typename DislocationNode<dim,corder>::VectorDim& newP)
    {
        // std::cout<<" Current position "<<this->get_P().transpose()<<std::endl;
        // std::cout<<" New position "<<newP.transpose()<<std::endl;
        VerboseDislocationNode(1, "  Setting P for Network Node " << this->tag() << std::endl;);

        if((SplineNodeType::get_P()-newP).norm()>FLT_EPSILON)   //Check with the base classs
        {
            SplineNodeType::set_P(newP);
            // this->updateConfinement(this->get_P());
            p_Simplex=get_includingSimplex(this->get_P(),p_Simplex); // update including simplex
            
        }
        
        return (this->get_P()-newP).norm()<FLT_EPSILON;
    }


    
    template <int dim, short unsigned int corder>
    const typename DislocationNode<dim,corder>::VectorDim& DislocationNode<dim,corder>::get_V() const
    {/*! The nodal velocity vector
      */
        return velocity;
    }

    // /**********************************************************************/
    // template <int dim, short unsigned int corder>
    // bool DislocationNode<dim,corder>::isZeroBurgersNode() const
    // {
    //     VerboseDislocationNode(4,"DislocationNode "<<this->tag()<<" isZeroBurgersNode "<<std::flush;);
    //     bool temp = true;
    //     for (const auto &neighborIter : this->neighbors())
    //     {
    //         temp *= std::get<1>(neighborIter.second)->hasZeroBurgers();
    //         if (!temp)
    //         {
    //             break;
    //         }
    //     }
        
    //     VerboseDislocationNode(4,temp<<std::endl;);

    //     return temp;
    // }

    // template <int dim, short unsigned int corder>
    // const typename DislocationNode<dim,corder>::VectorDim& DislocationNode<dim,corder>::get_P() const
    // {/*! The nodal velocity vector
    //   */
    //     VerboseDislocationNode(2, "  Getting P for Network Node " << this->tag() << std::endl;);

    //     const VectorDim netNodePosition(SplineNodeType::get_P());
    //     // check the consistency between the loop Nodes and networkNodes
    //     for (const auto& lNode : this->loopNodes())
    //     {
    //         if (!(this->network().simulationParameters.isPeriodicSimulation() && lNode->periodicPlanePatch()) && lNode->get_P().norm()>FLT_EPSILON)
    //         {
    //             const VectorDim patchShift(lNode->periodicPlanePatch() ? lNode->periodicPlanePatch()->shift : VectorDim::Zero());
    //             const VectorDim lNodePos(lNode->get_P() + patchShift);
    //             if ((lNodePos - netNodePosition).squaredNorm() > FLT_EPSILON)
    //             {
    //                 std::cout << " Mistmatch for networknode " << this->tag() << " at loop Node " << lNode->tag() << std::endl;
    //                 std::cout << "netNodePosition => lNodePos " << netNodePosition.transpose() << "=>" << lNodePos.transpose() << std::endl;
    //                 assert(false && " Position mismatch ");
    //             }
    //         }
            
    //     }
    //     return (SplineNodeType::get_P());
    // }

    template <int dim, short unsigned int corder>
    const double& DislocationNode<dim,corder>::velocityReduction() const
    {
        return velocityReductionCoeff;
    }
    
    template <int dim, short unsigned int corder>
    typename DislocationNode<dim,corder>::MeshLocation DislocationNode<dim,corder>::meshLocation() const
    {/*!\returns the position of *this relative to the bonudary:
      * 1 = inside mesh
      * 2 = on mesh boundary
      */
        MeshLocation temp = MeshLocation::outsideMesh;

        if(!isVirtualBoundaryNode())
        {
            if(isBoundaryNode())
            {
                temp=MeshLocation::onMeshBoundary;
            }
            else
            {
                if(isGrainBoundaryNode())
                {
                    temp=MeshLocation::onRegionBoundary;
                }
                else
                {
                    temp=MeshLocation::insideMesh;
                }
            }
        }
        return temp;
    }
    
    template <int dim, short unsigned int corder>
    const std::shared_ptr<typename DislocationNode<dim,corder>::NetworkNodeType>& DislocationNode<dim,corder>::virtualBoundaryNode() const
    {
        return virtualNode;
    }
    
    template <int dim, short unsigned int corder>
    typename DislocationNode<dim,corder>::VectorDim DislocationNode<dim,corder>::invariantDirectionOfMotion() const
    {/*!\returns the direction of alignment if all links connected to this node are geometrically aligned.
      * Otherwise it returns the zero vector.
      */
        VectorDim temp(this->neighbors().size()? std::get<1>(this->neighbors().begin()->second)->chord().normalized() : VectorDim::Zero());
        for (const auto& neighborIter : this->neighbors())
        {
            if(std::get<1>(neighborIter.second)->chord().normalized().cross(temp).norm()>FLT_EPSILON)
            {
                temp.setZero();
            }
        }
        return temp;
    }
    
    template <int dim, short unsigned int corder>
    bool DislocationNode<dim,corder>::isMovableTo(const VectorDim& X) const
    {
        for(const auto& loopNode : this->loopNodes())
        {
            const VectorDim shift(loopNode->periodicPlanePatch()? loopNode->periodicPlanePatch()->shift : VectorDim::Zero());
            
            if(!loopNode->isMovableTo(X-shift))
            {
                VerboseDislocationNode(2,"  DislocationNode "<<this->tag()<< " NOT movable "<<std::endl;);
                return false;
            }
        }
        return true;
        
//        if(this->network().simulationParameters.isPeriodicSimulation() && this->isBoundaryNode())
//        {// cannot move boundary nodes under periodic simulations
//            return (X-this->get_P()).squaredNorm()<FLT_EPSILON;
//        }
//        else
//        {
//            bool isMovable=true;
//
//            VerboseDislocationNode(4,"checking if PlanarDislocationNode "<<this->sID<< " isMovable:"<<std::endl;);
//
//            for(const auto& gp : this->glidePlanes())
//            {// X must be contained by all glidePlanes
//                isMovable*=gp->contains(X);
//            }
//            VerboseDislocationNode(4,"  meshPlanes contains X? "<<isMovable<<std::endl;);
//
//            if(isMovable)
//            {
//                for(const auto& pair : this->neighbors())
//                {
//                    if(std::get<1>(pair.second)->isSessile())
//                    {// sessile segments cannot change direction if this node is moved
//                        const double currentNorm((std::get<0>(pair.second)->get_P()-this->get_P()).norm());
//                        const double newNorm((std::get<0>(pair.second)->get_P()-X).norm());
//                        if(currentNorm>FLT_EPSILON && newNorm>FLT_EPSILON)
//                        {
//                            const bool sessileNeighborMovable=((std::get<0>(pair.second)->get_P()-X).cross(std::get<0>(pair.second)->get_P()-this->get_P()).norm()<FLT_EPSILON*currentNorm*newNorm);
//                            VerboseDislocationNode(4,"  sessileNeighbor "<<std::get<1>(pair.second)->tag()<< " movable?"<<sessileNeighborMovable<<std::endl;);
//                            isMovable=(isMovable&&sessileNeighborMovable);
//                            //                                isMovable*=sessileNeighborMovable;
//                            if(!isMovable)
//                            {
//                                break;
//                            }
//                        }
//                    }
//                }
//            }
//
//            if(isMovable && this->isOnBoundary())
//            {
//                isMovable*=this->boundingBoxSegments().contains(X);
//            }
//            return isMovable;
//        }
        
        
    }
    
    template <int dim, short unsigned int corder>
    const Simplex<dim,dim>* DislocationNode<dim,corder>::includingSimplex() const
    {/*!\returns A pointer to the const Simplex imcluding *this PlanarDislocationNode
      */
        return p_Simplex;
    }
    
    template <int dim, short unsigned int corder>
    bool DislocationNode<dim,corder>::isVirtualBoundaryNode() const
    {
        return masterNode && this->meshFaces().size()==0;
//        return false;
    }
    
    // template <int dim, short unsigned int corder>
    // bool DislocationNode<dim,corder>::isBoundaryNode() const
    // {
    //     return this->isOnExternalBoundary();
    // }

    template <int dim, short unsigned int corder>
    bool DislocationNode<dim,corder>::isBoundaryNode() const
    {
        return isOnExternalBoundary();
    }
    
    template <int dim, short unsigned int corder>
    bool DislocationNode<dim,corder>::isGrainBoundaryNode() const
    {
        return this->isOnInternalBoundary();
    }

    template <int dim, short unsigned int corder>
    typename DislocationNode<dim,corder>::GlidePlaneContainerType DislocationNode<dim,corder>::glidePlanes() const
    {
        GlidePlaneContainerType temp;
        for (const auto &ln : this->loopNodes())
        {
            if (ln->periodicPlanePatch())
            {
                temp.emplace(ln->periodicPlanePatch()->glidePlane.get());
            }
        }
        return temp;
    }

    template <int dim, short unsigned int corder>
    typename DislocationNode<dim,corder>::PlanarMeshFaceContainerType DislocationNode<dim,corder>::meshFaces() const
    {
        PlanarMeshFaceContainerType temp;
        for (const auto &ln : this->loopNodes())
        {
            if (ln->periodicPlaneEdge.first)
            {
                temp = ln->periodicPlaneEdge.first->meshIntersection->faces;
            }
            if (ln->periodicPlaneEdge.second)
            {
                for (const auto &face : ln->periodicPlaneEdge.second->meshIntersection->faces)
                {
                    temp.emplace(face);
                }
            }
        }
        return temp;
    }

    template <int dim, short unsigned int corder>
    bool DislocationNode<dim,corder>::isOnExternalBoundary() const
    { /*!\returns _isOnExternalBoundarySegment.
       */
        bool _isOnExternalBoundary(false);
        for (const auto &face : meshFaces())
        {
            if (face->regionIDs.first == face->regionIDs.second)
            {
                _isOnExternalBoundary = true;
            }
        }

        return _isOnExternalBoundary;
    }

    template <int dim, short unsigned int corder>
    bool DislocationNode<dim,corder>::isOnInternalBoundary() const
    {
        bool _isOnInternalBoundary(false);
        for (const auto &face : meshFaces())
        {
            if (face->regionIDs.first != face->regionIDs.second)
            {
                _isOnInternalBoundary = true;
            }
        }
        return _isOnInternalBoundary;
    }

    template <int dim, short unsigned int corder>
    bool DislocationNode<dim,corder>::isOnBoundary() const
    {
        return isOnExternalBoundary() || isOnInternalBoundary();
    }
    template <int dim, short unsigned int corder>
    typename DislocationNode<dim,corder>::VectorDim DislocationNode<dim,corder>::bndNormal() const
    {
        VectorDim _outNormal(VectorDim::Zero());
        for (const auto &face : meshFaces())
        {
            _outNormal += face->outNormal();
        }
        const double _outNormalNorm(_outNormal.norm());
        if (_outNormalNorm > FLT_EPSILON)
        {
            _outNormal /= _outNormalNorm;
        }
        else
        {
            _outNormal.setZero();
        }
        return _outNormal;
    }
    template <int dim, short unsigned int corder>
    typename DislocationNode<dim, corder>::VectorDim DislocationNode<dim, corder>::snapToGlidePlanesinPeriodic(const VectorDim &P)
    {
        GlidePlaneContainerType gps(glidePlanes());
        switch (gps.size())
        {
        case 0:
        {
            assert(false && "Glide plane size must be larger than 0");
            return P;
            break;
        }
        case 1:
        {
            return (*gps.begin())->snapToPlane(P);
            break;
        }
        case 2:
        {
            const PlanePlaneIntersection<dim> ppi(**gps.begin(), **gps.rbegin());
            assert(ppi.type == PlanePlaneIntersection<dim>::INCIDENT && "Intersection must be incident");
            return ppi.P + (P - ppi.P).dot(ppi.d) * ppi.d;
            break;
        }
        case 3:
        {
            const PlanePlaneIntersection<dim> ppi(**gps.begin(), **gps.rbegin());
            assert(ppi.type == PlanePlaneIntersection<dim>::INCIDENT && "Intersection must be incident");
            const auto iterP(++gps.begin());
            const PlaneLineIntersection<dim> pli((*iterP)->P, (*iterP)->unitNormal, ppi.P, ppi.d);
            assert(pli.type == PlaneLineIntersection<dim>::INCIDENT && "Plane line intersection must be incident");
            return pli.P;
            break;
        }
        default:
        {
            const auto iterPlane1(gps.begin());
            const auto iterPlane2(++gps.begin());
            auto iterPlane3(++(++gps.begin()));

            const PlanePlaneIntersection<dim> ppi(**iterPlane1, **iterPlane2);
            assert(ppi.type == PlanePlaneIntersection<dim>::INCIDENT && "Intersection must be incident");
            const PlaneLineIntersection<dim> pli((*iterPlane3)->P, (*iterPlane3)->unitNormal, ppi.P, ppi.d);
            const VectorDim snappedPos(pli.P);
            assert(pli.type == PlaneLineIntersection<dim>::INCIDENT && "Plane line intersection must be incident");
            // if (pli.type!=PlaneLineIntersection<dim>::INCIDENT)
            // {
            //     std::cout<<" IterPlane1 "<<(*iterPlane1)->P.transpose()<<"\t"<<(*iterPlane1)->unitNormal.transpose()<<std::endl;
            //     std::cout<<" IterPlane2 "<<(*iterPlane2)->P.transpose()<<"\t"<<(*iterPlane2)->unitNormal.transpose()<<std::endl;
            //     std::cout<<" IterPlane3 "<<(*iterPlane3)->P.transpose()<<"\t"<<(*iterPlane3)->unitNormal.transpose()<<std::endl;
            //     std::cout<<"ppi infor "<<std::endl;
            //     std::cout<<ppi.P.transpose()<<"\t"<<ppi.d.transpose()<<std::endl;
            //     std::cout<<" glide plane size "<<gps.size()<<" Printing glide planes "<<std::endl;
            //     for (const auto& gp : gps)
            //     {
            //         std::cout<<gp->P.transpose()<<"\t"<<gp->unitNormal.transpose()<<"\t"<<gp->planeIndex<<std::endl;
            //     }
            //     std::cout<<" intersection type "<<pli.type<<std::endl;
            // }
            while (++iterPlane3 != gps.end())
            {
                if (!(*iterPlane3)->contains(snappedPos))
                {
                    std::cout << " Glide Plane " << (*iterPlane3)->P.transpose() << " " << (*iterPlane3)->unitNormal.transpose() << std::endl;
                    assert(false && "Plane must contain the position for glide plane size >=3");
                }
            }
            return snappedPos;
            break;
        }
        }
    }
    template <int dim, short unsigned int corder>
    typename DislocationNode<dim, corder>::GrainContainerType DislocationNode<dim, corder>::grains() const
    {
        GrainContainerType temp;
        for (const auto &glidePlane : glidePlanes())
        {
            temp.insert(&glidePlane->grain);
        }
        return temp;
    }

   template class DislocationNode<3,0>;
    
}
#endif
