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
#include <PlanesIntersection.h>

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
//    /* init */,virtualNode(nullptr)
//    /* init */,masterNode(nullptr)
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
                throw std::runtime_error("DISLOCATION NODE OUTSIDE MESH.");
//                assert(0 && "DISLOCATION NODE OUTSIDE MESH.");
            }
        }
        return temp.second;
    }

    
    template <int dim, short unsigned int corder>
    void DislocationNode<dim,corder>::projectVelocity()
    {
        
        VectorOfNormalsType temp;
//        temp.push_back(VectorDim::UnitZ());

        
        for(const auto& loopNode : this->loopNodes())
        {
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

            GramSchmidt::orthoNormalize(temp);
            
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
        
        velocity=vNew;
        projectVelocity();

        if (this->network().capMaxVelocity && !this->network().timeIntegrator.isTimeStepControlling(*this))
        {
            //Need to cap the max velocity based on the time step
            const double vmax=this->network().timeIntegrator.dxMax/this->network().timeIntegrator.dtMax;
            assert(fabs(vmax)>FLT_EPSILON && "Max velocity should be greater than 0");
            if (velocity.norm()>=FLT_EPSILON && velocity.norm()>fabs(vmax))
            {
                const double velredFac(velocity.norm()/fabs(vmax));
                assert(velredFac>=1.0 && "Velocity reduction factor must be greater than 1");
                velocity/=velredFac;
                totalCappedNodes+=1;
            }
        }
        
        if(this->network().use_velocityFilter)
        {
            if(!isBoundaryNode())
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
            else
            {
                // TO DO BOUNDARY NODES NEED TO INTERPOLATE ALSO THE VELOCTY FILTER
            }
            
        }
    }


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
        VerboseDislocationNode(1, " Try  Setting P for Network Node " << this->tag()<<" to"<<newP.transpose()<< std::endl;);
        const VectorDim snappedPosition(this->snapToGlidePlanesinPeriodic(newP));
        VerboseDislocationNode(2, " snappedPosition= " << snappedPosition.transpose()<< std::endl;);
        if((this->get_P()-snappedPosition).norm()>FLT_EPSILON)
        {
            for (auto &loopNode : this->loopNodes())
            {
                loopNode->set_P(loopNode->periodicPlanePatch() ? snappedPosition - loopNode->periodicPlanePatch()->shift : snappedPosition);
            }
        }

        return true;
    }
    
    template <int dim, short unsigned int corder>
    bool DislocationNode<dim,corder>::set_P(const typename DislocationNode<dim,corder>::VectorDim& newP)
    {
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

    template <int dim, short unsigned int corder>
    std::set<typename DislocationNode<dim,corder>::LoopType*> DislocationNode<dim,corder>::sessileLoops() const
    {
        std::set<LoopType*> temp;
        for (const auto& loop : this->loops())
        {
            if (loop->loopType == DislocationLoopIO<dim>::SESSILELOOP)
            {
                temp.insert(loop);
            }
        }
        return temp;
    }


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

//        if(!isVirtualBoundaryNode())
//        {
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
//        }
        return temp;
    }
    
//    template <int dim, short unsigned int corder>
//    const std::shared_ptr<typename DislocationNode<dim,corder>::NetworkNodeType>& DislocationNode<dim,corder>::virtualBoundaryNode() const
//    {
//        return virtualNode;
//    }
    
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
        
    }
    
    template <int dim, short unsigned int corder>
    const Simplex<dim,dim>* DislocationNode<dim,corder>::includingSimplex() const
    {/*!\returns A pointer to the const Simplex imcluding *this PlanarDislocationNode
      */
        return p_Simplex;
    }
    
//    template <int dim, short unsigned int corder>
//    bool DislocationNode<dim,corder>::isVirtualBoundaryNode() const
//    {
//        return masterNode && this->meshFaces().size()==0;
////        return false;
//    }
    


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
                break;
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
                break;
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
    typename DislocationNode<dim, corder>::VectorDim DislocationNode<dim, corder>::snapToGlidePlanesinPeriodic(const VectorDim &x) const
    {
        GlidePlaneContainerType gps(glidePlanes());
        if(gps.size())
        {
            Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> N(Eigen::Matrix<double,dim,Eigen::Dynamic>::Zero(dim,gps.size()));
            Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> P(Eigen::Matrix<double,dim,Eigen::Dynamic>::Zero(dim,gps.size()));

            int k=0;
            for(const auto& plane : gps)
            {
                N.col(k)=plane->unitNormal;
                P.col(k)=plane->P;
                ++k;
            }
            PlanesIntersection<dim> pInt(N,P,FLT_EPSILON);

            const std::pair<bool,VectorDim> snapped(pInt.snap(x));
            if(snapped.first)
            {
                return snapped.second;
            }
            else
            {
                throw std::runtime_error("Cannot snap, glidePlanes dont intersect.");
                return snapped.second;
            }
        }
        else
        {
            if(sessileLoops().size()==this->loops().size())
            {
                return x;
            }
            else
            {
                throw std::runtime_error("All loops must be sessile if there are no glide planes.");
                return x;
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
    
    template <int dim, short unsigned int corder>
    int DislocationNode<dim,corder>::totalCappedNodes=0;
    
    template class DislocationNode<3,0>;
    
}
#endif
