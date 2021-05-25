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
#include <GramSchmidt.h>


namespace model
{
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    DislocationNode<dim,corder,InterpolationType>::DislocationNode(LoopNetworkType* const net,
                                                                   const VectorDim& P,
                                                                   const VectorDim& V,
                                                                   const double& vrc) :
    /* init */ NetworkNode<DislocationNode>(net)
    /* init */,SplineNodeType(P)
    /* init */,ConfinedDislocationObjectType(this->get_P())
    /* init */,p_Simplex(get_includingSimplex(this->get_P(),(const Simplex<dim,dim>*) NULL))
    /* init */,velocity(V)
    /* init */,vOld(velocity)
    /* init */,velocityReductionCoeff(vrc)
    /* init */,virtualNode(nullptr)
    /* init */,masterNode(nullptr)
    {

    }
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    std::shared_ptr<DislocationNode<dim,corder,InterpolationType>> DislocationNode<dim,corder,InterpolationType>::clone() const
    {
        return std::shared_ptr<DislocationNode<dim,corder,InterpolationType>>(new DislocationNode(this->p_network(),this->get_P(),get_V(),velocityReduction()));
    }
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    const Simplex<dim,dim>* DislocationNode<dim,corder,InterpolationType>::get_includingSimplex(const VectorDim& X,const Simplex<dim,dim>* const guess) const
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
                model::cout<<"PlanarDislocationNode "<<this->sID<<" @ "<<X.transpose()<<std::endl;
                std::cout<<"Simplex "<<temp.second->xID<<std::endl;
                model::cout<<"bary "<<temp.second->pos2bary(X)<<std::endl;
                std::cout<<"face of barymin is "<<temp.second->child(faceID).xID<<std::endl;
                model::cout<<"face of barymin is boundary Simplex? "<<temp.second->child(faceID).isBoundarySimplex()<<std::endl;
                model::cout<<"face of barymin is region-boundary Simplex? "<<temp.second->child(faceID).isRegionBoundarySimplex()<<std::endl;
                assert(0 && "DISLOCATION NODE OUTSIDE MESH.");
            }
        }
        return temp.second;
    }
    
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    void DislocationNode<dim,corder,InterpolationType>::projectVelocity()
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
            
            GramSchmidt::orthoNormalize(temp);
            
            for(const auto& vec : temp)
            {
                velocity-=velocity.dot(vec)*vec;
            }
            
        }
    }
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    void DislocationNode<dim,corder,InterpolationType>::set_V(const VectorDim& vNew)
    {
        vOld=velocity; // store current value of velocity before updating
        
        //            velocity=this->prjM*vNew; // kill numerical errors from the iterative solver
        velocity=vNew;
        projectVelocity();
        
        if(use_velocityFilter && !isBoundaryNode())
        {
            const double filterThreshold=0.05*velocity.norm()*vOld.norm()+FLT_EPSILON;
            
            if(velocity.dot(vOld)<-filterThreshold)
            {
                velocityReductionCoeff*=velocityReductionFactor;
            }
            else if(velocity.dot(vOld)>filterThreshold)
            {
                velocityReductionCoeff/=velocityReductionFactor;
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
    }
    


    template <int dim, short unsigned int corder, typename InterpolationType>
    void DislocationNode<dim,corder,InterpolationType>::initFromFile(const std::string& fileName)
    {
        use_velocityFilter=TextFileParser(fileName).readScalar<double>("use_velocityFilter",true);
        velocityReductionFactor=TextFileParser(fileName).readScalar<double>("velocityReductionFactor",true);
        assert(velocityReductionFactor>0.0 && velocityReductionFactor<=1.0);
        verboseDislocationNode=TextFileParser(fileName).readScalar<int>("verboseDislocationNode",true);
    }

    template <int dim, short unsigned int corder, typename InterpolationType>
    void DislocationNode<dim,corder,InterpolationType>::updateGeometry()
    {
        for(const auto& neighboor : this->neighbors())
        {
            std::get<1>(neighboor.second)->updateGeometry();
        }
    }

    template <int dim, short unsigned int corder, typename InterpolationType>
    void DislocationNode<dim,corder,InterpolationType>::trySet_P(const typename DislocationNode<dim,corder,InterpolationType>::VectorDim& newP)
    {
        if((this->get_P()-newP).norm()>FLT_EPSILON)
        {
                        for(auto& loopNode : this->loopNodes())
                        {
                            loopNode->set_P(loopNode->periodicPlanePatch()? newP-loopNode->periodicPlanePatch()->shift : newP);
                        }
        }
    }
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    bool DislocationNode<dim,corder,InterpolationType>::set_P(const typename DislocationNode<dim,corder,InterpolationType>::VectorDim& newP)
    {
        if((this->get_P()-newP).norm()>FLT_EPSILON)
        {
            SplineNodeType::set_P(newP);
            this->updateConfinement(this->get_P());
            p_Simplex=get_includingSimplex(this->get_P(),p_Simplex); // update including simplex
        }
        
        return (this->get_P()-newP).norm()<FLT_EPSILON;
    }
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    const typename DislocationNode<dim,corder,InterpolationType>::VectorDim& DislocationNode<dim,corder,InterpolationType>::get_V() const
    {/*! The nodal velocity vector
      */
        return velocity;
    }

    template <int dim, short unsigned int corder, typename InterpolationType>
    const double& DislocationNode<dim,corder,InterpolationType>::velocityReduction() const
    {
        return velocityReductionCoeff;
    }
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    typename DislocationNode<dim,corder,InterpolationType>::MeshLocation DislocationNode<dim,corder,InterpolationType>::meshLocation() const
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
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    const std::shared_ptr<typename DislocationNode<dim,corder,InterpolationType>::NetworkNodeType>& DislocationNode<dim,corder,InterpolationType>::virtualBoundaryNode() const
    {
        return virtualNode;
    }
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    typename DislocationNode<dim,corder,InterpolationType>::VectorDim DislocationNode<dim,corder,InterpolationType>::invariantDirectionOfMotion() const
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
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    bool DislocationNode<dim,corder,InterpolationType>::isMovableTo(const VectorDim& X) const
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
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    const Simplex<dim,dim>* DislocationNode<dim,corder,InterpolationType>::includingSimplex() const
    {/*!\returns A pointer to the const Simplex imcluding *this PlanarDislocationNode
      */
        return p_Simplex;
    }
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    bool DislocationNode<dim,corder,InterpolationType>::isVirtualBoundaryNode() const
    {
        return masterNode && this->meshFaces().size()==0;
//        return false;
    }
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    bool DislocationNode<dim,corder,InterpolationType>::isBoundaryNode() const
    {
        return this->isOnExternalBoundary();
    }
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    bool DislocationNode<dim,corder,InterpolationType>::isGrainBoundaryNode() const
    {
        return this->isOnInternalBoundary();
    }
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    bool DislocationNode<dim,corder,InterpolationType>::isRemovable(const double& Lmin,const double& relAreaThIn)
    {
        const double relAreaTh(this->loopNodes().size()>1? FLT_EPSILON : relAreaThIn);
        VerboseDislocationNode(2,"  Trying to remove DislocationNode "<<this->tag()<< ", relAreaThIn= "<<relAreaThIn<<std::endl;);
        VerboseDislocationNode(2,"  Trying to remove DislocationNode "<<this->tag()<< ", relAreaTh= "<<relAreaTh<<std::endl;);
        for(const auto& loopNode : this->loopNodes())
        {
            if(!loopNode->isRemovable(Lmin,relAreaTh))
            {
                VerboseDislocationNode(2,"  DislocationNode "<<this->tag()<< " NOT removable "<<std::endl;);
                return false;
            }
        }
        VerboseDislocationNode(2,"  DislocationNode "<<this->tag()<< " removable "<<std::endl;);
        return true;
    }

    
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    int DislocationNode<dim,corder,InterpolationType>::verboseDislocationNode=0;

    template <int dim, short unsigned int corder, typename InterpolationType>
    bool DislocationNode<dim,corder,InterpolationType>::use_velocityFilter=true;

    template <int dim, short unsigned int corder, typename InterpolationType>
    double DislocationNode<dim,corder,InterpolationType>::velocityReductionFactor=0.75;
    
}
#endif
