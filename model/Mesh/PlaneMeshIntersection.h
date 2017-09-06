/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2017 by Yinan Cui <cuiyinan@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PlaneMeshIntersection_H_
#define model_PlaneMeshIntersection_H_

#include <model/Utilities/StaticID.h>
#include <model/Utilities/NonCopyable.h>
#include <model/Mesh/SimplexTraits.h>
#include <model/Mesh/SimplexObserver.h>
#include <model/Mesh/SimplicialMesh.h>

namespace model
{
    
    /**************************************************************************/
    /**************************************************************************/
    template<int dim>
    struct PlaneMeshIntersection
    {
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef std::deque<std::pair<const Simplex<dim,dim-2>* const,VectorDim>> PlaneMeshIntersectionContainerType;
        typedef Simplex<dim,dim-2>   EdgeSimplexType;
        typedef std::set<const EdgeSimplexType*,SimplexCompare<dim,dim-2> >  SiblingsContainerType;
        typedef Simplex<dim,dim-1> ParentSimplexType;
        typedef std::pair<VectorDim,const Simplex<dim,0>* const> RootType;
        typedef std::deque<RootType> RootContainerType;
        
        static constexpr double meshIntersectionTol=FLT_EPSILON;

        
        const SimplicialMesh<dim>& mesh;
        
        /**********************************************************************/
        PlaneMeshIntersection(const SimplicialMesh<dim>& m) :
        /* init */ mesh(m)
        {
        
        }
        
        /**********************************************************************/
        PlaneMeshIntersectionContainerType reducedPlaneMeshIntersection(const VectorDim& P0,
                                                                        const VectorDim& nn,
                                                                        const int& rID)
        {
            const PlaneMeshIntersectionContainerType mpi=planeMeshIntersection(P0,nn,rID);
            
            const auto t0=std::chrono::system_clock::now();
            model::cout<<"Reducing plane/mesh intersection points "<<std::flush;

            PlaneMeshIntersectionContainerType temp;
            
            size_t tailID=mpi.size()-1;
            for(size_t k=0;k<mpi.size();++k)
            {
                const size_t k1=(k==mpi.size()-1? 0 : k+1);
            
                VectorDim dX=mpi[k].second-mpi[tailID].second;
                VectorDim dXnew=mpi[k1].second-mpi[k].second;
                
                const double dXnorm=dX.norm();
                const double dXnewNorm=dXnew.norm();
                assert(dXnorm > FLT_EPSILON);
                assert(dXnewNorm > FLT_EPSILON);
                
                if(dX.cross(dXnew).norm()/(dXnorm*dXnewNorm)>FLT_EPSILON )
                {// non zero angle
                    tailID=k;
                    temp.push_back(mpi[k]);
                }
            
            }
            
            model::cout<<"("<<mpi.size()<<"->"<<temp.size()<<")";
            model::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;

            return temp;
        }
        
        /**********************************************************************/
        PlaneMeshIntersectionContainerType planeMeshIntersection(const VectorDim& P0,
                                                                        const VectorDim& nn,
                                                                        const int& rID)
        {
            model::cout<<"Computing plane/mesh intersection "<<std::flush;
            const auto t0=std::chrono::system_clock::now();
            
            const double nNorm(nn.norm());
            assert(nNorm>FLT_EPSILON);
            const VectorDim N(nn/nNorm);
            
            PlaneMeshIntersectionContainerType temp;
            
                std::set<const EdgeSimplexType*> tested;
                
                // Find initial intersection
                const std::pair<const EdgeSimplexType*,RootContainerType> pEI=getFirstPlaneEdgeIntersection(P0,N,rID,tested);
                
                switch (pEI.second.size())
                {
                    case 1:
                    {
                        if(pEI.second[0].second==nullptr)
                        {// intersection within edge
                            temp.emplace_back(pEI.first,pEI.second[0].first);
                            intersectionStep(P0,N,validSiblings(*pEI.first,rID),tested,temp,rID);
                        }
                        else
                        {// intersection at node
                            temp.emplace_back(pEI.first,pEI.second[0].first);
                            markParents(*pEI.second[0].second,tested);
                            intersectionStep(P0,N,validUncles(*pEI.second[0].second,rID),tested,temp,rID);
                        }
                        
                        break;
                    }
                        
                    case 2:
                    {
                        temp.emplace_back(pEI.first,pEI.second[0].first);
                        markParents(*pEI.second[0].second,tested);
                        
                        temp.emplace_back(pEI.first,pEI.second[1].first);
                        markParents(*pEI.second[1].second,tested);
                        intersectionStep(P0,N,validUncles(*pEI.second[1].second,rID),tested,temp,rID);
                        
                        break;
                    }
                        
                    default:
                        assert(0 && "IMPOSSIBLE");
                        break;
                }
            
            model::cout<<" ("<<temp.size()<<" perimeter points)";
            model::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
            
            return temp;
        }
        
        /**********************************************************************/
        void markParents(const Simplex<dim,0>& vertex,
                                std::set<const EdgeSimplexType*>& tested)
        {
            for(const auto& parent : vertex.parents())
            {
                tested.insert(parent);
            }
        }
        
        
        /**********************************************************************/
        SiblingsContainerType validSiblings(const Simplex<dim,dim-2>& edge,
                                                   const int& rID)
        {
            
            SiblingsContainerType siblings;
            for(const auto& parent : edge.parents())
            {
                if(  (parent->isBoundarySimplex() || parent->isRegionBoundarySimplex())
                   && parent->isInRegion(rID))
                {
                    for(int c=0; c<ParentSimplexType::nFaces;++c)
                    {
                        siblings.insert(&parent->child(c));
                    }
                }
            }
            
            return siblings;
        }
        
        
        /**********************************************************************/
        SiblingsContainerType validUncles(const Simplex<dim,0>& vertex,
                                                 const int& rID)
        {
            
            // Collect triangles attached to vertex
            std::set<const Simplex<dim,dim-1>*> grandParents;
            for(const auto& parent : vertex.parents())
            {
                for(const auto& grandParent : parent->parents())
                {
                    grandParents.insert(grandParent);
                }
            }
            
            // Loop over triangles and collect valid edges
            SiblingsContainerType uncles;
            for(const auto& grandParent : grandParents)
            {
                if(  (grandParent->isBoundarySimplex() || grandParent->isRegionBoundarySimplex())
                   && grandParent->isInRegion(rID))
                {
                    for(int c=0; c<Simplex<dim,dim-1>::nFaces;++c)
                    {
                        uncles.insert(&grandParent->child(c));
                    }
                }
            }
            
            return uncles;
        }
        
        /**********************************************************************/
        void intersectionStep(const VectorDim& P0,
                                     const VectorDim& N,
                                     const SiblingsContainerType& edgeContainer,
                                     std::set<const EdgeSimplexType*>& tested,
                                     PlaneMeshIntersectionContainerType& temp,
                                     const int& rID)
        {
            
            std::deque<std::pair<const EdgeSimplexType*,RootType>> rootDeq;
            for(const auto& edge : edgeContainer)
            {
                if(tested.find(edge)==tested.end())
                {
                    RootContainerType root=planeEdgeIntersection(P0,N,*edge,tested);
                    
                    switch (root.size())
                    {
                        case 0:
                            break;
                            
                        case 1:
                            rootDeq.emplace_back(edge,root[0]);
                            break;
                            
                        default:
                            assert(0 && "uncles cannot have more than one root");
                            break;
                    }
                }
            }
            
            
            
            if(rootDeq.size()>0)
            {
                if(rootDeq.size()>1 )
                {
                    if(temp.size()>1) // not start of loop
                    {
                        for(const auto& pair : rootDeq)
                        {
                            assert((pair.second.first-rootDeq[0].second.first).norm()<FLT_EPSILON); // all intersections are the same, and they must happen at a node
                        }
                        
                        assert(rootDeq[0].second.second!=nullptr);
                        
                        temp.emplace_back(rootDeq[0].first,rootDeq[0].second.first);
                        markParents(*rootDeq[0].second.second,tested);
                        intersectionStep(P0,N,validUncles(*rootDeq[0].second.second,rID),tested,temp,rID);
                        
                    }
                    else
                    {
                        if(rootDeq[0].second.second==nullptr)
                        {// intersection along edge
                            temp.emplace_back(rootDeq[0].first,rootDeq[0].second.first);
                            intersectionStep(P0,N,validSiblings(*rootDeq[0].first,rID),tested,temp,rID);
                        }
                        else
                        {// intersection at node
                            temp.emplace_back(rootDeq[0].first,rootDeq[0].second.first);
                            markParents(*rootDeq[0].second.second,tested);
                            intersectionStep(P0,N,validUncles(*rootDeq[0].second.second,rID),tested,temp,rID);
                            
                        }
                        
                    }
                }
                else
                {
                    temp.emplace_back(rootDeq[0].first,rootDeq[0].second.first);
                    intersectionStep(P0,N,validSiblings(*rootDeq[0].first,rID),tested,temp,rID);
                }
            }
            
        }
        
        
        
        /**********************************************************************/
        std::pair<const Simplex<dim,dim-2>*,RootContainerType> getFirstPlaneEdgeIntersection(const VectorDim& P0,
                                                                                                    const VectorDim& n,
                                                                                                    const int& rID,
                                                                                                    std::set<const EdgeSimplexType*>& tested)
        {
            std::pair<const Simplex<dim,dim-2>*,RootContainerType> temp;
            
            for(const auto& edge : mesh.template observer<dim-2>() )
            {// loop over edges
                if(edge.second->isBoundarySimplex() || edge.second->isRegionBoundarySimplex())
                {
                    if(edge.second->isInRegion(rID))
                    {
                        temp=std::make_pair(edge.second,planeEdgeIntersection(P0,n,*edge.second,tested));
                        
                        if(temp.second.size())
                        { // first intersection found
                            break;
                        }
                    }
                    
                }
            }
            
            assert(temp.second.size() && "plane does not intersect mesh.");
            
            return temp;
        }
        
        /**********************************************************************/
        RootContainerType planeEdgeIntersection(const VectorDim& P0,
                                                       const VectorDim& n,
                                                       const Simplex<dim,dim-2>& edge,
                                                       std::set<const EdgeSimplexType*>& tested)
        {
            //            std::deque<VectorDim> temp;
            tested.insert(&edge);
            
            RootContainerType temp;
            
            const VectorDim& v0(edge.child(0).P0);
            const VectorDim& v1(edge.child(1).P0);
            
            // check intersection of v0->v1 with plane
            // x=v0+u(v1-v0)
            // (x-P0).n=0
            // (v0+u(v1-v0)-P0).n=0
            // u=(P0-v0).n/(v1-v0).n;
            const double edgeNorm=(v1-v0).norm();
            assert(edgeNorm>FLT_EPSILON && "mesh edge has zero norm.");
            const double den=(v1-v0).dot(n);
            const double num=(P0-v0).dot(n);
            const double P0v0norm=(P0-v0).norm();
            
            const double numCheck= (P0v0norm<FLT_EPSILON)? 0.0 : num/P0v0norm;
            
            if (fabs(den/edgeNorm)>meshIntersectionTol)
            {
                // edge intersects plane
                const double u=num/den;
                
                if(fabs(u)<meshIntersectionTol)
                {
                    temp.emplace_back(v0,&edge.child(0));
                }
                else if (u>=meshIntersectionTol && u<=1.0-meshIntersectionTol)
                {
                    temp.emplace_back((1.0-u)*v0 + u*v1,nullptr);
                }
                else if (fabs(1.0-u)<meshIntersectionTol)
                {
                    temp.emplace_back(v1,&edge.child(1));
                }
                else
                {// no roots
                    
                }
                
            }
            else
            {
                if (fabs(numCheck)>meshIntersectionTol)
                {// edge is parallel to plane, no intersection
                    
                }
                else
                {// edge is coplanar
                    temp.emplace_back(v0,&edge.child(0));
                    temp.emplace_back(v1,&edge.child(1));
                }
            }
            
            return temp;
            //            return std::make_pair(&edge,temp);
        }
        
    };
    
}	// close namespace
#endif