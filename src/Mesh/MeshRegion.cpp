/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_MeshRegion_cpp_
#define model_MeshRegion_cpp_

#include <set>
#include <deque>
#include <memory>
#include <assert.h>

#include <MeshModule.h>

namespace model
{
    
    
    /**********************************************************************/
    template<int dim>
    void MeshRegion<dim>::buildSingleFace(const Simplex<dim,dim-1>* newStart,
                                          std::shared_ptr<PlanarMeshFace<dim>> newFace,
                                          std::map<typename Simplex<dim,dim-1>::SimplexIDType,const Simplex<dim,dim-1>*>& allSimplices)
    {
        allSimplices.erase(newStart->xID);
        const auto boundaryNeighbors(newStart->boundaryNeighbors());
        for(const auto& neighbor : boundaryNeighbors)
        {
            const auto neighborRgnIDs(neighbor->regionIDs());
            if(   neighborRgnIDs.find(regionID)!=neighborRgnIDs.end()  // neighbor in region
               && allSimplices.find(neighbor->xID)!=allSimplices.end())   // neighbor found in allSimplices
            {
                if((newStart->outNormal()-neighbor->outNormal()).norm()<FLT_EPSILON)
                {// same plane
                    allSimplices.erase(neighbor->xID);
                    newFace->insert(neighbor);
                    buildSingleFace(neighbor,newFace,allSimplices);
                }
            }
        }
        faces().emplace(newFace->sID,newFace);
    }
    
    /**********************************************************************/
    template<int dim>
    void MeshRegion<dim>::buildFaces()
    {/*!Constructs and stores the external PlanarMeshFace(s) of this MeshRegion
      * This function supports non-convex regions.
      */
        const auto t0= std::chrono::system_clock::now();
        std::cout<<"MeshRegion "<<regionID<<" buiding PlanarMeshFace's "<<std::flush;
        
        
        faces().clear();
        std::map<typename Simplex<dim,dim-1>::SimplexIDType,const Simplex<dim,dim-1>*> allSimplices;
        for(const auto& simplex : simplices())
        {// collect each possible boundary/region_boundary simplices
            for(const auto& child : simplex->children())
            {
                if(child->isBoundarySimplex())
                {
                    allSimplices.emplace(child->xID,child.get());
                }
            }
        }
        
        while(allSimplices.size())
        {
            std::shared_ptr<PlanarMeshFace<dim>> newFace(new PlanarMeshFace<dim>(allSimplices.begin()->second));
            buildSingleFace(allSimplices.begin()->second,newFace,allSimplices);
        }
        
        for(auto& face : faces())
        {
            face.second->finalize();
        }
        std::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
    }
    
    
    /**********************************************************************/
    template<int dim>
    MeshRegion<dim>::MeshRegion(MeshRegionObserverType& ro,
                                const size_t& rID) :
    /* init */ regionObserver(ro),
    /* init */ regionID(rID)
    {
        const bool success=regionObserver.emplace(regionID,this).second;
        if(!success)
        {
            throw std::runtime_error("Cannot insert MeshRegion in MeshRegionObserver.");
        }
    }
    
    /**********************************************************************/
    template<int dim>
    MeshRegion<dim>::~MeshRegion()
    {
        regionObserver.erase(regionID);
//        const size_t n=regionObserver.erase(regionID);
//        if(n!=1)
//        {
//            throw std::runtime_error("Cannot erase MeshRegion from MeshRegionObserver.");
//        }
    }
    
    /**********************************************************************/
    template<int dim>
    void MeshRegion<dim>::update()
    {
        buildFaces();
    }
    
    /**********************************************************************/
    template<int dim>
    void MeshRegion<dim>::identifyParallelFaces(const std::set<int>& periodicFaceIDs)
    {
        for(const auto& face1 : faces())
        {
            for(const auto& face2 : faces())
            {
                if(face1.second.get()!=face2.second.get())
                {
                    if(abs(face1.second->outNormal().dot(face2.second->outNormal())+1.0)<FLT_EPSILON)
                    {
                        _parallelFaces.emplace(face1.first,face2.first);
                    }
                }
            }
        }
        
        for(const auto& pair : _parallelFaces)
        {
            if(periodicFaceIDs.find(pair.first)!=periodicFaceIDs.end() || periodicFaceIDs.find(pair.second)!=periodicFaceIDs.end())
            {// either pair.first or pair.second are declared as periodic
                const VectorDim shift(faces()[pair.first]->center()-faces()[pair.second]->center());
                std::set<const Simplex<dim,0>*> secondFaceVertices;
                
                for(const auto& secondFaceVertex : faces()[pair.second]->convexHull())
                {
                    secondFaceVertices.insert(secondFaceVertex);
                }
                
                for(const auto& firstFaceVertex : faces()[pair.first]->convexHull())
                {
//                    bool shiftFound(false);
                    for(const auto& secondFaceVertex : secondFaceVertices)
                    {
                        if((firstFaceVertex->P0-secondFaceVertex->P0-shift).norm()<FLT_EPSILON)
                        {
                            secondFaceVertices.erase(secondFaceVertex);
                            break;
                        }
                    }
                }
                
                if(secondFaceVertices.size()==0)
                {//all vertices have been paired
                    std::cout<<"Detected periodic face pair "<<pair.first<<"-"<<pair.second<<std::endl;
                    faces()[pair.first] ->periodicFacePair=std::make_pair(-shift,faces()[pair.second].get());
                    faces()[pair.second]->periodicFacePair=std::make_pair( shift,faces()[pair.first ].get());
                }

            }
        }
    }
    
    /**********************************************************************/
    template<int dim>
    const std::map<size_t,size_t>& MeshRegion<dim>::parallelFaces() const
    {
        return _parallelFaces;
    }
    
    /**********************************************************************/
    template<int dim>
    const std::set<const typename MeshRegion<dim>::SimplexType*>& MeshRegion<dim>::simplices() const
    {
        return *this;
    }
    
    /**********************************************************************/
    template<int dim>
    std::set<const typename MeshRegion<dim>::SimplexType*>& MeshRegion<dim>::simplices()
    {
        return *this;
    }
    
    /**********************************************************************/
    template<int dim>
    const typename MeshRegion<dim>::MeshFacesContainerType& MeshRegion<dim>::faces() const
    {
        return *this;
    }
    
    /**********************************************************************/
    template<int dim>
    typename MeshRegion<dim>::MeshFacesContainerType& MeshRegion<dim>::faces()
    {
        return *this;
    }
    
    /**********************************************************************/
    template<int dim>
    typename MeshRegion<dim>::VectorDim MeshRegion<dim>::outNormal(const std::set<size_t>& faceIDs) const
    {
        typename MeshRegion<dim>::VectorDim temp(MeshRegion<dim>::VectorDim::Zero());
        for(const auto& val : faceIDs)
        {
            temp+=faces().at(val)->outNormal();
        }
        const double tempNorm(temp.norm());
        return tempNorm>FLT_EPSILON? (temp/tempNorm).eval() :  MeshRegion<dim>::VectorDim::Zero();
    }
    

//    template class MeshRegion<1>;
//    template class MeshRegion<2>;
    template struct MeshRegion<3>;
    
    
}	// close namespace
#endif
