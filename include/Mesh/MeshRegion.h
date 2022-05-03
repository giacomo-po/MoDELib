/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_MeshRegion_H_
#define model_MeshRegion_H_

#include <set>
#include <deque>
#include <memory>
#include <assert.h>

#include <MeshModule.h>
//#include <MeshRegionObserver.h>
//#include <PlanarMeshFace.h>
//#include <Simplex.h>
//#include <TerminalColors.h>

namespace model
{
    
    /**************************************************************************/
    /**************************************************************************/
    template<int dim>
    struct MeshRegion : public std::set<const Simplex<dim,dim>*>
    /*               */,public std::map<int,std::shared_ptr<PlanarMeshFace<dim>>> // MeshRegionBoundary container
    {
        typedef Simplex<dim,dim> SimplexType;
//        static constexpr int dim=SimplexType::dim;
        typedef MeshRegion<dim> MeshRegionType;
        typedef MeshRegionObserver<MeshRegionType> MeshRegionObserverType;
        typedef std::map<int,std::shared_ptr<PlanarMeshFace<dim>>> MeshFacesContainerType;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        
    private:
        
        void buildSingleFace(const Simplex<dim,dim-1>* newStart,
                        std::shared_ptr<PlanarMeshFace<dim>> newFace,
                             std::map<typename Simplex<dim,dim-1>::SimplexIDType,const Simplex<dim,dim-1>*>& allSimplices);
        void buildFaces();
        
        MeshRegionObserverType& regionObserver;
        std::map<size_t,size_t> _parallelFaces;
        
    public:
        const size_t regionID;
        
        
        MeshRegion(MeshRegionObserverType& ro,
                   const size_t& rID);
        ~MeshRegion();
        void update();
        void identifyParallelFaces(const std::set<int>&);
        const std::map<size_t,size_t>& parallelFaces() const;
        const std::set<const SimplexType*>& simplices() const;
        std::set<const SimplexType*>& simplices();
        const MeshFacesContainerType& faces() const;
        MeshFacesContainerType& faces();
        VectorDim outNormal(const std::set<size_t>& faceIDs) const;
    };
    
}	// close namespace
#endif
