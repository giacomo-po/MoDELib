/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_MeshRegionBoundary_H_
#define model_MeshRegionBoundary_H_

#include <set>
#include <deque>
#include <assert.h>

#include <MeshRegionObserver.h>
#include <PlanarMeshFace.h>

namespace model
{

    template<int dim>
    struct MeshRegionBoundary : public std::set<const Simplex<dim,dim-1>*>
    /*                       */,public std::map<int,std::shared_ptr<PlanarMeshFace<dim>>>
    {
        typedef Simplex<dim,dim-1> SimplexType;
        typedef PlanarMeshFace<dim> PlanarMeshFaceType;
        typedef std::map<int,std::shared_ptr<PlanarMeshFaceType>> MeshFacesContainerType;
        
    private:
        
        void buildSingleFace(const Simplex<dim,dim-1>* newStart,
                             std::shared_ptr<PlanarMeshFace<dim>> newFace,
                             std::set<const Simplex<dim,dim-1>*>& allSimplices);
        void buildFaces();
                
    public:
        
        const std::pair<int,int> regionBndID;

        MeshRegionBoundary(const std::pair<size_t,size_t>& rbndID);
        void update();
        const std::set<const SimplexType*>& simplices() const;
        std::set<const SimplexType*>& simplices();
        const MeshFacesContainerType& faces() const;
        MeshFacesContainerType& faces();
        
    };
    
}	// close namespace
#endif
