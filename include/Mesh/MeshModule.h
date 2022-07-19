/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_MeshModule_H_
#define model_MeshModule_H_

#include <cfloat>
#include <Eigen/Dense>


namespace model
{
    template<short int dim,short int dmo>
    struct BoundarySimplex;
    
    class GmshReader;
    
    template <int dim>
    class MeshBoundarySegment;
    
    template <int dim>
    struct BoundingMeshSegments;
    
    template <int dim>
    struct MeshPlane;
    
    template<int dim>
    struct MeshRegion;
    
    template<int dim>
    struct MeshRegionBoundary;
    
    template<typename RegionType>
    struct MeshRegionObserver;

    template<int dim>
    struct PlanarMeshFace;
    
    template<short int dim, short int order>
    class Simplex;
    
    template<short int _dim, short int order>
    struct SimplexBase;
    
    template<short int dim,short int order>
    class SimplexChild;
    
    template<short int dim>
    struct SimplexObserver;
    
    template<short int dim,short int order>
    struct SimplexParent;
    
    template<int dim>
    struct SimplexReaderBase;
    
    template<int dim>
    struct SimplexReader;
    
    template<short int nVertices>
    struct SimplexID;
    
    template<short int dim, short int order>
    struct SimplexTraits;
    
    template<short int dim, short int order>
    struct SimplexVolume;

    template<int _dim>
    class SimplicialMesh;
    
    class TriangularMesh;

    template<short int _dim>
    struct BarycentricTraits;
    
}

#include <BoundarySimplex.h>
#include <GmshReader.h>
#include <MeshBoundarySegment.h>
#include <MeshPlane.h>
#include <MeshRegion.h>
#include <MeshRegionBoundary.h>
#include <MeshRegionObserver.h>
#include <PlanarMeshFace.h>
#include <Simplex.h>
#include <SimplexBase.h>
#include <SimplexChild.h>
#include <SimplexObserver.h>
#include <SimplexParent.h>
#include <SimplexReader.h>
#include <SimplexTraits.h>
#include <SimplexVolume.h>
#include <SimplicialMesh.h>
#include <TriangularMesh.h>
#include <BarycentricTraits.h>


#endif
