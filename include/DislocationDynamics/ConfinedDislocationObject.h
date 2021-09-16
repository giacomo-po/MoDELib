/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2017 by Yinan Cui <cuiyinan@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_ConfinedDislocationObject_H_
#define model_ConfinedDislocationObject_H_

#include <set>
#include <map>
#include <vector>
#include <algorithm>
#include <Eigen/Dense>
#include <GlidePlane.h>
#include <Grain.h>
#include <PlanarMeshFace.h>
#include <FiniteLineSegment.h>
#include <LineLineIntersection.h>

//#include <BoundingFiniteLineSegments.h>

namespace model
{
    
    
    
    template <int dim>
    struct ConfinedDislocationObject : private std::set<const GlidePlane<dim>*>
    /*                              */,private std::set<const PlanarMeshFace<dim>*>
    /*                              */,public BoundingMeshSegments<dim>
    {
        
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef Grain<dim> GrainType;
        typedef std::set<const GrainType*> GrainContainerType;
        typedef GlidePlane<dim> GlidePlaneType;
        typedef std::set<const GlidePlaneType*> GlidePlaneContainerType;
        typedef PlanarMeshFace<dim> PlanarMeshFaceType;
        typedef std::set<const PlanarMeshFaceType*> PlanarMeshFaceContainerType;
        typedef MeshBoundarySegment<dim> MeshBoundarySegmentType;
        typedef std::vector<VectorDim> PositionCointainerType;
        
        
    private:
        
        /**********************************************************************/
        void updateBoundingBoxWithMeshFace(const PlanarMeshFace<dim>& face);

        
        
        
        //        GlidePlaneObserver<dim>& gpObserver;
        PositionCointainerType posCointainer;
        std::unique_ptr<FiniteLineSegment<dim>> _glidePlaneIntersections;
        bool _isOnExternalBoundary;
        bool _isOnInternalBoundary;
        VectorDim _outNormal;
        
        
    public:
        
        /**********************************************************************/
        ConfinedDislocationObject(const PositionCointainerType& temp) ;
        /**********************************************************************/
        ConfinedDislocationObject(const VectorDim& P0) ;
        
        /**********************************************************************/
        ConfinedDislocationObject(const VectorDim& P0,const VectorDim& P1) ;
        
        
        /**********************************************************************/
        ConfinedDislocationObject(const ConfinedDislocationObject<dim>& A,
                                  const ConfinedDislocationObject<dim>& B) ;
        
        ConfinedDislocationObject& confinedObject();

        const ConfinedDislocationObject& confinedObject() const;
        
        /**********************************************************************/
        void clear();
        
        /**********************************************************************/
        const GlidePlaneContainerType& glidePlanes() const;
        
        GlidePlaneContainerType& glidePlanes();
        
        /**********************************************************************/
        const PlanarMeshFaceContainerType& meshFaces() const;
        PlanarMeshFaceContainerType& meshFaces();
        
        /**********************************************************************/
        std::set<size_t> meshFaceIDs() const;
        /**********************************************************************/
        bool isOnMeshFaces(const std::set<size_t>& faceIDs) const;
        
        /**********************************************************************/
        VectorDim snapToGlidePlanes(const VectorDim& P);

        /**********************************************************************/
        VectorDim snapToGlidePlanesinPeriodic(const VectorDim& P);
        
        /**********************************************************************/
        const std::unique_ptr<FiniteLineSegment<dim>>& glidePlaneIntersections() const;
        
        /**********************************************************************/
        const bool& isOnExternalBoundary() const;
        
        /**********************************************************************/
        const bool& isOnInternalBoundary() const;
        
        /**********************************************************************/
        bool isOnBoundary() const;
        
        /**********************************************************************/
        const VectorDim& bndNormal() const;
        
        /**********************************************************************/
        void updateConfinement(const VectorDim& P0);
        
        /**********************************************************************/
        void updateConfinement(const VectorDim& P0,const VectorDim& P1);
        
        void updateConfinement(const PositionCointainerType& temp);
        
        /**********************************************************************/
        void addGlidePlane(const GlidePlaneType* const lastGlidePlane);
        /**********************************************************************/
        void updateConfinement();
        
        /**********************************************************************/
        GrainContainerType grains() const;
        
        /**********************************************************************/
        std::vector<std::pair<const GlidePlane<dim>* const,const GlidePlane<dim>* const>> parallelAndCoincidentGlidePlanes(const GlidePlaneContainerType& other) const;
        
        
    };
}
#endif