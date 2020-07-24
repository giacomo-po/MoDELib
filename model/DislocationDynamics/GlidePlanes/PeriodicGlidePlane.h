/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PeriodicGlidePlane_H_
#define model_PeriodicGlidePlane_H_


#include <memory>
#include <string>
#include <list>

#include <GlidePlane.h>
#include <GlidePlaneFactory.h>
//#include <PeriodicGlidePlaneFactory.h>
#include <TerminalColors.h>
#include <DiophantineSolver.h>
#include <SharedPtrFactories.h>



namespace model
{
    
    template<int dim>
    class PeriodicGlidePlaneFactory;
    
    template<int dim>
    class PeriodicPlanePatch;
    
    template<int dim>
    class PeriodicPlaneNode;
    
    template<int dim>
    class PeriodicGlidePlane;
    
//    template<int dim>
//    class PeriodicPlaneEdge;
    
    /**********************************************************************/
    /**********************************************************************/
    template<int dim>
    struct PeriodicPlaneEdge
    {
        
        const PeriodicPlanePatch<dim>* const patch;
        const std::shared_ptr<PeriodicPlaneNode<dim>> source;
        const std::shared_ptr<PeriodicPlaneNode<dim>>   sink;
        const std::shared_ptr<const MeshBoundarySegment<dim>> meshIntersection;
        PeriodicPlaneEdge<dim>* next;
        PeriodicPlaneEdge<dim>* prev;
        PeriodicPlaneEdge<dim>* twin;
        
        /**********************************************************************/
        PeriodicPlaneEdge(const PeriodicPlanePatch<dim>* const patch_in,
                          const std::shared_ptr<PeriodicPlaneNode<dim>>& source_in,
                          const std::shared_ptr<PeriodicPlaneNode<dim>>& sink_in,
                          const std::shared_ptr<const MeshBoundarySegment<dim>> meshIntersection_in)
        ;
        ~PeriodicPlaneEdge();
        
        std::string tag() const;
    };
    
    /**********************************************************************/
    /**********************************************************************/
    template<int dim>
    struct NodalConnectivity
    {
        typedef PeriodicPlaneEdge<dim> PeriodicPlaneEdgeType;
        
        PeriodicPlaneEdgeType* inEdge;
        PeriodicPlaneEdgeType* outEdge;
        
        NodalConnectivity() :
        /* init */ inEdge(nullptr)
        /* init */,outEdge(nullptr)
        {
            
        }
        
    };
    
    /**********************************************************************/
    /**********************************************************************/
    template<int dim>
    class PeriodicPlaneNode : public StaticID<PeriodicPlaneNode<dim>>
    /*                     */,public Eigen::Matrix<double,dim-1,1>
    {
        
    public:
        
        typedef Eigen::Matrix<double,dim-1,1> VectorLowerDim;
        typedef PeriodicPlaneNode<dim> PeriodicPlaneNodeType;
        typedef PeriodicPlaneEdge<dim> PeriodicPlaneEdgeType;
        typedef NodalConnectivity<dim> NodalConnectivityType;
        typedef std::map<size_t,NodalConnectivity<dim>> NodalConnectivityContainerType;
        typedef std::set<PeriodicPlaneEdgeType*> InOutEdgeContainerType;
        
    private:
        
        NodalConnectivityContainerType _patchConnectivities;
        NodalConnectivityContainerType _neighborConnectivities;
        
        
    public:
        
        PeriodicPlaneNode(const VectorLowerDim& pos);
        void addLink(PeriodicPlaneEdgeType* const link);
        void removeLink(PeriodicPlaneEdgeType* const link);
        const NodalConnectivityContainerType& patchConnectivities() const;
        NodalConnectivityContainerType& patchConnectivities();
        const NodalConnectivityContainerType& neighborConnectivities() const;
        NodalConnectivityContainerType& neighborConnectivities();
        InOutEdgeContainerType inEdges() const;
        InOutEdgeContainerType outEdges() const;
    };
    
    /**********************************************************************/
    /**********************************************************************/
    template<int dim>
    struct PeriodicPlanePatch : public StaticID<PeriodicPlanePatch<dim>>
    /*                      */, private std::vector<std::shared_ptr<PeriodicPlaneEdge<dim>>>
    {
        
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef VectorDim KeyType;
        typedef Eigen::Matrix<double,dim-1,1> VectorLowerDim;
        typedef CompareVectorsByComponent<double,dim,float> CompareType;

        PeriodicGlidePlane<dim>* const periodicPlane;
        const VectorDim shift;
        const std::shared_ptr<GlidePlane<dim>> glidePlane;
        typedef std::vector<std::shared_ptr<PeriodicPlaneEdge<dim>>> PeriodicPlaneEdgeContainerType;
        
        PeriodicPlanePatch(PeriodicGlidePlane<dim>& periodicPlane_in,
                           const VectorDim& shift_in);
        ~PeriodicPlanePatch();
        static bool isRightHandedBoundary(const BoundingMeshSegments<dim>& bnd,const Plane<dim>& plane);
        void addMeshIntersections(const BoundingMeshSegments<dim>& bms);
        const PeriodicPlaneEdgeContainerType& edges() const;
        int contains(const VectorLowerDim& test);
        
        
    };
    
    
    
    
    
    template<int dim>
    struct PeriodicGlidePlaneBase : public StaticID<PeriodicGlidePlaneBase<dim>>
    /*                           */,private std::map<Eigen::Matrix<double,dim-1,1>,const std::weak_ptr<PeriodicPlaneNode<dim>>,CompareVectorsByComponent<double,dim-1,float>>
    /*                           */,private std::set<const PeriodicPlaneEdge<dim>*>
    
    {
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef Eigen::Matrix<double,dim-1,1> VectorLowerDim;
        typedef std::map<Eigen::Matrix<double,dim-1,1>,const std::weak_ptr<PeriodicPlaneNode<dim>>,CompareVectorsByComponent<double,dim-1,float>> NodeCointainerType;
        typedef std::set<const PeriodicPlaneEdge<dim>*> UntwinnedEdgeContainerType;
        typedef std::vector<const PeriodicPlaneEdge<dim>*> BoundaryContainerType;
        typedef std::vector<BoundaryContainerType>    BoundariesContainerType;
        
        BoundariesContainerType _outerBoundaries;
        BoundariesContainerType _innerBoundaries;
        GlidePlaneFactory<dim>& glidePlaneFactory;
        const std::shared_ptr<GlidePlane<dim>> referencePlane;
        const MatrixDim L2G;

        
        static MatrixDim getL2G(VectorDim z);
        PeriodicGlidePlaneBase(GlidePlaneFactory<dim>& glidePlaneFactory_in,const GlidePlaneKey<dim>& referencePlaneKey);
        GlidePlaneKey<dim> getGlidePlaneKey(const VectorDim& shift);
        std::shared_ptr<GlidePlane<dim>> getGlidePlane(const VectorDim& shift);
        BoundariesContainerType& outerBoundaries();
        const BoundariesContainerType& outerBoundaries() const;
        BoundariesContainerType& innerBoundaries();
        const BoundariesContainerType& innerBoundaries() const;
        const UntwinnedEdgeContainerType& untwinnedEdges() const;
        UntwinnedEdgeContainerType& untwinnedEdges();
        void addUntwinnedEdge(const PeriodicPlaneEdge<dim>* link);
        void removeUntwinnedEdge(const PeriodicPlaneEdge<dim>* link);
        NodeCointainerType& nodes();
        const NodeCointainerType& nodes() const;
        int isInsideOuterBoundary(const VectorLowerDim& test);
        bool isCompact() const;
        static VectorDim rightHandedNormal(const std::vector<const PeriodicPlaneEdge<dim>*>& bnd);
        double rightHandedArea(const std::vector<const PeriodicPlaneEdge<dim>*>& bnd);
        bool isRightHandedBoundary(const std::vector<const PeriodicPlaneEdge<dim>*>& bnd);
        VectorLowerDim getLocalPosition(const VectorDim& point,const VectorDim& shift) const;
        VectorDim getGlobalPosition(const VectorLowerDim& point) const;
        std::shared_ptr<PeriodicPlaneNode<dim> > getSharedNode(const VectorDim& pointDim,const VectorDim& shift);
    };
    


    
    /**********************************************************************/
    /**********************************************************************/
    // ,
    template<int dim>
    class PeriodicGlidePlane : public PeriodicGlidePlaneBase<dim>
    /*                      */,public KeyConstructableSharedPtrFactory<PeriodicGlidePlane<dim>,PeriodicPlanePatch<dim>> // container of patches
    {
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef Eigen::Matrix<long int, dim, 1> VectorDimI;
        typedef Eigen::Matrix<double,dim-1,1> VectorLowerDim;
        typedef KeyConstructableSharedPtrFactory<PeriodicGlidePlane<dim>,PeriodicPlanePatch<dim>> PatchContainerType;
        typedef std::set<const PeriodicPlaneEdge<dim>*> UntwinnedEdgeContainerType;
        typedef std::vector<const PeriodicPlaneEdge<dim>*> BoundaryContainerType;
        typedef std::vector<BoundaryContainerType>    BoundariesContainerType;
        typedef PeriodicGlidePlaneFactory<dim> PeriodicGlidePlaneFactoryType;
        typedef GlidePlane<dim> GlidePlaneType;
        typedef typename GlidePlaneType::KeyType GlidePlaneKeyType;

    public:
        
        PeriodicGlidePlaneFactoryType& periodicGlidePlaneFactory;

        
        PeriodicGlidePlane(PeriodicGlidePlaneFactoryType& pgpf,const GlidePlaneKeyType& key_in);
        const PatchContainerType& patches() const;
        PatchContainerType& patches();
        void createNewBoundary(const PeriodicPlaneEdge<dim>* currentEdge,UntwinnedEdgeContainerType& untwinnedCopy);
        void updateBoundaries();
        VectorDim getShift(const PeriodicPlaneEdge<dim>& edge) const;
        void fillHoles();
        void addPatchesContainingPolygon(const std::vector<VectorDim>& polyPoints,const VectorDim& shift);
        void addPatchesContainingPolygon(const std::vector<VectorDim>& polyPoints);
        void addPatchesContainingPolygon(const std::vector<VectorLowerDim>& polyPoints);
        void print();
        std::shared_ptr<PeriodicPlanePatch<dim>> getPatch(const VectorDim& shift);
        VectorDim getGlidePlaneShiftfromReferencePlane(const GlidePlane<dim> *gp) const;
    };
    
    template<int dim>
    struct PeriodicPlanePatchIO
    {
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        size_t glidePlaneID;
        size_t referencePlaneID;
        VectorDim shift;
        
        PeriodicPlanePatchIO();
        PeriodicPlanePatchIO(const PeriodicPlanePatch<dim>& patch);
        PeriodicPlanePatchIO(std::stringstream& ss);
//        template <class T>
//        friend T& operator << (T& os, const PeriodicPlanePatchIO<dim>& ds);
        
    };
    

}
#endif
