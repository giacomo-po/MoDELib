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
#include <GramSchmidt.h>


namespace model
{
    
    template<int dim>
    struct PeriodicGlidePlaneFactory;
    
    template<int dim>
    struct PeriodicPlanePatch;
    
    template<int dim>
    class PeriodicPlaneNode;
   
    template<int dim>
    struct PeriodicPatchBoundary;
    
    template<int dim>
    class PeriodicGlidePlane;
    
//    template<int dim>
//    class PeriodicPlaneEdge;
    
    /**********************************************************************/
    /**********************************************************************/
    template<int dim>
    struct PeriodicPlaneEdge
    {
        
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        const PeriodicPlanePatch<dim>* const patch;
        const std::shared_ptr<PeriodicPlaneNode<dim>> source;
        const std::shared_ptr<PeriodicPlaneNode<dim>>   sink;
        const std::shared_ptr<const MeshBoundarySegment<dim>> meshIntersection;
        const short int edgeID;
        const VectorDim deltaShift;
        
        PeriodicPlaneEdge<dim>* next;
        PeriodicPlaneEdge<dim>* prev;
        PeriodicPlaneEdge<dim>* twin;
        
        //static VectorDim getDeltaShift(const PeriodicPlanePatch<dim>* const patch,const std::shared_ptr<const MeshBoundarySegment<dim>> meshIntersection);
        
        /**********************************************************************/
        PeriodicPlaneEdge(const PeriodicPlanePatch<dim>* const patch_in,
                          const std::shared_ptr<PeriodicPlaneNode<dim>>& source_in,
                          const std::shared_ptr<PeriodicPlaneNode<dim>>& sink_in,
                          const std::shared_ptr<const MeshBoundarySegment<dim>> meshIntersection_in,
                          const short int&)
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
        typedef VectorLowerDim KeyType;
        typedef CompareVectorsByComponent<double,dim-1,float> CompareType;
        
    private:
        
        NodalConnectivityContainerType _patchConnectivities;
        NodalConnectivityContainerType _neighborConnectivities;
        
        
    public:
        
        PeriodicPlaneNode(PeriodicPatchBoundary<dim>* const,const VectorLowerDim& pos);
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
        typedef std::vector<std::shared_ptr<PeriodicPlaneEdge<dim>>> PeriodicPlaneEdgeContainerType;

        PeriodicPatchBoundary<dim>* const patchBoundary;
        const VectorDim shift;
        const std::shared_ptr<GlidePlane<dim>> glidePlane;
        
        PeriodicPlanePatch(PeriodicPatchBoundary<dim>* const,const VectorDim&);
        ~PeriodicPlanePatch();
        static bool isRightHandedBoundary(const BoundingMeshSegments<dim>&, const Plane<dim>&);
        void addMeshIntersections(const BoundingMeshSegments<dim>&);
        const PeriodicPlaneEdgeContainerType& edges() const;
        std::shared_ptr<PeriodicPlaneEdge<dim>> edgeOnFaces(const std::set<const PlanarMeshFace<dim>*>& faces) const;
        int contains(const VectorLowerDim& test) const;
    };
    
    
    
    
    
    template<int dim>
    struct PeriodicPatchBoundary : public StaticID<PeriodicPatchBoundary<dim>>
//    /*                          */,private std::map<Eigen::Matrix<double,dim-1,1>,const std::weak_ptr<PeriodicPlaneNode<dim>>,CompareVectorsByComponent<double,dim-1,float>>
    /*                          */,public KeyConstructableWeakPtrFactory<PeriodicPatchBoundary<dim>,PeriodicPlaneNode<dim>,typename PeriodicPlaneNode<dim>::CompareType>
    /*                          */,private std::set<const PeriodicPlaneEdge<dim>*>
    /*                          */,public KeyConstructableWeakPtrFactory<PeriodicPatchBoundary<dim>,PeriodicPlanePatch<dim>,typename PeriodicPlanePatch<dim>::CompareType>
    
    {
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef Eigen::Matrix<double,dim-1,1> VectorLowerDim;
        typedef KeyConstructableWeakPtrFactory<PeriodicPatchBoundary<dim>,PeriodicPlaneNode<dim>,typename PeriodicPlaneNode<dim>::CompareType> NodeCointainerType;
        //        typedef std::map<Eigen::Matrix<double,dim-1,1>,const std::weak_ptr<PeriodicPlaneNode<dim>>,CompareVectorsByComponent<double,dim-1,float>> NodeCointainerType;
        typedef std::set<const PeriodicPlaneEdge<dim>*> UntwinnedEdgeContainerType;
        typedef std::vector<const PeriodicPlaneEdge<dim>*> BoundaryContainerType;
        typedef std::vector<BoundaryContainerType>    BoundariesContainerType;
        typedef KeyConstructableWeakPtrFactory<PeriodicPatchBoundary<dim>,PeriodicPlanePatch<dim>,typename PeriodicPlanePatch<dim>::CompareType> PatchContainerType;

        BoundariesContainerType _outerBoundaries;
        BoundariesContainerType _innerBoundaries;
        GlidePlaneFactory<dim>& glidePlaneFactory;
        const std::shared_ptr<GlidePlane<dim>> referencePlane;
//        const MatrixDim L2G;

        
//        static MatrixDim getL2G(VectorDim z);
        PeriodicPatchBoundary(GlidePlaneFactory<dim>&, const std::shared_ptr<GlidePlane<dim>>& referencePlane_in);
//        GlidePlaneKey<dim> getGlidePlaneKey(const VectorDim& shift);
//        std::shared_ptr<GlidePlane<dim>> getGlidePlane(const VectorDim& shift);
        void createNewBoundary(const PeriodicPlaneEdge<dim>* currentEdge,UntwinnedEdgeContainerType& untwinnedCopy);
        void updateBoundaries();
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
        std::shared_ptr<PeriodicPlanePatch<dim>> getPatch(const VectorDim& shift);
        GlidePlaneKey<dim> getGlidePlaneKey(const VectorDim& shift);
        std::shared_ptr<GlidePlane<dim>> getGlidePlane(const VectorDim& shift);
        const PatchContainerType& patches() const;
        PatchContainerType& patches();


    };
    


    
    /**********************************************************************/
    /**********************************************************************/
    // ,
    template<int dim>
    class PeriodicGlidePlane : public PeriodicPatchBoundary<dim>
//    /*                      */,public KeyConstructableSharedPtrFactory<PeriodicGlidePlane<dim>,PeriodicPlanePatch<dim>> // container of patches
    {


    public:
        
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef Eigen::Matrix<long int, dim, 1> VectorDimI;
        typedef Eigen::Matrix<double,dim-1,1> VectorLowerDim;
//        typedef KeyConstructableSharedPtrFactory<PeriodicGlidePlane<dim>,PeriodicPlanePatch<dim>> PatchContainerType;
        typedef std::set<const PeriodicPlaneEdge<dim>*> UntwinnedEdgeContainerType;
        typedef std::vector<const PeriodicPlaneEdge<dim>*> BoundaryContainerType;
        typedef std::vector<BoundaryContainerType>    BoundariesContainerType;
        typedef PeriodicGlidePlaneFactory<dim> PeriodicGlidePlaneFactoryType;
        typedef GlidePlane<dim> GlidePlaneType;
        typedef typename GlidePlaneType::KeyType KeyType;
        
        GlidePlaneFactory<dim>& glidePlaneFactory;
        PeriodicGlidePlaneFactoryType& periodicGlidePlaneFactory;
        const std::shared_ptr<GlidePlane<dim>> referencePlane;

        
        PeriodicGlidePlane(PeriodicGlidePlaneFactoryType* const pgpf,const KeyType& key_in);
        void print();
//        VectorDim getGlidePlaneShiftfromReferencePlane(const GlidePlane<dim> *gp) const;

        template<typename T>
        std::vector<std::tuple<VectorLowerDim,VectorDim,std::pair<short int,short  int>,std::map<VectorDim,std::set<std::pair<short int,short  int>>,CompareVectorsByComponent<double,dim,float>>,size_t ,const T* const>> polygonPatchIntersection(const std::vector<std::pair<VectorDim,const T* const>>& polyPoints,const bool& intersectInternalNodes=false);
        template<typename T>
        std::vector<std::tuple<VectorLowerDim,VectorDim,std::pair<short int,short  int>,std::map<VectorDim,std::set<std::pair<short int,short  int>>,CompareVectorsByComponent<double,dim,float>>,size_t ,const T* const>> polygonPatchIntersection(const std::vector<std::pair<VectorLowerDim,const T* const>>& polyPoints,const bool& intersectInternalNodes=false);
        
        
        VectorDim findPatch(const VectorLowerDim&,const VectorDim&);
        
        std::vector<std::shared_ptr<PeriodicPlanePatch<dim>>> filledPatches(const std::vector<VectorDim>& patchShifts);

    };
    
    
    
    template<int dim>
    struct PeriodicPlanePatchIO
    {
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        size_t glidePlaneID;
        size_t patchBoundaryID;
        VectorDim shift;
        
        PeriodicPlanePatchIO();
        PeriodicPlanePatchIO(const PeriodicPlanePatch<dim>& patch);
        PeriodicPlanePatchIO(std::stringstream& ss);
//        template <class T>
//        friend T& operator << (T& os, const PeriodicPlanePatchIO<dim>& ds);
        
    };

template <int dim,class T>
T& operator << (T& os, const PeriodicPlanePatchIO<dim>& ds)
{
    os  << ds.glidePlaneID<<"\t"
    /**/<< ds.patchBoundaryID<<"\t"
    /**/<< std::setprecision(15)<<std::scientific<<ds.shift.transpose();
    return os;
}

template<int dim>
template<typename T>
 std::vector<std::tuple<typename PeriodicGlidePlane<dim>::VectorLowerDim,typename PeriodicGlidePlane<dim>::VectorDim,std::pair<short int,short  int>,std::map<typename PeriodicGlidePlane<dim>::VectorDim,std::set<std::pair<short int,short  int>>,CompareVectorsByComponent<double,dim,float>>,size_t, const T* const>> PeriodicGlidePlane<dim>::polygonPatchIntersection(const std::vector<std::pair<VectorDim,const T* const>>& polyPoints, const bool& intersectInternalNodes)
{
    std::vector<std::pair<VectorLowerDim,const T* const>> lowerPolyPoints;
    for(const auto& point : polyPoints)
    {
        assert(this->referencePlane->contains(point.first) && "reference plane does not cointain point");
        lowerPolyPoints.push_back(std::make_pair(referencePlane->localPosition(point.first-VectorDim::Zero()),point.second));
    }
    return polygonPatchIntersection(lowerPolyPoints,intersectInternalNodes);
}

template<int dim>
template<typename T>
std::vector<std::tuple<typename PeriodicGlidePlane<dim>::VectorLowerDim,typename PeriodicGlidePlane<dim>::VectorDim,std::pair<short int,short  int>,std::map<typename PeriodicGlidePlane<dim>::VectorDim,std::set<std::pair<short int,short  int>>,CompareVectorsByComponent<double,dim,float>>,size_t, const T* const>> PeriodicGlidePlane<dim>::polygonPatchIntersection(const std::vector<std::pair<VectorLowerDim, const T* const>>& polyPoints, const bool& intersectInternalNodes)
{
std::vector<std::tuple<typename PeriodicGlidePlane<dim>::VectorLowerDim,typename PeriodicGlidePlane<dim>::VectorDim,std::pair<short int,short  int>, std::map<VectorDim,std::set<std::pair<short int,short  int>>,CompareVectorsByComponent<double,dim,float>> ,size_t,const T* const>> ppiPPedges; // 2dPos,shift,std::pair<edgeIDs>,LoopNodeType*
// 2dPos,shift,std::pair<edgeIDs>,std::set<std::pair<edgeIDs>>,size_t,LoopNodeType*
//This is done to take care of all the edges that the two nodes are intersecting, size_t indicates the remaining intersection
if (polyPoints.size())
{
    PeriodicPatchBoundary<dim> patchBoundary(glidePlaneFactory, referencePlane);
    auto lastPatch(patchBoundary.getPatch(findPatch(polyPoints[0].first, VectorDim::Zero())));
    // const size_t loopID((*polyPoints.begin()).second->loop()->sID);
    // const size_t runID((*polyPoints.begin()).second->network().simulationParameters.runID);

    //                std::ofstream polyFile("poly.txt");
    //                for(const auto& node : polyPoints)
    //                {
    //                    polyFile<<"    "<<node.transpose()<<std::endl;
    //                }

    //        int removeMe=0;

    //        ppi.emplace_back(polyPoints[0],lastPatch,nullptr);

    // Find other points
    for (size_t k = 0; k < polyPoints.size(); ++k)
    {
        std::vector<std::tuple<typename PeriodicGlidePlane<dim>::VectorLowerDim,typename PeriodicGlidePlane<dim>::VectorDim,std::pair<short int,short int>,const T* const>> ppi; // 2dPos,shift,std::pair<edgeIDs>,LoopNodeType*
        std::map<VectorDim,std::set<std::pair<short int,short  int>>,CompareVectorsByComponent<double,dim,float>> edgeMaps; //This will populate the set of edges that this link has crossed along with the patchShift
        std::vector<std::pair<short int,short  int>> edgeVector; //This will just give the number of edges intersected
        std::vector<std::shared_ptr<PeriodicPlanePatch<dim>>> currentPatches;
        //            std::cout<<"k point ="<<k<<std::endl;
        currentPatches.push_back(lastPatch); // keep patches alive durin current line segment search
        const auto startPoint(polyPoints[k]);
        const auto endPoint(k == polyPoints.size() - 1 ? polyPoints[0] : polyPoints[k + 1]);
        if(startPoint.second->periodicPlanePatch() != endPoint.second->periodicPlanePatch() || intersectInternalNodes)
        {
            while (true)
            {
                std::map<const PeriodicPlanePatch<dim> *, std::vector<std::pair<const VectorLowerDim, const PeriodicPlaneEdge<dim> *const>>> crossdEdges;

                for (const auto &bndEdge : patchBoundary.untwinnedEdges())
                { // loop over outer boundaries and holes
                    SegmentSegmentDistance<dim - 1> ssd(startPoint.first, endPoint.first, *bndEdge->source, *bndEdge->sink);
                    // std::cout << " Determining intersection between " << startPoint.second->tag() << " and " << endPoint.second->tag()<<" dMin "<<ssd.dMin <<" den is "<< std::flush;
                    // const double alignemntAngle(((endPoint.first - startPoint.first).normalized()).dot((*bndEdge->sink - *bndEdge->source).normalized()));
                    // const double maxAngleAlignment(cos(1.0e-4 * M_PI / 180.0)); // Alignment angle is fixed to 10 degrees

                    if (ssd.isIntersecting(FLT_EPSILON))
                    { // intersection with current boundary found
                        // std::cout << std::setprecision(15) << " Intersection found " << ssd.t << " " << ssd.u << std::endl;
                        crossdEdges[bndEdge->patch].emplace_back(0.5 * (ssd.x0 + ssd.x1), bndEdge);
                    }
                    //Changed to allow for the intersection of the near parallel segment with the boundary
                    // else
                    // {
                    //     if (startPoint.second && endPoint.second)
                    //     {
                    //         // Check the alignment of the links with the bndEdge Source
                    //         //                        std::cout<<"Coming here to check for 3D intersection"<<std::endl;
                    //         const double maxAngleAlignment(cos(5.0 * M_PI / 180.0)); //Alignment angle is fixed to 10 degrees
                    //         // std::cout<<" Alignemnt Angle =>maxAngleAlignment"<<alignemntAngle<<" "<<maxAngleAlignment<<" for nodes "<<startPoint.second->tag()
                    //         // <<"=>"<<endPoint.second->tag()<<" ssdDmin "<<ssd.dMin<<std::endl;

                    //         if (fabs(alignemntAngle) > maxAngleAlignment && ssd.dMin < 1000*FLT_EPSILON) //If very close to the boundary and aligned only then check in 3D
                    //         {
                    //             //The value of 1000 is determined by trial and error
                    //             //Determine the patches of startPoint and endPoint
                    //             if (startPoint.second->periodicPlanePatch() != endPoint.second->periodicPlanePatch())
                    //             {
                    //                 //The two loop nodes belong to different patches, a node is needed to be inserted on the edge
                    //                 //This is only valid with the if for the alignment angle and the ssd.dMin
                    //                 crossdEdges[bndEdge->patch].emplace_back(0.5 * (ssd.x0 + ssd.x1), bndEdge);
                    //             }
                    //         }
                    //     }
                    // }
                    
                }

                if (crossdEdges.size() == 0)
                { // No untwinned edges have been crossed, break and move to next pair of polyPoints
                    ppi.emplace_back(endPoint.first, lastPatch->shift, std::make_pair(-1,-1), endPoint.second);
                    // edgeMaps.emplace(VectorDim::Zero(),std::make_pair(-1,-1));
                    //                    printIntersections(ppi,currentPatches,lastPatch,removeMe,patchBoundary.untwinnedEdges());
                    break;
                }
                else
                {
                    for (const auto &pair : crossdEdges)
                    { // loop over crossed patches
                        //VectorDim shift(pair.first->shift);
                        //                        std::cout<<"crossed patches="<<crossdEdges.size()<<std::endl;
                        assert((pair.second.size()>0 && pair.second.size()<=2) && "At max two periodic plane edges can be intersected (2D geometry constraint)");

                        // if (pair.second.size()==2)
                        // {
                        //     std::cout<<" Two intersection case "<<std::endl;
                        //     std::cout<<" First "<<pair.second.begin()->first.transpose()<<std::endl;
                        //     std::cout<<" second "<<pair.second.rbegin()->first.transpose()<<std::endl;
                        //     std::cout<<" difference "<<(pair.second.begin()->first-pair.second.rbegin()->first).norm()<<std::endl;
                        //     std::cout<<" difference "<<((pair.second.begin()->first-pair.second.rbegin()->first).norm()<FLT_EPSILON)<<std::endl;
                        // }
                        
                        if (pair.second.size()==2 && (pair.second.begin()->first-pair.second.rbegin()->first).norm()<FLT_EPSILON)
                        {
                            // pair.second.size() cannot be 4 even for diagonally opposite ends since we are looping over the untwineed edges for the intersection.
                            //Need to check that two subsequen links are at the same 2D position

                            //The two intersection points are the same... Insert at diagonally opposite end

                            // std::cout<<"Coming in two intersection case "<<std::endl;
                            assert((pair.second.begin()->second!=pair.second.rbegin()->second) && "Two intersection must be different");
                            ppi.emplace_back(pair.second.begin()->first, pair.first->shift, std::make_pair(pair.second.begin()->second->edgeID, pair.second.rbegin()->second->edgeID), nullptr);
                            const auto edgeSetIter1 (edgeMaps.find(pair.first->shift));
                            edgeVector.push_back(std::make_pair(pair.second.begin()->second->edgeID, pair.second.rbegin()->second->edgeID));
                            
                            if (edgeSetIter1==edgeMaps.end())
                            {
                                std::set<std::pair<short int,short  int>> edgeSet;
                                edgeSet.insert(std::make_pair(pair.second.begin()->second->edgeID, pair.second.rbegin()->second->edgeID));
                                // edgeMaps.emplace(std::piecewise_construct,
                                //                 std::forward_as_tuple(pair.first->shift),
                                //                 std::forward_as_tuple(pair.second.begin()->second->edgeID, pair.second.rbegin()->second->edgeID));
                                edgeMaps.emplace(pair.first->shift,edgeSet);
                            }
                            else
                            {
                                edgeSetIter1->second.insert(std::make_pair(pair.second.begin()->second->edgeID, pair.second.rbegin()->second->edgeID));
                            }

                            const VectorDim localShift(pair.first->shift + pair.second.begin()->second->deltaShift+pair.second.rbegin()->second->deltaShift);
                            //here three patches are needed to be inserted... two at periodically opposite faces and one at diagonally opposite faces
                            lastPatch = patchBoundary.patches().getFromKey(pair.first->shift + pair.second.begin()->second->deltaShift);
                            currentPatches.push_back(lastPatch); // keep patches alive durin current line segment search
                            lastPatch = patchBoundary.patches().getFromKey(pair.first->shift + pair.second.rbegin()->second->deltaShift);
                            currentPatches.push_back(lastPatch); // keep patches alive durin current line segment search
                            lastPatch = patchBoundary.patches().getFromKey(localShift);
                            currentPatches.push_back(lastPatch); // keep patches alive durin current line segment search
                            assert(pair.second.begin()->second->twin); //Now the twin for both the edges must exist
                            assert(pair.second.rbegin()->second->twin);
                            assert((pair.second.begin()->second->twin!=pair.second.rbegin()->second->twin) && "Two intersection must be different after twinned link");
                            //Insert a new node at the twinned edge of the current pair of edges
                            //Need to determine the correct edgeIDs
                            short int beginTwinEdgeID(-1);
                            short rbeginTwinEdgeID(-1);
                            //We need to get the edgeID corresponding to the diagonally opposite patch
                            if (pair.second.begin()->second==pair.second.rbegin()->second->prev)
                            {
                                assert(pair.second.begin()->second->twin->prev->twin && "A twin in the diagonally opposite patch must exist");
                                assert(pair.second.rbegin()->second->twin->next->twin && "A twin in the diagonally opposite patch must exist");
                                beginTwinEdgeID=pair.second.begin()->second->twin->prev->twin->edgeID;
                                rbeginTwinEdgeID=pair.second.rbegin()->second->twin->next->twin->edgeID;
                            }
                            else if (pair.second.begin()->second==pair.second.rbegin()->second->next)
                            {
                                assert(pair.second.begin()->second->twin->next->twin && "A twin in the diagonally opposite patch must exist");
                                assert(pair.second.rbegin()->second->twin->prev->twin && "A twin in the diagonally opposite patch must exist");

                                beginTwinEdgeID=pair.second.begin()->second->twin->next->twin->edgeID;
                                rbeginTwinEdgeID=pair.second.rbegin()->second->twin->prev->twin->edgeID;
                            }
                            else
                            {
                                assert(false && "Edges must be consecutive");
                            }
                            
                            assert(beginTwinEdgeID!=-1 && "A twin ID must exist");
                            assert(rbeginTwinEdgeID!=-1 && "A twin ID must exist");

                            ppi.emplace_back(pair.second.begin()->first, localShift, std::make_pair(beginTwinEdgeID, rbeginTwinEdgeID), nullptr);
                            edgeVector.push_back(std::make_pair(beginTwinEdgeID, rbeginTwinEdgeID));

                            const auto edgeSetIter2 (edgeMaps.find(localShift));
                            if (edgeSetIter2==edgeMaps.end())
                            {
                                std::set<std::pair<short int,short  int>> edgeSet;
                                edgeSet.insert(std::make_pair(beginTwinEdgeID, rbeginTwinEdgeID));
                                // edgeMaps.emplace(std::piecewise_construct,
                                //                 std::forward_as_tuple(localShift),
                                //                 std::forward_as_tuple(pair.second.begin()->second->twin->edgeID, pair.second.rbegin()->second->twin->edgeID));
                                edgeMaps.emplace(localShift,edgeSet);
                            }
                            else
                            {
                                edgeSetIter2->second.insert(std::make_pair(beginTwinEdgeID, rbeginTwinEdgeID));
                            }
                        }
                        else
                        {
                            for (const auto &edge : pair.second)
                            { // loop over crossed edges
                                ppi.emplace_back(edge.first, pair.first->shift, std::make_pair(edge.second->edgeID,-1), nullptr);
                                edgeVector.push_back(std::make_pair(edge.second->edgeID,-1));

                                const auto edgeSetIter1(edgeMaps.find(pair.first->shift));
                                if (edgeSetIter1 == edgeMaps.end())
                                {
                                    std::set<std::pair<short int,short  int>> edgeSet;
                                    edgeSet.insert(std::make_pair(edge.second->edgeID,-1));
                                    // edgeMaps.emplace(std::piecewise_construct,
                                    //             std::forward_as_tuple(pair.first->shift),
                                    //             std::forward_as_tuple(edge.second->edgeID,-1));
                                    edgeMaps.emplace(pair.first->shift, edgeSet);
                                }
                                else
                                {
                                    edgeSetIter1->second.insert(std::make_pair(edge.second->edgeID,-1));
                                }
                                //  shift += edge.second->deltaShift;
                                const VectorDim localShift(pair.first->shift + edge.second->deltaShift);
                                // Insert patches corresponding to the individual intersections as well
                                lastPatch = patchBoundary.patches().getFromKey(localShift);
                                currentPatches.push_back(lastPatch); // keep patches alive durin current line segment search
                                assert(edge.second->twin);
                                ppi.emplace_back(edge.first, edge.second->twin->patch->shift, std::make_pair(edge.second->twin->edgeID,-1), nullptr);
                                edgeVector.push_back(std::make_pair(edge.second->twin->edgeID,-1));

                                const auto edgeSetIter2(edgeMaps.find(edge.second->twin->patch->shift));
                                if (edgeSetIter2 == edgeMaps.end())
                                {
                                    std::set<std::pair<short int,short  int>> edgeSet;
                                    edgeSet.insert(std::make_pair(edge.second->twin->edgeID,-1));
                                    // edgeMaps.emplace(std::piecewise_construct,
                                    //             std::forward_as_tuple(edge.second->twin->patch->shift),
                                    //             std::forward_as_tuple(edge.second->twin->edgeID,-1));
                                    // edgeMaps.emplace(edge.second->twin->patch->shift, std::make_pair(edge.second->twin->edgeID,-1));
                                    edgeMaps.emplace(edge.second->twin->patch->shift, edgeSet);
                                }
                                else
                                {
                                    edgeSetIter2->second.insert(std::make_pair(edge.second->twin->edgeID,-1));
                                }
                            }
                        }
                    }
                    if (lastPatch->contains(endPoint.first))
                    {
                        ppi.emplace_back(endPoint.first, lastPatch->shift, std::make_pair(-1,-1), endPoint.second);
                        // edgeMaps.emplace(VectorDim::Zero(), std::make_pair(-1, -1));
                        //                        printIntersections(ppi,currentPatches,lastPatch,removeMe,patchBoundary.untwinnedEdges());
                        break;
                    }
                }
            }
        }
        else
        {
            ppi.emplace_back(endPoint.first, lastPatch->shift, std::make_pair(-1,-1), endPoint.second);
        }

        //Here populate the ppiPPedges
        size_t crossCount=0;
        for (const auto& ppiIter : ppi )
        {
            // //LoopNode* is the last
            // if (edgeVector.size()==0)
            // {
            //     ppiPPedges.emplace_back(std::get<0>(ppiIter),std::get<1>(ppiIter),std::get<2>(ppiIter),edgeMaps,0,std::get<3>(ppiIter));
            // }
            // else
            // {
            size_t edgesStillRemaining(0);
            if (std::get<3>(ppiIter) == nullptr)
            {
                //That mean intersection with the edge is here
                crossCount++;
                edgesStillRemaining = edgeVector.size() - crossCount;
            }
            ppiPPedges.emplace_back(std::get<0>(ppiIter), std::get<1>(ppiIter), std::get<2>(ppiIter), edgeMaps, edgesStillRemaining, std::get<3>(ppiIter));
            // }
        }
        if (edgeVector.size() != 0)
        {
            if ((edgeVector.size() - crossCount ) != 0)
            {
                std::cout<<"Edge VectorSize "<<edgeVector.size()<<std::endl;
                std::cout<<"Cross count "<<crossCount<<std::endl;
            }
            assert((edgeVector.size() - crossCount) == 0 && "No edges should be remaining after completion of polygonPatchIntersection");
        }
        // if (true)
        // {
        //     std::ofstream edgesFile("Debug/edges_"+ std::to_string(runID) + "_" + std::to_string(loopID) + ".txt", std::ofstream::out | std::ofstream::app);
        //     std::ofstream pointsFile("Debug/points_"+ std::to_string(runID) + "_" + std::to_string(loopID) + ".txt",std::ofstream::out | std::ofstream::app);

        //     for (const auto &edge : lastPatch->edges())
        //     {
        //         edgesFile << edge->source->sID << " " << edge->sink->sID << " " << lastPatch->sID << std::endl;
        //         pointsFile << "    " << edge->source->sID << " " << (*edge->source.get()).transpose() << std::endl;
        //     }
        //     edgesFile.close();
        //     pointsFile.close();
        // }
    }

    // if (true)
    // {

    //     std::ofstream polyFile("Debug/poly_"+std::to_string(runID) + "_" + std::to_string(loopID) + ".txt");
    //     //                    std::ofstream poly3DFile("poly3D.txt");
    //     for (const auto &tup : ppi)
    //     {
    //         //                        if(!node.second.expired())
    //         //                        {
    //         polyFile << "    " << std::get<0>(tup).transpose() << std::endl;
    //         //                        }
    //     }
    //     polyFile.close();
    // }
}


return ppiPPedges;
}

template<typename T1,typename T2, typename T3,typename T4>
void printIntersections(const T1& ppi,const T2& usedPatches,const T3& lastPatch,int& removeMe, const T4& untwinnedEdges)
{
std::ofstream pointsFile("points_"+std::to_string(removeMe)+".txt");
for(const auto& tup : ppi)
{
    //                        if(!node.second.expired())
    //                        {
    pointsFile<<"    "<<std::get<0>(tup).transpose()<<std::endl;
    //                        }
}

std::ofstream edgesFile("edges_"+std::to_string(removeMe)+".txt");
for(const auto& patch : usedPatches)
{
    for(const auto& edge : patch->edges())
    {
        if(untwinnedEdges.find(edge.get())!=untwinnedEdges.end())
        {
            edgesFile<<patch->sID<<" "<<edge->source->transpose()<<" "<<edge->sink->transpose()<<" 1"<<std::endl;
        }
        else
        {
            edgesFile<<patch->sID<<" "<<edge->source->transpose()<<" "<<edge->sink->transpose()<<" 0"<<std::endl;
        }
            
    }
}

std::ofstream lastPatchFile("lastPatch_"+std::to_string(removeMe)+".txt");
for(const auto& edge : lastPatch->edges())
{
    lastPatchFile<<lastPatch->sID<<" "<<edge->source->transpose()<<" "<<edge->sink->transpose()<<std::endl;
}
removeMe++;
}
    

}
#endif
