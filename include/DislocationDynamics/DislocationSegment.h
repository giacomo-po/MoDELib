/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po       <gpo@ucla.edu>.
 * Copyright (C) 2011 by Benjamin Ramirez <ramirezbrf@gmail.com>.
 * Copyright (C) 2011 by Tamer Crsoby     <tamercrosby@gmail.com>,
 * Copyright (C) 2011 by Can Erel         <canerel55@gmail.com>,
 * Copyright (C) 2011 by Mamdouh Mohamed  <msm07d@fsu.edu>
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationSegment_h_
#define model_DislocationSegment_h_

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <DislocationDynamicsModule.h>
#include <TypeTraits.h>
#include <NetworkLink.h>
#include <SplineSegment.h>
#include <ConfinedDislocationObject.h>
#include <DislocationQuadraturePoint.h>
#include <DislocationLoopIO.h>
#include <PlaneLineIntersection.h>

#ifndef NDEBUG
#define VerboseDislocationSegment(N,x) if(verboseDislocationSegment>=N){std::cout<<x;}
#else
#define VerboseDislocationSegment(N,x)
#endif

namespace model
{
    template <int dim, short unsigned int corder>
    class DislocationSegment : public NetworkLink<DislocationSegment<dim,corder>>
    /*                      */,public SplineSegment<dim,corder>
    /*                      */,public DislocationQuadraturePointContainer<dim,corder>
    {
        
    public:
        
        typedef TypeTraits<DislocationSegment<dim,corder>> TraitsType;
        typedef typename TraitsType::LoopNetworkType LoopNetworkType;
        typedef typename TraitsType::LoopType LoopType;
        typedef typename TraitsType::LoopNodeType LoopNodeType;
        typedef typename TraitsType::LoopLinkType LoopLinkType;
        typedef typename TraitsType::NetworkNodeType NetworkNodeType;
        typedef typename TraitsType::NetworkLinkType NetworkLinkType;
        typedef typename TraitsType::FlowType FlowType;
        typedef typename TraitsType::VectorDim VectorDim;
        typedef typename TraitsType::MatrixDim MatrixDim;
        typedef typename TraitsType::MeshLocation MeshLocation;
        typedef ConfinedDislocationObject<dim> ConfinedDislocationObjectType;
        typedef SplineSegment<dim,corder> SplineSegmentType;
        typedef typename SplineSegmentType::VectorNdof VectorNdof;
        typedef typename SplineSegmentType::MatrixNdof MatrixNdof;
        typedef typename SplineSegmentType::VectorNcoeff VectorNcoeff;
        static constexpr int Ndof  = SplineSegmentType::Ndof;        // number of Hermite dofs
        static constexpr int Ncoeff  = SplineSegmentType::Ncoeff;        // number of Hermite dofs
        typedef Grain<dim> GrainType;
        typedef std::set<const GrainType *> GrainContainerType;
        typedef GlidePlane<dim> GlidePlaneType;
        typedef std::set<const GlidePlaneType *> GlidePlaneContainerType;
        typedef PlanarMeshFace<dim> PlanarMeshFaceType;
        typedef std::set<const PlanarMeshFaceType *> PlanarMeshFaceContainerType;
        typedef typename DislocationQuadraturePoint<dim,corder>::QuadratureDynamicType QuadratureDynamicType;
        
//        std::map<size_t,
//        /*    */ std::pair<VectorNcoeff,VectorDim>,
//        /*    */ std::less<size_t>
//        /*    */ > h2posMap;
//        Eigen::Matrix<double, Ndof, Eigen::Dynamic> Mseg;
        MatrixNdof Kqq; //! Segment Stiffness Matrix
        VectorNdof Fq; //! Segment Nodal Force Vector
        VectorDim Burgers; //! The Burgers vector
        double BurgersNorm;
        StraightDislocationSegment<dim> straight;
        std::shared_ptr<SlipSystem> _slipSystem;

        
        static const MatrixDim I;
        static const VectorDim zeroVector;
        static double quadPerLength;
//        static bool assembleWithTangentProjection;
        static int verboseDislocationSegment;

        
        DislocationSegment(LoopNetworkType* const net,
                         const std::shared_ptr<NetworkNodeType>&,
                         const std::shared_ptr<NetworkNodeType>&
                         );
        
        void addLoopLink(LoopLinkType* const pL);
        void removeLoopLink(LoopLinkType* const pL);
        void updateSlipSystem();
        static void initFromFile(const std::string&);

        MeshLocation meshLocation() const;
        bool isBoundarySegment() const;
        bool isGrainBoundarySegment() const;
        void updateGeometry();
        const std::shared_ptr<SlipSystem>& slipSystem() const;
        bool hasZeroBurgers() const;
//        bool isVirtualBoundarySegment() const;
        bool isGlissile() const;        
        bool isSessile() const;
        void assembleGlide(const bool&);
        const VectorDim& burgers() const;
        const VectorDim& glidePlaneNormal() const;
        void addToGlobalAssembly(std::deque<Eigen::Triplet<double> >& kqqT,Eigen::VectorXd& FQ) const;
        void updateQuadraturePointsSeg();
        int velocityGroup(const double &maxVelocity, const std::set<int> &subcyclingBins) const;
        GlidePlaneContainerType glidePlanes() const;
        PlanarMeshFaceContainerType meshFaces() const;
        bool isOnExternalBoundary() const;
        bool isOnInternalBoundary() const;
        bool isOnBoundary() const;
        VectorDim bndNormal() const;
        GrainContainerType grains() const;
        std::vector<std::pair<const GlidePlane<dim> *const, const GlidePlane<dim> *const>> parallelAndCoincidentGlidePlanes(const GlidePlaneContainerType &other) const;
        VectorDim snapToGlidePlanes(const VectorDim &P) const;
        
    };
}
#endif
