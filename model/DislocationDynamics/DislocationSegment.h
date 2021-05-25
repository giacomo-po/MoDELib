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
#include <TypeTraits.h>
#include <NetworkLink.h>
#include <SplineSegment.h>
#include <ConfinedDislocationObject.h>
#include <DislocationQuadraturePoint.h>

#ifndef NDEBUG
#define VerboseDislocationSegment(N,x) if(verboseDislocationSegment>=N){model::cout<<x;}
#else
#define VerboseDislocationSegment(N,x)
#endif

namespace model
{
    template <int dim, short unsigned int corder, typename InterpolationType>
    class DislocationSegment : public NetworkLink<DislocationSegment<dim,corder,InterpolationType>>
    /*                      */,public SplineSegment<dim,corder>
    /*                      */,public ConfinedDislocationObject<dim>
    /*                      */,public DislocationQuadraturePointContainer<dim,corder>
    {
        
    public:
        
        typedef TypeTraits<DislocationSegment<dim,corder,InterpolationType>> TraitsType;
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

        std::map<size_t,
        /*    */ std::pair<VectorNcoeff,VectorDim>,
        /*    */ std::less<size_t>
        /*    */ > h2posMap;
        Eigen::Matrix<double, Ndof, Eigen::Dynamic> Mseg;
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
        bool isVirtualBoundarySegment() const;
        bool isGlissile() const;        
        bool isSessile() const;
        void assembleGlide();
        const VectorDim& burgers() const;
        const VectorDim& glidePlaneNormal() const;
        void addToGlobalAssembly(std::deque<Eigen::Triplet<double> >& kqqT,Eigen::VectorXd& FQ) const;

    };
//    {
//        typedef typename TypeTraits<DislocationSegment<dim,corder,InterpolationType>>::NodeType NodeType;
//        typedef PlanarDislocationSegment<DislocationSegment<dim,corder,InterpolationType>> PlanarDislocationSegmentType;
//        static constexpr int Ncoeff=PlanarDislocationSegmentType::Ncoeff;
//        static constexpr int pOrder=PlanarDislocationSegmentType::pOrder;
//        typedef typename PlanarDislocationSegmentType::VectorDim VectorDim;
//        typedef typename PlanarDislocationSegmentType::MatrixDim MatrixDim;
//        typedef DislocationParticle<dim> DislocationParticleType;
//
//#ifdef _MODEL_GREATWHITE_
//#include <DislocationSegmentGreatWhite.h>
//#endif
//
//        /**********************************************************************/
//        DislocationSegment(const std::shared_ptr<NodeType>& nI,
//                           const std::shared_ptr<NodeType>& nJ) :
//        /* init */ PlanarDislocationSegmentType(nI,nJ)
//        {
//
//        }
//
//        /**********************************************************************/
//        const MatrixDim& midPointStress() const __attribute__ ((deprecated))
//        {/*!\returns The stress matrix for the centre point over this segment.*/
//            return this->quadraturePoints().size()? quadraturePoint(this->quadraturePoints().size()/2).stress : MatrixDim::Zero();
//        }
//    };
}
#endif
