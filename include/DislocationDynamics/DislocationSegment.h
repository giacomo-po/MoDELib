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
    // /*                      */,public ConfinedDislocationObject<dim>
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
        void updateQuadraturePointsSeg();
        /**********************************************************************/
        GlidePlaneContainerType glidePlanes() const
        {
            GlidePlaneContainerType temp;
            for (const auto &ln : this->loopLinks())
            {
                if (ln->periodicPlanePatch())
                {
                    temp.emplace(ln->periodicPlanePatch()->glidePlane.get());
                }
            }
            return temp;
        }

        /**********************************************************************/
        PlanarMeshFaceContainerType meshFaces() const
        {
            PlanarMeshFaceContainerType temp;
            const PlanarMeshFaceContainerType sourceMeshFaces(this->source->meshFaces());
            const PlanarMeshFaceContainerType sinkMeshFaces(this->sink->meshFaces());

            for (const auto &sourceMeshFace : sourceMeshFaces)
            {
                if (sinkMeshFaces.find(sourceMeshFace) != sinkMeshFaces.end())
                {
                    temp.emplace(sourceMeshFace);
                }
            }
            return temp;
        }

        /**********************************************************************/
        bool isOnExternalBoundary() const
        { /*!\returns _isOnExternalBoundarySegment.
           */
            bool _isOnExternalBoundary(false);
            for (const auto &face : meshFaces())
            {
                if (face->regionIDs.first == face->regionIDs.second)
                {
                    _isOnExternalBoundary = true;
                }
            }

            return _isOnExternalBoundary;
        }

        /**********************************************************************/
        bool isOnInternalBoundary() const
        {
            bool _isOnInternalBoundary(false);
            for (const auto &face : meshFaces())
            {
                if (face->regionIDs.first != face->regionIDs.second)
                {
                    _isOnInternalBoundary = true;
                }
            }
            return _isOnInternalBoundary;
        }

        // /**********************************************************************/
        bool isOnBoundary() const
        {
            return isOnExternalBoundary() || isOnInternalBoundary();
        }

        /**********************************************************************/
        VectorDim bndNormal() const
        {
            VectorDim _outNormal(VectorDim::Zero());
            for (const auto &face : meshFaces())
            {
                _outNormal += face->outNormal();
            }
            const double _outNormalNorm(_outNormal.norm());
            if (_outNormalNorm > FLT_EPSILON)
            {
                _outNormal /= _outNormalNorm;
            }
            else
            {
                _outNormal.setZero();
            }
            return _outNormal;
        }
        /**********************************************************************/
        GrainContainerType grains() const
        {
            GrainContainerType temp;
            for (const auto &glidePlane : glidePlanes())
            {
                temp.insert(&glidePlane->grain);
            }
            return temp;
        }
        /**********************************************************************/
        std::vector<std::pair<const GlidePlane<dim> *const, const GlidePlane<dim> *const>> parallelAndCoincidentGlidePlanes(const GlidePlaneContainerType &other) const
        {
            std::vector<std::pair<const GlidePlane<dim> *const, const GlidePlane<dim> *const>> pp;

            for (const auto &plane : glidePlanes())
            {
                for (const auto &otherPlane : other)
                {
                    //                    if(plane!=otherPlane && plane->n.cross(otherPlane->n).squaredNorm()==0)
                    if (plane->n.cross(otherPlane->n).squaredNorm() == 0)
                    { // parallel planes
                        pp.emplace_back(plane, otherPlane);
                    }
                }
            }
            return pp;
        }
        /**********************************************************************/
        VectorDim snapToGlidePlanes(const VectorDim &P)
        {
            GlidePlaneContainerType gps(glidePlanes());
            switch (gps.size())
            {
            case 0:
            {
                assert(false && "Glide plane size must be larger than 0");
                return P;
                break;
            }
            case 1:
            {
                return (*gps.begin())->snapToPlane(P);
                break;
            }
            case 2:
            {
                const PlanePlaneIntersection<dim> ppi(**gps.begin(), **gps.rbegin());
                assert(ppi.type == PlanePlaneIntersection<dim>::INCIDENT && "Intersection must be incident");
                return ppi.P + (P - ppi.P).dot(ppi.d) * ppi.d;
                break;
            }
            case 3:
            {
                const PlanePlaneIntersection<dim> ppi(**gps.begin(), **gps.rbegin());
                assert(ppi.type == PlanePlaneIntersection<dim>::INCIDENT && "Intersection must be incident");
                const auto iterP(++gps.begin());
                const PlaneLineIntersection<dim> pli((*iterP)->P, (*iterP)->unitNormal, ppi.P, ppi.d);
                assert(pli.type == PlaneLineIntersection<dim>::INCIDENT && "Plane line intersection must be incident");
                return pli.P;
                break;
            }
            default:
            {
                const auto iterPlane1(gps.begin());
                const auto iterPlane2(++gps.begin());
                auto iterPlane3(++(++gps.begin()));

                const PlanePlaneIntersection<dim> ppi(**iterPlane1, **iterPlane2);
                assert(ppi.type == PlanePlaneIntersection<dim>::INCIDENT && "Intersection must be incident");
                const PlaneLineIntersection<dim> pli((*iterPlane3)->P, (*iterPlane3)->unitNormal, ppi.P, ppi.d);
                const VectorDim snappedPos(pli.P);
                assert(pli.type == PlaneLineIntersection<dim>::INCIDENT && "Plane line intersection must be incident");
                // if (pli.type!=PlaneLineIntersection<dim>::INCIDENT)
                // {
                //     std::cout<<" IterPlane1 "<<(*iterPlane1)->P.transpose()<<"\t"<<(*iterPlane1)->unitNormal.transpose()<<std::endl;
                //     std::cout<<" IterPlane2 "<<(*iterPlane2)->P.transpose()<<"\t"<<(*iterPlane2)->unitNormal.transpose()<<std::endl;
                //     std::cout<<" IterPlane3 "<<(*iterPlane3)->P.transpose()<<"\t"<<(*iterPlane3)->unitNormal.transpose()<<std::endl;
                //     std::cout<<"ppi infor "<<std::endl;
                //     std::cout<<ppi.P.transpose()<<"\t"<<ppi.d.transpose()<<std::endl;
                //     std::cout<<" glide plane size "<<gps.size()<<" Printing glide planes "<<std::endl;
                //     for (const auto& gp : gps)
                //     {
                //         std::cout<<gp->P.transpose()<<"\t"<<gp->unitNormal.transpose()<<"\t"<<gp->planeIndex<<std::endl;
                //     }
                //     std::cout<<" intersection type "<<pli.type<<std::endl;
                // }
                while (++iterPlane3 != gps.end())
                {
                    if (!(*iterPlane3)->contains(snappedPos))
                    {
                        std::cout << " Glide Plane " << (*iterPlane3)->P.transpose() << " " << (*iterPlane3)->unitNormal.transpose() << std::endl;
                        assert(false && "Plane must contain the position for glide plane size >=3");
                    }
                }
                return snappedPos;
                break;
            }
            }
        }
    };
//    {
//        typedef typename TypeTraits<DislocationSegment<dim,corder>>::NodeType NodeType;
//        typedef PlanarDislocationSegment<DislocationSegment<dim,corder>> PlanarDislocationSegmentType;
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
