/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po         <gpo@ucla.edu>.
 * Copyright (C) 2011 by Benjamin Ramirez   <ramirezbrf@gmail.com>.
 * Copyright (C) 2011 by Tamer Crsoby       <tamercrosby@gmail.com>,
 * Copyright (C) 2011 by Can Erel           <canerel55@gmail.com>,
 * Copyright (C) 2011 by Mamdouh Mohamed    <msm07d@fsu.edu>
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationNode_h_
#define model_DislocationNode_h_

//#include <PlanarDislocationNode.h>

#include <DislocationDynamicsModule.h>
#include <PlaneLineIntersection.h>
#ifndef NDEBUG
#define VerboseDislocationNode(N,x) if(this->network().verboseDislocationNode>=N){std::cout<<x;}
#else
#define VerboseDislocationNode(N,x)
#endif

namespace model
{
    
    template <int dim, short unsigned int corder>
    class DislocationNode : public NetworkNode<DislocationNode<dim,corder>>
    /*                   */,public SplineNode<DislocationNode<dim,corder>,dim,corder,Hermite>
    // /*                   */,public ConfinedDislocationObject<dim>
    {
        
        public:
        
        typedef TypeTraits<DislocationLoop<dim,corder>> TraitsType;
        typedef typename TraitsType::LoopNetworkType LoopNetworkType;
        typedef typename TraitsType::LoopType LoopType;
        typedef typename TraitsType::LoopNodeType LoopNodeType;
        typedef typename TraitsType::LoopLinkType LoopLinkType;
        typedef typename TraitsType::NetworkNodeType NetworkNodeType;
        typedef typename TraitsType::NetworkLinkType NetworkLinkType;
        typedef typename TraitsType::FlowType FlowType;
        typedef typename TraitsType::VectorDim VectorDim;
        typedef SplineNode<DislocationNode<dim,corder>,dim,corder,Hermite> SplineNodeType;
        typedef ConfinedDislocationObject<dim> ConfinedDislocationObjectType;
        typedef typename TypeTraits<NetworkNodeType>::MeshLocation MeshLocation;
        typedef std::vector<VectorDim> VectorOfNormalsType;
        typedef GlidePlane<dim> GlidePlaneType;
        typedef std::set<const GlidePlaneType *> GlidePlaneContainerType;
        typedef PlanarMeshFace<dim> PlanarMeshFaceType;
        typedef std::set<const PlanarMeshFaceType *> PlanarMeshFaceContainerType;
        typedef Grain<dim> GrainType;
        typedef std::set<const GrainType*> GrainContainerType;
        // static bool use_velocityFilter;
        // static double velocityReductionFactor;
        // static int verboseDislocationNode;

        const Simplex<dim,dim>* p_Simplex;
        VectorDim velocity;
        VectorDim vOld;         //! The previous velocity vector of *this PlanarDislocationNode
        double velocityReductionCoeff;
        std::shared_ptr<NetworkNodeType> virtualNode;
        NetworkNodeType* const masterNode;

        
        DislocationNode(LoopNetworkType* const,const VectorDim&,const VectorDim&,const double&);
        ~DislocationNode();
        std::shared_ptr<DislocationNode> clone() const;
        const Simplex<dim,dim>* get_includingSimplex(const VectorDim&,const Simplex<dim,dim>* const) const;
        const Simplex<dim,dim>* includingSimplex() const;
        bool isMovableTo(const VectorDim&) const;
        const double& velocityReduction() const;
        bool isVirtualBoundaryNode() const __attribute__ ((deprecated));
        bool isBoundaryNode() const;
        bool isGrainBoundaryNode() const;
        void updateGeometry();
        bool set_P(const VectorDim&);
        bool trySet_P(const VectorDim&);
        const VectorDim& get_V() const;
        // const VectorDim& get_P() const;
        MeshLocation meshLocation() const;
        void set_V(const VectorDim& vNew);
        void projectVelocity();
        // size_t uniqueLoopNodes() const;
        // void addLoopNode(LoopNodeType* const);
        // void removeLoopNode(LoopNodeType* const);
        // bool isRemovable(const double& Lmin,const double& relAreaThIn);
        VectorDim invariantDirectionOfMotion() const;
        const std::shared_ptr<NetworkNodeType>& virtualBoundaryNode() const;
        // static void initFromFile(const std::string&);
        // bool isZeroBurgersNode() const;
        /**********************************************************************/
        GlidePlaneContainerType glidePlanes() const
        {
            GlidePlaneContainerType temp;
            for (const auto &ln : this->loopNodes())
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
            for (const auto&ln : this->loopNodes() )
            {
                if (ln->periodicPlaneEdge.first)
                {
                  temp=ln->periodicPlaneEdge.first->meshIntersection->faces;  
                }
                if (ln->periodicPlaneEdge.second)
                {
                    for (const auto& face : ln->periodicPlaneEdge.second->meshIntersection->faces)
                    {
                          temp.emplace(face);  
                    }
                }
            }
            return temp;
        }
        /**********************************************************************/
        bool isOnExternalBoundary() const
        { /*!\returns _isOnExternalBoundarySegment.
           */
            bool _isOnExternalBoundary(false);
            for(const auto& face : meshFaces())
            {                
                if(face->regionIDs.first==face->regionIDs.second)
                {
                    _isOnExternalBoundary=true;
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
            return  isOnExternalBoundary() || isOnInternalBoundary();
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
        VectorDim snapToGlidePlanesinPeriodic(const VectorDim& P)
        {
            GlidePlaneContainerType gps(glidePlanes());
            switch(gps.size())
            {
                case 0 :
                {
                    assert (false && "Glide plane size must be larger than 0");
                    return P;
                    break;
                }
                case 1 :
                {
                    return (*gps.begin())->snapToPlane(P);
                    break;
                }
                case 2 : 
                {
                    const PlanePlaneIntersection<dim> ppi(**gps.begin(),**gps.rbegin());
                    assert(ppi.type==PlanePlaneIntersection<dim>::INCIDENT && "Intersection must be incident");
                    return ppi.P+(P-ppi.P).dot(ppi.d)*ppi.d;
                    break;
                }
                case 3 : 
                {
                    const PlanePlaneIntersection<dim> ppi(**gps.begin(),**gps.rbegin());
                    assert(ppi.type==PlanePlaneIntersection<dim>::INCIDENT && "Intersection must be incident");
                    const auto iterP(++gps.begin());
                    const PlaneLineIntersection<dim> pli((*iterP)->P,(*iterP)->unitNormal,ppi.P,ppi.d);
                    assert(pli.type==PlaneLineIntersection<dim>::INCIDENT && "Plane line intersection must be incident");
                    return pli.P;
                    break;

                }
                default :
                {
                    const auto iterPlane1(gps.begin());
                    const auto iterPlane2(++gps.begin());
                    auto iterPlane3(++(++gps.begin()));

                    const PlanePlaneIntersection<dim> ppi(**iterPlane1,**iterPlane2);
                    assert(ppi.type==PlanePlaneIntersection<dim>::INCIDENT && "Intersection must be incident");
                    const PlaneLineIntersection<dim> pli((*iterPlane3)->P,(*iterPlane3)->unitNormal,ppi.P,ppi.d);
                    const VectorDim snappedPos(pli.P);
                    assert(pli.type==PlaneLineIntersection<dim>::INCIDENT && "Plane line intersection must be incident");
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
                    while (++iterPlane3 !=gps.end())
                    {
                        if (!(*iterPlane3)->contains(snappedPos))
                        {
                            std::cout<<" Glide Plane "<<(*iterPlane3)->P.transpose()<<" "<<(*iterPlane3)->unitNormal.transpose()<<std::endl;
                            assert(false && "Plane must contain the position for glide plane size >=3");
                        }
                    }
                    return snappedPos;
                    break;
                }
            }
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

    };
    
}
#endif
