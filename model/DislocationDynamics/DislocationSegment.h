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


#ifndef model_DISLOCATIONSEGMENT_H
#define model_DISLOCATIONSEGMENT_H

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include <set>


#include <model/Network/Operations/EdgeExpansion.h>


#include <model/Quadrature/Quadrature.h>
#include <model/Quadrature/QuadPow.h>
#include <model/DislocationDynamics/DislocationNetworkTraits.h>
#include <model/DislocationDynamics/DislocationConsts.h>
#include <model/Geometry/Splines/SplineSegmentBase.h>
#include <model/DislocationDynamics/Materials/Material.h>
#include <model/DislocationDynamics/Materials/CrystalOrientation.h>

#include <model/DislocationDynamics/DislocationSharedObjects.h>
#include <model/DislocationDynamics/GlidePlanes/GlidePlaneObserver.h>
#include <model/DislocationDynamics/GlidePlanes/GlidePlane.h>
#include <model/Geometry/Splines/Intersection/PlanarSplineImplicitization.h>
#include <model/Geometry/Splines/Coeff2Hermite.h>

//#include <model/DislocationDynamics/NearestNeighbor/DislocationQuadratureParticle.h>
#include <model/DislocationDynamics/NearestNeighbor/DislocationParticle.h>
#include <model/ParticleInteraction/ParticleSystem.h>

#include <model/DislocationDynamics/DislocationLocalReference.h>
#include <model/DislocationDynamics/Junctions/DislocationSegmentIntersection.h>

//#include <model/DislocationDynamics/CrossSlip/CrossSlipSegment.h>
//#include <model/DislocationDynamics/DislocationMobility.h>
#include <model/Utilities/UniqueOutputFile.h>

namespace model
{
    
    
    template <short unsigned int dim>
    struct PlanarDislocationSegment
    {/*! Class template used ad a base of DislocationSegment used to
      * initialize glidePlaneNormal and sessilePlaneNormal before the base
      * NetworkLink.
      */
        
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef LatticeVector<dim> LatticeVectorType;
        
        //! The glide plane unit normal vector
        const LatticePlane   glidePlane;
        const LatticePlane sessilePlane;
        
        VectorDim   glidePlaneNormal;
        VectorDim sessilePlaneNormal;
        
        
        /**********************************************************************/
        template <typename LinkType>
        static const LatticePlaneBase& find_glidePlane(const LatticeVectorType& sourceL,
                                         const LatticeVectorType& sinkL,
                                         const LatticeVectorType& Burgers,
                                         const EdgeRef<LinkType>& ee)
        {
        
            const LatticeVectorType& linkSourceL=ee.E.source->get_L();
            const LatticeVectorType& linkSinkL=ee.E.sink->get_L();
            
            const ReciprocalLatticeVector<dim> triagleNormal=(sourceL-linkSourceL).cross(sinkL-linkSourceL)+(sourceL-linkSinkL).cross(sinkL-linkSinkL);
            
            
            return (triagleNormal.cross(ee.E.glidePlane.n).squaredNorm()) ?  CrystalOrientation<dim>::find_glidePlane(sinkL-sourceL,Burgers) :
            /*                                                           */  ee.E.glidePlane.n;
        }
        
        /**********************************************************************/
        template <typename LinkType>
        static const LatticePlaneBase& find_sessilePlane(const LatticeVectorType& sourceL,
                                           const LatticeVectorType& sinkL,
                                           const LatticeVectorType& Burgers,
                                           const EdgeRef<LinkType>& ee)
        {
            
            const LatticeVectorType& linkSourceL=ee.E.source->get_L();
            const LatticeVectorType& linkSinkL=ee.E.sink->get_L();
            
            const ReciprocalLatticeVector<dim> triagleNormal=(sourceL-linkSourceL).cross(sinkL-linkSourceL)+(sourceL-linkSinkL).cross(sinkL-linkSinkL);
            
            
            return (triagleNormal.cross(ee.E.glidePlane.n).squaredNorm()) ?  CrystalOrientation<dim>::find_sessilePlane(sinkL-sourceL,Burgers) :
            /*                                                           */  ee.E.sessilePlane.n;
        }
        
        /**********************************************************************/
        PlanarDislocationSegment(const LatticeVectorType& sourceL,
                                 const LatticeVectorType& sinkL,
                                 const LatticeVectorType& Burgers) :
        /* init list       */   glidePlane(sourceL,CrystalOrientation<dim>::find_glidePlane(sinkL-sourceL,Burgers)),
        /* init list       */ sessilePlane(sourceL,CrystalOrientation<dim>::find_sessilePlane(sinkL-sourceL,Burgers)),
        /* init list       */ glidePlaneNormal(glidePlane.n.cartesian().normalized()),
        /* init list       */ sessilePlaneNormal(sessilePlane.n.cartesian().normalized())
        {
            
        }
        
        /**********************************************************************/
        template <typename LinkType>
        PlanarDislocationSegment(const LatticeVectorType& sourceL,
                                 const LatticeVectorType& sinkL,
                                 const LatticeVectorType& Burgers,
                                 const EdgeRef<LinkType>& ee) :
        /* init list       */   glidePlane(sourceL,  find_glidePlane(sourceL,sinkL,Burgers,ee)),
        /* init list       */ sessilePlane(sourceL,find_sessilePlane(sourceL,sinkL,Burgers,ee)),
        /* init list       */ glidePlaneNormal(glidePlane.n.cartesian().normalized()),
        /* init list       */ sessilePlaneNormal(sessilePlane.n.cartesian().normalized())
        {
        
        }

        
    };
    
    /******************************************************************/
    /******************************************************************/
    template <short unsigned int _dim, short unsigned int corder, typename InterpolationType,
    /*	   */ template <short unsigned int, size_t> class QuadratureRule>
    class DislocationSegment : public PlanarDislocationSegment<_dim>,
    /*	                    */ public SplineSegmentBase<DislocationSegment<_dim,corder,InterpolationType,QuadratureRule>,
    /*                                              */ _dim, corder>,
    /*	                    */ public GlidePlaneObserver<DislocationSegment<_dim,corder,InterpolationType,QuadratureRule> >
    {
        
        
    public:
        
        enum{dim=_dim}; // make dim available outside class
        //        enum{qOrder=16};
        
        
        typedef DislocationSegment<dim,corder,InterpolationType,QuadratureRule> Derived; 		// Define "Derived" so that NetworkTypedefs.h can be used
#include <model/Network/NetworkTypedefs.h>
#include <model/Geometry/Splines/SplineEnums.h>
        
        typedef PlanarDislocationSegment<dim> PlanarSegmentType;
        typedef SplineSegmentBase<Derived,dim,corder> SegmentBaseType;
        typedef std::map<size_t,LinkType* const> AddressMapType;
        typedef typename AddressMapType::iterator AddressMapIteratorType;
        //        typedef Eigen::Matrix<double,dim,qOrder>	MatrixDimQorder;
        typedef Eigen::Matrix<double,dim,Eigen::Dynamic>	MatrixDimQorder;
        
        //        typedef Eigen::Matrix<double, 1, qOrder>	VectorQorder;
        typedef Eigen::Matrix<double, 1, Eigen::Dynamic>	VectorQorder;
        
        //        typedef QuadPow<Ncoeff-1,qOrder,QuadratureRule> QuadPowType;
        typedef typename TypeTraits<Derived>::QuadPowDynamicType QuadPowDynamicType;
        typedef typename TypeTraits<Derived>::QuadratureDynamicType QuadratureDynamicType;
        
        typedef typename GlidePlaneObserver<LinkType>::GlidePlaneType GlidePlaneType;
        typedef typename GlidePlaneObserver<LinkType>::GlidePlaneSharedPtrType GlidePlaneSharedPtrType;
        
        
        
        typedef DislocationParticle<dim> DislocationParticleType;
        typedef std::vector<DislocationParticleType*> QuadratureParticleContainerType;
        typedef LatticeVector<dim> LatticeVectorType;
        typedef ReciprocalLatticeDirection<dim> ReciprocalLatticeDirectionType;
        
        /******************************************************************/
    private: //  data members
        /******************************************************************/
        //! Parametric tangents at the quadrature points
        MatrixDimQorder rugauss;
        //! Scalar jacobian corrersponding to the quadrature points
        VectorQorder jgauss;
        //! Tangents corrersponding to the quadrature points
        MatrixDimQorder rlgauss;
        //! Matrix of shape functions at the quadrature points
        //        Eigen::Matrix<double,qOrder,Ncoeff> SFgauss;
        Eigen::Matrix<double,Eigen::Dynamic,Ncoeff> SFgauss;
        
        
        //enum {Nslips=MaterialType::Nslips};
        
        std::set<size_t> segmentDOFs;
        Eigen::Matrix<double, Ndof, Eigen::Dynamic> Mseg;
        
        //! Matrix of PK-force at quadrature points
        //		MatrixDimQorder pkGauss;
        //! Segment Stiffness Matrix
        MatrixNdof Kqq;
        //! Segment Nodal Force Vector
        VectorNdof Fq;
        //! The identity matrix
        static const Eigen::Matrix<double,_dim,_dim> I;
        
        //        //! A static vector of zeros
        //        static const Eigen::Matrix<double,_dim,1> zeroDim;
        
        //        VectorDim boundaryLoopNormal;
        
        
        /******************************************************************/
    public: //  data members
        /******************************************************************/
        
        
        QuadratureParticleContainerType quadratureParticleContainer;
        
        
        //! The Burgers vector
        const VectorDim Burgers;
        const bool isSessile;
        const std::deque<const LatticePlaneBase*> conjugatePlaneNormals;
        
        //        const LatticeVectorType& Burgers;
        
        
        
        DislocationSharedObjects<dim> shared;
        
        
        //! A shared pointer to the GlidePlane of this segment
        const GlidePlaneSharedPtrType pGlidePlane;
        
        
        static double quadPerLength;
        size_t qOrder;
        
        static double virtualSegmentDistance;
        
        
        //! Positions corrersponding to the quadrature points
        MatrixDimQorder rgauss;
        
        
        //        std::array<MatrixDim, qOrder> stressGauss;
        std::deque<MatrixDim> stressGauss;
        
        //! PK force corrersponding to the quadrature points
        MatrixDimQorder pkGauss;
        
    private:
        
        /**********************************************************************/
        MatrixNdof stiffness_integrand(const int& k) const
        { /*! @param[in] k the current quadrature point
           *  The stiffness matrix integrand evaluated at the k-th quadrature point.
           *	\f[
           *		\mathbf{K}^* = \mathbf{N}^T \mathbf{B} \mathbf{N} \frac{dl}{du}
           *	\f]
           */
            const MatrixDimNdof SFEx(SFgaussEx(k));
            //			return temp.transpose()*Material<Isotropic>::B*temp*jgauss(k);
            return SFEx.transpose()*SFEx*jgauss(k); // inverse mobility law
        }
        
        /**********************************************************************/
        MatrixDimNdof SFgaussEx(const int& k) const
        { /*! The MatrixDimNdof matrix of shape functions at the k-th quadrature point
           *  TO DO: EXPRESSION NEEDS TO BE GENERALIZED
           */
            return (MatrixDimNdof()<<I*SFgauss(k,0),I*SFgauss(k,1),I*SFgauss(k,2),I*SFgauss(k,3)).finished();
        }
        
        /**********************************************************************/
        VectorNdof PKintegrand(const int& k) const
        { /*! The force vector integrand evaluated at the k-th quadrature point.
           *  @param[in] k the current quadrature point
           */
            const VectorDim glideForce = pkGauss.col(k)-pkGauss.col(k).dot(this->glidePlaneNormal)*this->glidePlaneNormal;
            const double glideForceNorm(glideForce.norm());
            VectorDim vv=VectorDim::Zero();
            if(glideForceNorm>FLT_EPSILON)
            {
                double v =  Material<Isotropic>::velocity(stressGauss[k],Burgers,rlgauss.col(k),this->glidePlaneNormal);
                assert(v>= 0.0 && "Velocity must be a positive scalar");
                const bool useNonLinearVelocity=true;
                if(useNonLinearVelocity && v>FLT_EPSILON)
                {
                    v= 1.0-std::exp(-v);
                }
                                
                vv= v * glideForce/glideForceNorm;
            }
            //return temp.transpose()*pkGauss.col(k)*jgauss(k);
            return SFgaussEx(k).transpose()*vv*jgauss(k); // inverse mobility law
//            return SFgaussEx(k).transpose()*radiativeVel(pkGauss.col(k))*jgauss(k); // inverse mobility law
            //            return temp.transpose()*dm.getVelocity(stressGauss[k],rlgauss.col(k))*jgauss(k); // inverse mobility law
        }
        
//        /**********************************************************************/
//        VectorDim radiativeVel(const VectorDim& pkF) const
//        {
////            const VectorDim v0(Material<Isotropic>::Binv*pkF);
//            const VectorDim v0(Material<Isotropic>::velocity(pkF));
//            const double v0N(v0.norm());
//            const double csf(0.7*Material<Isotropic>::cs);
//            return (v0N>FLT_EPSILON)? csf*(1.0-std::exp(-v0N/csf))*v0/v0N : v0;
//        }
        
        
        //#ifdef UserStressFile
        //#include UserStressFile
        //#endif
        
        /******************************************************************/
    public: // member functions
        /******************************************************************/
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        
        DislocationSegment(const std::pair<NodeType*,NodeType*> nodePair, const LatticeVectorType& Fin) :
        /* base class initialization */ PlanarSegmentType(nodePair.first->get_L(),nodePair.second->get_L(),Fin),
        /* base class initialization */ SegmentBaseType(nodePair,Fin),
        /* init list       */ Burgers(this->flow.cartesian() * Material<Isotropic>::b),
        /* init list       */ isSessile(this->flow.dot(this->glidePlane.n)!=0),
        /* init list       */ conjugatePlaneNormals(CrystalOrientation<dim>::conjugatePlaneNormal(this->flow,this->glidePlane.n)),
        //        /* init list       */ boundaryLoopNormal(this->glidePlaneNormal),
        /* init list       */ pGlidePlane(this->findExistingGlidePlane(this->glidePlaneNormal,this->source->get_P().dot(this->glidePlaneNormal))), // change this
        /* init list       */ qOrder(QuadPowDynamicType::lowerOrder(quadPerLength*this->chord().norm()))
        {/*! Constructor with pointers to source and sink, and flow
          *  @param[in] NodePair_in the pair of source and sink pointers
          *  @param[in] Flow_in the input flow
          */
            
            DislocationEnergyRules<dim>::template findEdgeConfiguration<NodeType>(*this->source); // This should not be called in edge expansion or contraction
            this->source->make_T();
            
            DislocationEnergyRules<dim>::template findEdgeConfiguration<NodeType>(*this->sink); // This should not be called in edge expansion or contraction
            this->sink->make_T();
            
            assert(this->flow.squaredNorm()>0.0);
            pGlidePlane->addToGLidePlane(this);
            pkGauss.setZero(dim,qOrder); // necessary if this is not assembled
        }
        
        /* Constructor from EdgeExpansion) ************************************/
        DislocationSegment(const std::pair<NodeType*,NodeType*> nodePair, const EdgeRef<LinkType>& ee) :
//        /* base class initialization */ PlanarSegmentType(nodePair.first->get_L(),nodePair.second->get_L(),ee.E.flow),
        /* base class initialization */ PlanarSegmentType(nodePair.first->get_L(),nodePair.second->get_L(),ee.E.flow,ee),
        /* base class initialization */ SegmentBaseType::SplineSegmentBase(nodePair,ee),
        /* init list       */ Burgers(this->flow.cartesian() * Material<Isotropic>::b),
        /* init list       */ isSessile(this->flow.dot(this->glidePlane.n)!=0),
        /* init list       */ conjugatePlaneNormals(CrystalOrientation<dim>::conjugatePlaneNormal(this->flow,this->glidePlane.n)),
        //        /* init list       */ boundaryLoopNormal(this->glidePlaneNormal),
        /* init list       */ pGlidePlane(this->findExistingGlidePlane(this->glidePlaneNormal,this->source->get_P().dot(this->glidePlaneNormal))), // change this
        /* init list       */ qOrder(QuadPowDynamicType::lowerOrder(quadPerLength*this->chord().norm()))
        {/*! Constructor with pointers to source and sink, and EdgeRef
          *  @param[in] NodePair_in the pair of source and sink pointers
          *  @param[in] ee the expanding edge
          */
            
            DislocationEnergyRules<dim>::template findEdgeConfiguration<NodeType>(*this->source); // This should not be called in edge expansion or contraction
            this->source->make_T();
            
            DislocationEnergyRules<dim>::template findEdgeConfiguration<NodeType>(*this->sink); // This should not be called in edge expansion or contraction
            this->sink->make_T();
            
            assert(this->flow.squaredNorm()>0.0);
            pGlidePlane->addToGLidePlane(this);
            pkGauss.setZero(dim,qOrder); // necessary if this is not assembled
        }
        
        /* Destructor *********************************************************/
        ~DislocationSegment()
        {/*! Destructor
          */
            //! Add this to static VirtualBoundarySlipContainer
            //			if (shared.boundary_type==softBoundary && shared.use_bvp)
            //            {
            //				if(is_boundarySegment() && this->chord().norm()>FLT_EPSILON)
            //                {
            //                    shared.vbsc.add(*this);
            //				}
            //			}
            //! Removes this from *pGlidePlane
            pGlidePlane->removeFromGlidePlane(this);
            //            quadratureParticleContainer.clear();
        }
        
        
        
        /**********************************************************************/
        void updateQuadraturePoints(ParticleSystem<DislocationParticleType>& particleSystem)
        {/*! @param[in] particleSystem the ParticleSystem of DislocationParticle
          *  Computes all geometric properties at the k-th quadrature point
          */
            
            
            quadratureParticleContainer.clear();
            qOrder=QuadPowDynamicType::lowerOrder(quadPerLength*this->chord().norm());
            quadratureParticleContainer.reserve(qOrder);
            
            const MatrixNcoeff  SFCH(this->get_SFCH());
            const MatrixNcoeffDim qH(this->get_qH());
            SFgauss.setZero(qOrder,Ncoeff);
            rgauss.setZero(dim,qOrder);
            rugauss.setZero(dim,qOrder);
            jgauss.setZero(qOrder);
            rlgauss.setZero(dim,qOrder);
            pkGauss.setZero(dim,qOrder);
            
            // Compute geometric quantities
            for (unsigned int k=0;k<qOrder;++k)
            {
                //                SFgauss.row(k)=QuadPowType::uPow.row(k)*SFCH; // WHY ARE WE LOOPING TO DO THIS MATRIX MULTIPLICATION???? THIS SHOULD BE STORED IN QUADRATURE PARTICLE
                SFgauss.row(k)=QuadPowDynamicType::uPow(qOrder).row(k)*SFCH; // WHY ARE WE LOOPING TO DO THIS MATRIX MULTIPLICATION???? THIS SHOULD BE STORED IN QUADRATURE PARTICLE
                rgauss.col(k)=SFgauss.row(k)*qH; // WHY ARE WE LOOPING TO DO THIS MATRIX MULTIPLICATION???? THIS SHOULD BE STORED IN QUADRATURE PARTICLE
                //                rugauss.col(k)=QuadPowType::duPow.row(k)*SFCH.template block<Ncoeff-1,Ncoeff>(1,0)*qH; // WHY ARE WE LOOPING TO DO THIS MATRIX MULTIPLICATION???? THIS SHOULD BE STORED IN QUADRATURE PARTICLE
                rugauss.col(k)=QuadPowDynamicType::duPow(qOrder).row(k)*SFCH.template block<Ncoeff-1,Ncoeff>(1,0)*qH; // WHY ARE WE LOOPING TO DO THIS MATRIX MULTIPLICATION???? THIS SHOULD BE STORED IN QUADRATURE PARTICLE
                jgauss(k)=rugauss.col(k).norm();
                rlgauss.col(k)=rugauss.col(k)/jgauss(k);
            }
            
            if(!is_boundarySegment())
            {
                for (unsigned int k=0;k<qOrder;++k)
                {
                    quadratureParticleContainer.push_back(particleSystem.addParticle(rgauss.col(k),
                                                                                     this->source->sID,this->sink->sID,k,
                                                                                     rugauss.col(k),
                                                                                     Burgers,
                                                                                     QuadratureDynamicType::abscissa(qOrder,k),
                                                                                     QuadratureDynamicType::weight(qOrder,k),
                                                                                     true,true,  // stressSource enabled, stressField enabled,
                                                                                     true,true,  //   dispSource enabled,   dispField enabled,
                                                                                     true,true));//   enrgSource enabled,   enrgField enabled,
                }
                
            }
            else // bonudary segment
            {
                if(shared.use_bvp) // using FEM correction
                {
                    if(shared.use_virtualSegments)
                    {
                        for (unsigned int k=0;k<qOrder;++k)
                        {
                            quadratureParticleContainer.push_back(particleSystem.addParticle(rgauss.col(k),
                                                                                             this->source->sID,this->sink->sID,k,
                                                                                             rugauss.col(k),Burgers,
                                                                                             QuadratureDynamicType::abscissa(qOrder,k),
                                                                                             QuadratureDynamicType::weight(qOrder,k),
                                                                                             false,true,  // stressSource disabled, stressField enabled,
                                                                                             false,true,  //   dispSource disabled,   dispField enabled,
                                                                                             false,true));//   enrgSource disabled,   enrgField enabled,
                        }
                        
                        const VectorDim P1(this->source->get_P());
                        const VectorDim P2(P1+this->source->bndNormal()*virtualSegmentDistance);
                        const VectorDim P3(this->sink->get_P());
                        const VectorDim P4(P3+this->sink->bndNormal()*virtualSegmentDistance);
                        
                        const size_t qOrder12=QuadPowDynamicType::lowerOrder(quadPerLength*virtualSegmentDistance);
                        
                        // Place Quadrature-particles on P1->P2
                        if(!this->source->isPureBoundaryNode())
                        {
                            const VectorDim d12=(P2-P1).normalized();
                            
                            for (unsigned int k=0;k<qOrder12;++k)
                            {
                                
                                particleSystem.addParticle(P1+d12*virtualSegmentDistance*QuadratureDynamicType::abscissa(qOrder12,k),
                                                           this->source->sID,this->sink->sID,qOrder+k,
                                                           d12*virtualSegmentDistance,
                                                           Burgers,
                                                           QuadratureDynamicType::abscissa(qOrder12,k),
                                                           QuadratureDynamicType::weight(qOrder12,k),
                                                           true,false,  // stressSource enabled, stressField disabled,
                                                           true,false,   //   dispSource  enabled,   dispField disabled,
                                                           true,false);
                            }
                            
                        }
                        
                        // Place Quadrature-particles on P2->P4
                        const double L24=(P4-P2).norm();
                        const size_t qOrder24=QuadPowDynamicType::lowerOrder(quadPerLength*L24);
                        const VectorDim d24=(P4-P2)/L24;
                        for (unsigned int k=0;k<qOrder24;++k)
                        {
                            particleSystem.addParticle(P2+d24*L24*QuadratureDynamicType::abscissa(qOrder24,k),
                                                       this->source->sID,this->sink->sID,qOrder+qOrder12+k,
                                                       d24*L24,
                                                       Burgers,
                                                       QuadratureDynamicType::abscissa(qOrder24,k),
                                                       QuadratureDynamicType::weight(qOrder24,k),
                                                       true,false,  // stressSource enabled, stressField disabled,
                                                       true,false,   //   dispSource  enabled,   dispField disabled,
                                                       true,false);
                        }
                        
                        // Place Quadrature-particles on P4->P3
                        if(!this->sink->isPureBoundaryNode())
                        {
                            const VectorDim d43=(P3-P4).normalized();
                            for (unsigned int k=0;k<qOrder12;++k)
                            {
                                
                                particleSystem.addParticle(P4+d43*virtualSegmentDistance*QuadratureDynamicType::abscissa(qOrder12,k),
                                                           this->source->sID,this->sink->sID,qOrder+qOrder12+qOrder24+k,
                                                           d43*virtualSegmentDistance,
                                                           Burgers,
                                                           QuadratureDynamicType::abscissa(qOrder12,k),
                                                           QuadratureDynamicType::weight(qOrder12,k),
                                                           true,false,  // stressSource enabled, stressField disabled,
                                                           true,false,   //   dispSource  enabled,   dispField disabled,
                                                           true,false);
                            }
                            
                        }
                        
                    }
                    else // bnd segment, with bvp, without virtual segments
                    {
                        for (unsigned int k=0;k<qOrder;++k)
                        {
                            quadratureParticleContainer.push_back(particleSystem.addParticle(rgauss.col(k),
                                                                                             this->source->sID,this->sink->sID,k,
                                                                                             rugauss.col(k),Burgers,
                                                                                             QuadratureDynamicType::abscissa(qOrder,k),
                                                                                             QuadratureDynamicType::weight(qOrder,k),
                                                                                             true,true,  // stressSource enabled, stressField enabled,
                                                                                             true,true,  //   dispSource enabled,   dispField enabled,
                                                                                             true,true));//   enrgSource enabled,   enrgField enabled,
                        }
                    }
                }
                else // bonudary segment without bvp, do not place quadrature particles
                {
                    
                }
                
            }
            
            
            
            //            const MatrixNcoeff  SFCH(this->get_SFCH());
            //            const MatrixNcoeffDim qH(this->get_qH());
            //            SFgauss.setZero(qOrder,Ncoeff);
            //            rgauss.setZero(dim,qOrder);
            //            rugauss.setZero(dim,qOrder);
            //            jgauss.setZero(qOrder);
            //            rlgauss.setZero(dim,qOrder);
            //            pkGauss.setZero(dim,qOrder);
            //
            //            for (unsigned int k=0;k<qOrder;++k)
            //            {
            //                //                SFgauss.row(k)=QuadPowType::uPow.row(k)*SFCH; // WHY ARE WE LOOPING TO DO THIS MATRIX MULTIPLICATION???? THIS SHOULD BE STORED IN QUADRATURE PARTICLE
            //                SFgauss.row(k)=QuadPowDynamicType::uPow(qOrder).row(k)*SFCH; // WHY ARE WE LOOPING TO DO THIS MATRIX MULTIPLICATION???? THIS SHOULD BE STORED IN QUADRATURE PARTICLE
            //                rgauss.col(k)=SFgauss.row(k)*qH; // WHY ARE WE LOOPING TO DO THIS MATRIX MULTIPLICATION???? THIS SHOULD BE STORED IN QUADRATURE PARTICLE
            //                //                rugauss.col(k)=QuadPowType::duPow.row(k)*SFCH.template block<Ncoeff-1,Ncoeff>(1,0)*qH; // WHY ARE WE LOOPING TO DO THIS MATRIX MULTIPLICATION???? THIS SHOULD BE STORED IN QUADRATURE PARTICLE
            //                rugauss.col(k)=QuadPowDynamicType::duPow(qOrder).row(k)*SFCH.template block<Ncoeff-1,Ncoeff>(1,0)*qH; // WHY ARE WE LOOPING TO DO THIS MATRIX MULTIPLICATION???? THIS SHOULD BE STORED IN QUADRATURE PARTICLE
            //                jgauss(k)=rugauss.col(k).norm();
            //                rlgauss.col(k)=rugauss.col(k)/jgauss(k);
            //
            //
            //                if (!is_boundarySegment())
            //                {// not a boundary segment: place gauss points as both source and field points
            //                    quadratureParticleContainer.push_back(particleSystem.addParticle(rgauss.col(k),
            //                                                                                     this->source->sID,this->sink->sID,k,
            //                                                                                     rugauss.col(k),Burgers,
            //                                                                                     QuadratureDynamicType::abscissa(qOrder,k),
            //                                                                                     QuadratureDynamicType::weight(qOrder,k),
            //                                                                                     true,true,  // stressSource enabled, stressField enabled,
            //                                                                                     true,true,  //   dispSource enabled,   dispField enabled,
            //                                                                                     true,true));//   enrgSource enabled,   enrgField enabled,
            //                }
            //                else
            //                {// a boundary segment
            //                    if(shared.use_bvp)
            //                    {
            //                        if(shared.use_virtualSegments)
            //                        {// a boundary segment with bvp and with virtual segments: place gauss points only as field
            //                            quadratureParticleContainer.push_back(particleSystem.addParticle(rgauss.col(k),
            //                                                                                             this->source->sID,this->sink->sID,k,
            //                                                                                             rugauss.col(k),Burgers,
            //                                                                                             QuadratureDynamicType::abscissa(qOrder,k),
            //                                                                                             QuadratureDynamicType::weight(qOrder,k),
            //                                                                                             false,true,  // stressSource disabled, stressField enabled,
            //                                                                                             false,true,   //   dispSource  enabled,   dispField enabled,
            //                                                                                             false,true)); //   enrgSource  enabled,   enrgField enabled,
            //                            // first triangle is P1->P2->P3, second triangle is P2->P4->P3
            //
            //
            //                        }
            //                        else
            //                        {// a boundary segment with bvp and without virtual segments: place gauss points as both source and field points
            //                            // disable stress source, enable displacement and energy source
            //                            quadratureParticleContainer.push_back(particleSystem.addParticle(rgauss.col(k),
            //                                                                                             this->source->sID,this->sink->sID,k,
            //                                                                                             rugauss.col(k),Burgers,
            //                                                                                             QuadratureDynamicType::abscissa(qOrder,k),
            //                                                                                             QuadratureDynamicType::weight(qOrder,k),
            //                                                                                             true,true,  // stressSource enabled, stressField enabled,
            //                                                                                             true,true,  //   dispSource enabled,   dispField enabled,
            //                                                                                             true,true));//   enrgSource enabled,   enrgField enabled,
            //                        }
            //                    }
            //                    else
            //                    {// a boundary segment without bvp: don't place gauss points
            //
            //                    }
            //                }
            //            }
            
            
        }
        
        
        
        /**********************************************************************/
        MatrixDim stressAtQuadrature(const size_t & k) const
        {/*!@param[in] k the k-th quandrature point
          *\returns the stress field at the k-th quandrature point
          */
            
            MatrixDim temp(quadratureParticleContainer[k]->stress()+shared.extStressController.externalStress());
            if(shared.use_bvp)
            {
                //                temp += shared.bvpSolver.stress(quadratureParticleContainer[k]->P,this->source->includingSimplex());
                temp += shared.bvpSolver.stress(rgauss.col(k),this->source->includingSimplex());
                
            }
                        
            return temp;
            //            //            return (shared.use_bvp) ? ((quadratureParticleContainer[k]->stress(this->source->bvpStress,this->sink->bvpStress)+shared.vbsc.stress(quadratureParticleContainer[k]->P)+shared.externalStress)*Burgers).cross(rlgauss.col(k))
            //            return (shared.use_bvp) ? (quadratureParticleContainer[k]->stress()+shared.externalStress+shared.bvpSolver.stress(quadratureParticleContainer[k]->P,this->source->includingSimplex()))
            //			/*                   */ : (quadratureParticleContainer[k]->stress()+shared.externalStress);
            
        }
        
        /**********************************************************************/
        MatrixDimQorder get_pkGauss() const
        {/*!\returns the matrix of PK force at the quadrature points
          */
            return pkGauss;
        }
        
        /**********************************************************************/
        void assemble()
        {/*!Computes the following:
          * - edge stiffness matrix Kqq
          * - nodal force vector Fq
          * - edge-to-component matrix Mseg
          */
            
            //! 1- Compute and store stress and PK-force at quadrature points
            stressGauss.clear();
            if(!shared.use_bvp && is_boundarySegment())
            {
                pkGauss.setZero(dim,qOrder);
            }
            else
            {
                for (unsigned int k=0;k<qOrder;++k)
                {
                    stressGauss.push_back(stressAtQuadrature(k));
                    pkGauss.col(k)=(stressGauss[k]*Burgers).cross(rlgauss.col(k));
                }
            }
            
            /*! 2- Assemble the force vector of this segment
             *	\f[
             *		\mathbf{K} = int_0^1 \mathbf{K}^*(u) du
             *	\f]
             */
            Fq.setZero();
            //            Quadrature<1,qOrder,QuadratureRule>::integrate(this,Fq,&LinkType::PKintegrand);
            QuadratureDynamicType::integrate(qOrder,this,Fq,&LinkType::PKintegrand);
            
            /*! 3- Assembles the stiffness matrix of this segment.
             *	\f[
             *		\mathbf{K} = int_0^1 \mathbf{K}^*(u) du
             *	\f]
             */
            Kqq.setZero();
            //            Quadrature<1,qOrder,QuadratureRule>::integrate(this,Kqq,&LinkType::stiffness_integrand);
            QuadratureDynamicType::integrate(qOrder,this,Kqq,&LinkType::stiffness_integrand);
            
            //            //! 3-
            //            ortC.setZero();
            //
            //            Quadrature<1,qOrder,QuadratureRule>::integrate(this,ortC,&LinkType::ortC_integrand);
            
            
            
            //! 4-
            
            //            std::set<size_t> segmentDOFs;
            segmentDOFs.clear();
            
            const Eigen::VectorXi sourceDOFs(this->source->dofID());
            for(int k=0;k<sourceDOFs.rows();++k)
            {
                segmentDOFs.insert(sourceDOFs(k));
            }
            
            const Eigen::VectorXi   sinkDOFs(this->sink->dofID());
            for(int k=0;k<sinkDOFs.rows();++k)
            {
                segmentDOFs.insert(sinkDOFs(k));
            }
            
            //const size_t N(segmentDOFs.size());
            
            // Eigen::Matrix<double, Ndof, Eigen::Dynamic> Mseg(Eigen::Matrix<double, Ndof, Eigen::Dynamic>::Zero(Ndof,N));
            Mseg.setZero(Ndof,segmentDOFs.size());
            
            Eigen::Matrix<double, Ndof/2, Eigen::Dynamic> Mso(this->source->W2H());
            //            Eigen::Matrix<double, Ndof/2, Eigen::Dynamic> Mso(this->source->W2Ht());
            Mso.block(dim,0,dim,Mso.cols())*=this->sourceTfactor;
            
            Eigen::Matrix<double, Ndof/2, Eigen::Dynamic> Msi(this->sink->W2H());
            //            Eigen::Matrix<double, Ndof/2, Eigen::Dynamic> Msi(this->sink->W2Ht());
            Msi.block(dim,0,dim,Msi.cols())*=-this->sinkTfactor;
            
            
            for (int k=0;k<Mso.cols();++k)
            {
                const std::set<size_t>::const_iterator f(segmentDOFs.find(sourceDOFs(k)));
                assert(f!=segmentDOFs.end());
                unsigned int curCol(std::distance(segmentDOFs.begin(),f));
                Mseg.template block<Ndof/2,1>(0,curCol)=Mso.col(k);
            }
            
            for (int k=0;k<Msi.cols();++k)
            {
                const std::set<size_t>::const_iterator f(segmentDOFs.find(sinkDOFs(k)));
                assert(f!=segmentDOFs.end());
                unsigned int curCol(std::distance(segmentDOFs.begin(),f));
                Mseg.template block<Ndof/2,1>(Ndof/2,curCol)=Msi.col(k);
            }
        }
        
        
        //        Eigen::Matrix<double,1,Ndof> ortC_integrand(const int& k) const {
        //            return rugauss.col(k).transpose()*SFgaussEx(k)*Quadrature<1,qOrder,QuadratureRule>::abscissa(k)*(1.0-Quadrature<1,qOrder,QuadratureRule>::abscissa(k));
        //            //                                              *Quadrature<1,qOrder,QuadratureRule>::abscissa(k)*(1.0-Quadrature<1,qOrder,QuadratureRule>::abscissa(k));
        //            //return rugauss.col(k).transpose()*SFgaussEx(k);
        //        }
        
        
        
        /**********************************************************************/
        void addToGlobalAssembly(std::deque<Eigen::Triplet<double> >& kqqT,
                                 Eigen::VectorXd& FQ) const
        {/*!\param[in] kqqT the stiffness matrix of the network component
          * \param[in] FQ the force vector of the network component
          */
            
            const Eigen::MatrixXd tempKqq(Mseg.transpose()*Kqq*Mseg); // Create the temporaty stiffness matrix and push into triplets
            for (unsigned int i=0;i<segmentDOFs.size();++i)
            {
                std::set<size_t>::const_iterator iterI(segmentDOFs.begin());
                std::advance(iterI,i);
                for (unsigned int j=0;j<segmentDOFs.size();++j)
                {
                    std::set<size_t>::const_iterator iterJ(segmentDOFs.begin());
                    std::advance(iterJ,j);
                    if (std::fabs(tempKqq(i,j))>FLT_EPSILON)
                    {
                        kqqT.push_back(Eigen::Triplet<double>(*iterI,*iterJ,tempKqq(i,j)));
                    }
                }
            }
            
            const Eigen::VectorXd tempFq(Mseg.transpose()*Fq); // Create temporary force vector and add to global FQ
            for (unsigned int i=0;i<segmentDOFs.size();++i)
            {
                std::set<size_t>::const_iterator iterI(segmentDOFs.begin());
                std::advance(iterI,i);
                FQ(*iterI)+=tempFq(i);
            }
        }
        
        /**********************************************************************/
        const MatrixNdof& get_Kqq() const
        {/*!\returns the stiffness matrix of this segment
          */
            return Kqq;
        }
        
        /**********************************************************************/
        const VectorNdof& get_Fq() const
        {/*!\returns the nodal force vector for this segment
          */
            return Fq;
        }
        
        /**********************************************************************/
        MatrixDim plasticDistortionRate() const
        {/*!\returns the plastic strain rate generated by *this segment:
          *	\f[
          *		\beta^P_{ij}=\int_{\mathbb{R}^3}\int d\ell dV=-b_i\int\epsilon_{jkl}w_kd\ell_l
          *	\f]
          */
            //			VectorDim temp(VectorDim::Zero());
            //			Quadrature<1,qOrder,QuadratureRule>::integrate(this,temp,&LinkType::plasticStrainRateIntegrand);
            
            //\todo this integral should be calculated using shape functions
            
            const VectorDim V((this->source->get_V().template segment<dim>(0)+this->sink->get_V().template segment<dim>(0))*0.5);
            return -Burgers*V.cross(this->chord()).transpose()*(!is_boundarySegment());
        }
        
        /**********************************************************************/
        MatrixDim plasticStrainRate() const
        {
            const VectorDim temp(plasticDistortionRate());
            return (temp+temp.transpose())*0.5;
        }
        
        /**********************************************************************/
        bool is_boundarySegment() const
        {
            return (this->source->isBoundaryNode() && this->sink->isBoundaryNode() );
        }
        
        /**********************************************************************/
        Eigen::Matrix<double,dim-1,Ncoeff> hermiteLocalCoefficient() const
        {
            const MatrixDim G2L(DislocationLocalReference<dim>::global2local(this->chord(),this->glidePlaneNormal));
            Eigen::Matrix<double,dim-1,Ncoeff> HrCf = Eigen::Matrix<double,dim-1,Ncoeff>::Zero();
            HrCf.col(1)= (G2L*this->sourceT()*this->chordParametricLength()).template segment<dim-1>(0);
            HrCf.col(2)= (G2L*(this->sink->get_P()-this->source->get_P())).template segment<dim-1>(0);
            HrCf.col(3)= (G2L*this->sinkT()*this->chordParametricLength()).template segment<dim-1>(0);
            return HrCf;
        }
        
        /**********************************************************************/
        Eigen::Matrix<double,dim-1,Ncoeff> polynomialLocalCoeff() const //change name polynomialCoeff
        {
            return Coeff2Hermite<pOrder>::template h2c<dim-1>(hermiteLocalCoefficient());
        }
        
        /**********************************************************************/
        VectorDim integralPK() const
        {/*!\returns The integral of the PK force over this segment.
          *!\todo FINISH HERE
          */
            return pkGauss.col(qOrder/2)*this->chord().norm();
        }
        

        
        
        /**********************************************************************/
        int lineTriangleIntersection(const VectorDim& X,const VectorDim& s,
                                     const VectorDim& P1,const VectorDim& P2,const VectorDim& P3) const
        {
            
            int temp=0;
            
            const VectorDim v12(P2-P1);
            const VectorDim v13(P3-P1);
            
            VectorDim n(v12.cross(v13)); // right-handed norm to triangle P1->P2->P3
            const double nNorm(n.norm());
            assert(nNorm>FLT_EPSILON && "n has zero norm");
            n/=nNorm; // normalize
            const double nDots(n.dot(s));
            
            if(std::fabs(nDots)>FLT_EPSILON) // s in not parallel to triangle
            {
                MatrixDim M(MatrixDim::Zero());
                M.col(0)=v12;
                M.col(1)=v13;
                M.col(2)=-s;
                const VectorDim a=M.inverse()*(X-P1);
                
                if(   a(0)>=0.0 && a(0)<=1.0      // intersection with trignale exists
                   && a(1)>=0.0 && a(0)+a(1)<=1.0 // intersection with trignale exists
                   && a(2)>=0.0) //  intersection along positive S-vector
                {
                    if(nDots>=0.0)
                    {
                        temp=-1;
                    }
                    else
                    {
                        temp=+1;
                    }
                }
            }
//            else // s is parallel to triangle
//            {
//                assert(0 && "IMPLEMENT INTERSECTION AT INFINITY");
////                if((P1-X).dot(n)>=0.0) // X is below the trianlge. Positive intersection at infinity
////                {
////                    temp=-1;
////                }
////                else // X is above the trianlge. Negative intersection at infinity
////                {
////                    temp=+1;
////                }
//            }
            
            return temp;
        }
        
        /**********************************************************************/
        void addToSolidAngleJump(const VectorDim& Pf, const VectorDim& Sf, VectorDim& dispJump) const
        {
            
            if(is_boundarySegment() && shared.use_virtualSegments)
            {
                
                // first triangle is P1->P2->P3, second triangle is P2->P4->P3
                const VectorDim P1(this->source->get_P());
                const VectorDim P2(P1+this->source->bndNormal()*virtualSegmentDistance);
                const VectorDim P3(this->sink->get_P());
                const VectorDim P4(P3+this->sink->bndNormal()*virtualSegmentDistance);
                
                dispJump += Burgers*lineTriangleIntersection(Pf,Sf,P1,P2,P3);
                dispJump += Burgers*lineTriangleIntersection(Pf,Sf,P2,P4,P3);
                
                
//                if (dispJump.norm()>0.0)
//                {
//                    std::cout<<this->source->sID<<"->"<<this->sink->sID<<std::endl;
//                    std::cout<<Pf.transpose()<<std::endl;
//                    std::cout<<Sf.transpose()<<std::endl;
//                    std::cout<<P1.transpose()<<std::endl;
//                    std::cout<<P2.transpose()<<std::endl;
//                    std::cout<<P3.transpose()<<std::endl;
//                    std::cout<<P4.transpose()<<std::endl;
//                }
                
            }
        }
        
        /**********************************************************************/
        double arcLength() const
        {
            return SegmentBaseType::template arcLength<16,QuadratureRule>();
        }
        
        /**********************************************************************/
        VectorDim velocity(const double& u) const
        {
            return this->source->get_V().template segment<dim>(0)*(1.0-u)+this->sink->get_V().template segment<dim>(0)*u;
        }
        
        /*************************************************************/
        VectorDim integratedVelocity() const
        {
            VectorDim temp(VectorDim::Zero());
            Quadrature<1,16,QuadratureRule>::integrate(this,temp,&Derived::rm_integrand);
            return temp;
        }
        
        VectorDim integratedVelocityKernel(const double& u) const
        {
            return velocity(u)*this->get_j(u);
        }

        
        /**********************************************************************/
        template <class T>
        friend T& operator << (T& os, const Derived& ds)
        {
            os  << ds.source->sID<<"\t"<< ds.sink->sID<<"\t"
            /**/<< std::setprecision(15)<<std::scientific<<ds.Burgers.transpose()<<"\t"
            /**/<< std::setprecision(15)<<std::scientific<<ds.glidePlaneNormal.transpose()<<"\t"
            /**/<< ds.sourceTfactor<<"\t"
            /**/<< ds.sinkTfactor<<"\t"
            /**/<< ds.pSN()->sID;
            return os;
        }
        
    };
    
    // Static Data
    template <short unsigned int dim, short unsigned int corder, typename InterpolationType,
    /*	   */ template <short unsigned int, size_t> class QuadratureRule>
    const Eigen::Matrix<double,dim,dim> DislocationSegment<dim,corder,InterpolationType,QuadratureRule>::I=Eigen::Matrix<double,dim,dim>::Identity();
    
    template <short unsigned int dim, short unsigned int corder, typename InterpolationType,
    /*	   */ template <short unsigned int, size_t> class QuadratureRule>
    double DislocationSegment<dim,corder,InterpolationType,QuadratureRule>::quadPerLength=0.2;
    
    template <short unsigned int dim, short unsigned int corder, typename InterpolationType,
    /*	   */ template <short unsigned int, size_t> class QuadratureRule>
    double DislocationSegment<dim,corder,InterpolationType,QuadratureRule>::virtualSegmentDistance=200.0;
    
    
} // namespace model
#endif

