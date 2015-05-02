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
#include <model/Utilities/CompareVectorsByComponent.h>

//#include <model/DislocationDynamics/NearestNeighbor/DislocationQuadratureParticle.h>
#include <model/DislocationDynamics/NearestNeighbor/DislocationParticle.h>
#include <model/ParticleInteraction/ParticleSystem.h>

#include <model/DislocationDynamics/DislocationLocalReference.h>
#include <model/DislocationDynamics/Junctions/DislocationSegmentIntersection.h>

//#include <model/DislocationDynamics/CrossSlip/CrossSlipSegment.h>
#include <model/DislocationDynamics/DislocationMobility.h>
#include <model/Utilities/UniqueOutputFile.h>


//#include <model/DislocationDynamics/BVP/VirtualBoundarySlipContainer.h>


namespace model
{
    
    
    template <short unsigned int dim>
    struct PlanarDislocationSegment
    {/*! Class template used ad a base of DislocationSegment used to
      * initialize glidePlaneNormal and sessilePlaneNormal before the base
      * NetworkLink.
      */

        typedef Eigen::Matrix<double,dim,1> VectorDim;
//        typedef ReciprocalLatticeDirection<dim> ReciprocalLatticeDirectionType;
        typedef LatticeVector<dim> LatticeVectorType;

        //! The glide plane unit normal vector
//        const LatticePlaneBase& glidePlaneReciprocal;
        const LatticePlane   glidePlane;
        const LatticePlane sessilePlane;

        VectorDim   glidePlaneNormal;
         VectorDim sessilePlaneNormal;

        
        PlanarDislocationSegment(const LatticeVectorType& sourceL,
                                 const LatticeVectorType& sinkL,
                                 const LatticeVectorType& Burgers) :
        /* init list       */   glidePlane(sourceL,CrystalOrientation<dim>::find_glidePlane(sinkL-sourceL,Burgers)),
        /* init list       */ sessilePlane(sourceL,CrystalOrientation<dim>::find_sessilePlane(sinkL-sourceL,Burgers)),
        /* init list       */ glidePlaneNormal(glidePlane.n.cartesian().normalized()),
        /* init list       */ sessilePlaneNormal(sessilePlane.n.cartesian().normalized())
        {
        
        }
        
//        template<typename LinkType>
//        PlanarDislocationSegment(const ExpandingEdge<LinkType>& ee) :
//        /* init list       */ glidePlaneNormal(ee.E.glidePlaneNormal),
//        /* init list       */ sessilePlaneNormal(ee.E.sessilePlaneNormal)
//        {
//            THIS DOES NOT WORK FOR CROSS SLIP
//        }
    };
    
    template <short unsigned int _dim, short unsigned int corder, typename InterpolationType,
    /*	   */ short unsigned int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule>
    class DislocationSegment : public PlanarDislocationSegment<_dim>,
    /*	                    */ public SplineSegmentBase<DislocationSegment<_dim,corder,InterpolationType,qOrder,QuadratureRule>,
    /*                                              */ _dim, corder>,
    /*	                    */ public GlidePlaneObserver<DislocationSegment<_dim,corder,InterpolationType,qOrder,QuadratureRule> >
    {
        
        
    public:
        
        enum{dim=_dim}; // make dim available outside class
        
        
        typedef DislocationSegment<dim,corder,InterpolationType,qOrder,QuadratureRule> Derived; 		// Define "Derived" so that NetworkTypedefs.h can be used
#include <model/Network/NetworkTypedefs.h>
#include <model/Geometry/Splines/SplineEnums.h>
        
        typedef PlanarDislocationSegment<dim> PlanarSegmentType;
        typedef SplineSegmentBase<Derived,dim,corder> SegmentBaseType;
        typedef std::map<size_t,LinkType* const> AddressMapType;
        typedef typename AddressMapType::iterator AddressMapIteratorType;
        typedef Eigen::Matrix<double,dim,qOrder>	MatrixDimQorder;
        typedef Eigen::Matrix<double, 1, qOrder>	VectorQorder;
        typedef QuadPow<Ncoeff-1,qOrder,QuadratureRule> QuadPowType;
        
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
        Eigen::Matrix<double,qOrder,Ncoeff> SFgauss;
        
        
        //enum {Nslips=MaterialType::Nslips};
        
        std::set<size_t> segmentDOFs;
        Eigen::Matrix<double, Ndof, Eigen::Dynamic> Mseg;
        
        //! The static MaterialType material
        //static MaterialType material;
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
        
        
//        static UniqueOutputFile<'K'> k_file;
        
        /******************************************************************/
    public: //  data members
        /******************************************************************/
        
        QuadratureParticleContainerType quadratureParticleContainer;
        
        
        //! The Burgers vector
//        const VectorDim Burgers;
        const VectorDim Burgers;
        const bool isSessile;
        const std::deque<const LatticePlaneBase*> conjugatePlaneNormals;
        
//        const LatticeVectorType& Burgers;
        
        //! The glide plane unit normal vector
//        const VectorDim   glidePlaneNormal;
//        
//        const VectorDim sessilePlaneNormal;
        
        VectorDim boundaryLoopNormal;

        
        //		static double coreLsquared; //!todo remove this
        
        DislocationSharedObjects<Derived> shared;
        
        
        //! A shared pointer to the GlidePlane of this segment
        const GlidePlaneSharedPtrType pGlidePlane;
        
        const DislocationMobility<dim> dm;
        
        //! Positions corrersponding to the quadrature points
        MatrixDimQorder rgauss;
        
        
        std::array<MatrixDim, qOrder> stressGauss;
        
        //! PK force corrersponding to the quadrature points
        MatrixDimQorder pkGauss;
        
    private:
        

        
        /* stiffness_integrand ************************************************/
        MatrixNdof stiffness_integrand(const int& k) const
        { /*! @param[in] k the current quadrature point
           *  The stiffness matrix integrand evaluated at the k-th quadrature point.
           *	\f[
           *		\mathbf{K}^* = \mathbf{N}^T \mathbf{B} \mathbf{N} \frac{dl}{du}
           *	\f]
           */
            MatrixDimNdof temp(SFgaussEx(k));
            //			return temp.transpose()*Material<Isotropic>::B*temp*jgauss(k);
            return temp.transpose()*temp*jgauss(k); // inverse mobility law
        }
        
        /* SFgaussEx **********************************************************/
        MatrixDimNdof SFgaussEx(const int& k) const
        { /*! The MatrixDimNdof matrix of shape functions at the k-th quadrature point
           *  TO DO: EXPRESSION NEEDS TO BE GENERALIZED
           */
            return (MatrixDimNdof()<<I*SFgauss(k,0),I*SFgauss(k,1),I*SFgauss(k,2),I*SFgauss(k,3)).finished();
        }
        
        
        
        /* PKintegrand ********************************************************/
        VectorNdof PKintegrand(const int& k) const
        { /*! The force vector integrand evaluated at the k-th quadrature point.
           *  @param[in] k the current quadrature point
           */
            MatrixDimNdof temp(SFgaussEx(k));
            //return temp.transpose()*pkGauss.col(k)*jgauss(k);
            return temp.transpose()*radiativeVel(pkGauss.col(k))*jgauss(k); // inverse mobility law
            //            return temp.transpose()*dm.getVelocity(stressGauss[k],rlgauss.col(k))*jgauss(k); // inverse mobility law
        }
        
        /**********************************************************************/
        VectorDim radiativeVel(const VectorDim& pkF) const
        {
            const VectorDim v0(Material<Isotropic>::Binv*pkF);
            const double v0N(v0.norm());
            const double csf(0.7*Material<Isotropic>::cs);
            return (v0N>FLT_EPSILON)? csf*(1.0-std::exp(-v0N/csf))*v0/v0N : v0;
        }
        
        
        //#ifdef UserStressFile
        //#include UserStressFile
        //#endif
        
        /******************************************************************/
    public: // member functions
        /******************************************************************/
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        
        /* Constructor with Nodes and FLow ************************************/
//        DislocationSegment(const std::pair<NodeType*,NodeType*> nodePair, const VectorDim & Fin) :
        DislocationSegment(const std::pair<NodeType*,NodeType*> nodePair, const LatticeVectorType& Fin) :
//        /* base class initialization */ PlanarSegmentType(nodePair.second->get_P()-nodePair.first->get_P(),Fin),
        /* base class initialization */ PlanarSegmentType(nodePair.first->get_L(),nodePair.second->get_L(),Fin),
        /* base class initialization */ SegmentBaseType(nodePair,Fin),
        /* init list       */ Burgers(this->flow.cartesian() * Material<Isotropic>::b),
        /* init list       */ isSessile(this->flow.dot(this->glidePlane.n)!=0),
        /* init list       */ conjugatePlaneNormals(CrystalOrientation<dim>::conjugatePlaneNormal(this->flow,this->glidePlane.n)),
        /* init list       */ boundaryLoopNormal(this->glidePlaneNormal),
        /* init list       */ pGlidePlane(this->findExistingGlidePlane(this->glidePlaneNormal,this->source->get_P().dot(this->glidePlaneNormal))), // change this
        /* init list       */ dm(this->glidePlaneNormal,Burgers)
        {/*! Constructor with pointers to source and sink, and flow
          *  @param[in] NodePair_in the pair of source and sink pointers
          *  @param[in] Flow_in the input flow
          */
            
//            this->source->make_planeNormals(); // unnecessary
            DislocationEnergyRules<dim>::template findEdgeConfiguration<NodeType>(*this->source); // This should not be called in edge expansion or contraction
            this->source->make_T();
            
//            this->sink->make_planeNormals();
            DislocationEnergyRules<dim>::template findEdgeConfiguration<NodeType>(*this->sink); // This should not be called in edge expansion or contraction
            this->sink->make_T();
            
            assert(this->flow.squaredNorm()>0.0);
            pGlidePlane->addToGLidePlane(this);
            pkGauss.setZero(); // necessary if this is not assembled
        }
        
        /* Constructor from EdgeExpansion) ************************************/
        DislocationSegment(const std::pair<NodeType*,NodeType*> nodePair, const ExpandingEdge<LinkType>& ee) :
//        /* base class initialization */ PlanarSegmentType(ee),
//        /* base class initialization */ PlanarSegmentType(nodePair.second->get_P()-nodePair.first->get_P(),ee.E.flow),
        /* base class initialization */ PlanarSegmentType(nodePair.first->get_L(),nodePair.second->get_L(),ee.E.flow),
        /* base class initialization */ SegmentBaseType::SplineSegmentBase(nodePair,ee),
        /* init list       */ Burgers(this->flow.cartesian() * Material<Isotropic>::b),
        /* init list       */ isSessile(this->flow.dot(this->glidePlane.n)!=0),
                /* init list       */ conjugatePlaneNormals(CrystalOrientation<dim>::conjugatePlaneNormal(this->flow,this->glidePlane.n)),
        /* init list       */ boundaryLoopNormal(this->glidePlaneNormal),
        /* init list       */ pGlidePlane(this->findExistingGlidePlane(this->glidePlaneNormal,this->source->get_P().dot(this->glidePlaneNormal))), // change this
        /* init list       */ dm(this->glidePlaneNormal,Burgers)
        {/*! Constructor with pointers to source and sink, and ExpandingEdge
          *  @param[in] NodePair_in the pair of source and sink pointers
          *  @param[in] ee the expanding edge
          */
            
//            this->source->make_planeNormals();
            DislocationEnergyRules<dim>::template findEdgeConfiguration<NodeType>(*this->source); // This should not be called in edge expansion or contraction
            this->source->make_T();
            
//            this->sink->make_planeNormals();
            DislocationEnergyRules<dim>::template findEdgeConfiguration<NodeType>(*this->sink); // This should not be called in edge expansion or contraction
            this->sink->make_T();
            
            assert(this->flow.squaredNorm()>0.0);
            pGlidePlane->addToGLidePlane(this);
            pkGauss.setZero(); // necessary if this is not assembled
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
            quadratureParticleContainer.reserve(qOrder);
            
            const MatrixNcoeff  SFCH(this->get_SFCH());
            const MatrixNcoeffDim qH(this->get_qH());
            for (unsigned int k=0;k<qOrder;++k)
            {
                //                    const MatrixNcoeff  SFCH(this->get_SFCH());
                //                    const MatrixNcoeffDim qH(this->get_qH());
                SFgauss.row(k)=QuadPowType::uPow.row(k)*SFCH; // WHY ARE WE LOOPING TO DO THIS MATRIX MULTIPLICATION???? THIS SHOULD BE STORED IN QUADRATURE PARTICLE
                rgauss.col(k)=SFgauss.row(k)*qH; // WHY ARE WE LOOPING TO DO THIS MATRIX MULTIPLICATION???? THIS SHOULD BE STORED IN QUADRATURE PARTICLE
                rugauss.col(k)=QuadPowType::duPow.row(k)*SFCH.template block<Ncoeff-1,Ncoeff>(1,0)*qH; // WHY ARE WE LOOPING TO DO THIS MATRIX MULTIPLICATION???? THIS SHOULD BE STORED IN QUADRATURE PARTICLE
                jgauss(k)=rugauss.col(k).norm();
                rlgauss.col(k)=rugauss.col(k)/jgauss(k);
            }

            // Add the quadrature particles
            if (is_boundarySegment() && shared.use_bvp && shared.use_virtualSegments)
            {
                
                const VectorDim C((this->source->get_P()+this->sink->get_P())*0.5);
                const VectorDim vec(this->source->get_P()-C);
                const VectorDim n((this->source->bndNormal()+this->sink->bndNormal()).normalized());
                boundaryLoopNormal=this->glidePlaneNormal;
                const MatrixDim R(Eigen::AngleAxisd(0.5*M_PI,boundaryLoopNormal));
                const double dldu=vec.norm()*M_PI;
                //double sign_theta=1.0;
                if((R*vec).dot(n)<0.0)
                {
                    boundaryLoopNormal*=-1.0;
                }
                
                for (unsigned int k=0;k<qOrder;++k)
                {
                    const MatrixDim Rk(Eigen::AngleAxisd(Quadrature<1,qOrder,QuadratureRule>::abscissas(k)*M_PI,boundaryLoopNormal));
                    const VectorDim rgauss_temp=C+Rk*vec;
                    const VectorDim rugauss_temp=boundaryLoopNormal.cross(Rk*vec).normalized()*dldu;
                    
//                    k_file<<rgauss_temp.transpose()<<" "<<rugauss_temp.transpose()<<"\n";
                    
                    quadratureParticleContainer.push_back(particleSystem.addParticle(rgauss_temp,rugauss_temp,Burgers,
                                                                                     Quadrature<1,qOrder,QuadratureRule>::abscissas(k),
                                                                                     Quadrature<1,qOrder,QuadratureRule>::weights(k)));
                }
            }
            else
            { // segment inside mesh
                boundaryLoopNormal=this->glidePlaneNormal;

                for (unsigned int k=0;k<qOrder;++k)
                {
                    quadratureParticleContainer.push_back(particleSystem.addParticle(rgauss.col(k),rugauss.col(k),Burgers,
                                                                                     Quadrature<1,qOrder,QuadratureRule>::abscissas(k),
                                                                                     Quadrature<1,qOrder,QuadratureRule>::weights(k)));
                }
            }
            

            
        }
        
//        /**********************************************************************/
//        void updateQuadraturePoints(ParticleSystem<DislocationParticleType>& particleSystem)
//        {/*! @param[in] particleSystem the ParticleSystem of DislocationParticle
//          *  Computes all geometric properties at the k-th quadrature point
//          */
//            
//            quadratureParticleContainer.clear();
//            quadratureParticleContainer.reserve(qOrder);
//            
//            for (unsigned int k=0;k<qOrder;++k)
//            {
//                MatrixNcoeff  SFCH(this->get_SFCH());
//                MatrixNcoeffDim qH(this->get_qH());
//                SFgauss.row(k)=QuadPowType::uPow.row(k)*SFCH; // WHY ARE WE LOOPING TO DO THIS MATRIX MULTIPLICATION???? THIS SHOULD BE STORED IN QUADRATURE PARTICLE
//                rgauss.col(k)=SFgauss.row(k)*qH; // WHY ARE WE LOOPING TO DO THIS MATRIX MULTIPLICATION???? THIS SHOULD BE STORED IN QUADRATURE PARTICLE
//                rugauss.col(k)=QuadPowType::duPow.row(k)*SFCH.template block<Ncoeff-1,Ncoeff>(1,0)*qH; // WHY ARE WE LOOPING TO DO THIS MATRIX MULTIPLICATION???? THIS SHOULD BE STORED IN QUADRATURE PARTICLE
//                jgauss(k)=rugauss.col(k).norm();
//                rlgauss.col(k)=rugauss.col(k)/jgauss(k);
//                
//                
//                // Add the quadrature particles
//                if (is_boundarySegment() && shared.use_bvp && shared.use_virtualSegments)
//                {
//                    // if bvp with virtualSegments is used, add particles with zero Burgers for segments on the boundary
////                    quadratureParticleContainer.push_back(particleSystem.addParticle(rgauss.col(k),rugauss.col(k),zeroDim,
////                                                                                     Quadrature<1,qOrder,QuadratureRule>::abscissas(k),
////                                                                                     Quadrature<1,qOrder,QuadratureRule>::weights(k)));
//                    asseert(0 && "FINISH HERE. PLACE GAUSS POINTS IN CIRCLE");
//                    const VectorDim C((this->source->get_P()+this->sink->get_P())*0.5);
//                    const VectorDim n((this->source->boundaryNormal+this->sink->boundaryNormal).normalized());
//                    const MatrixDim R(Eigen::AngleAxisf(Quadrature<1,qOrder,QuadratureRule>::abscissas(k)*M_PI, glidePlaneNormal));
//                
//                }
//                else
//                {
//                    quadratureParticleContainer.push_back(particleSystem.addParticle(rgauss.col(k),rugauss.col(k),Burgers,
//                                                                                     Quadrature<1,qOrder,QuadratureRule>::abscissas(k),
//                                                                                     Quadrature<1,qOrder,QuadratureRule>::weights(k)));
//                }
//            }
//            
//        }
        
        //		/**********************************************************************/
        //		VectorDim pkForce(const size_t & k)
        //        {/*!@param[in] k the k-th quandrature point
        //          *\returns the PK force at the k-th quandrature point
        //          */
        //            //            return (shared.use_bvp) ? ((quadratureParticleContainer[k]->stress(this->source->bvpStress,this->sink->bvpStress)+shared.vbsc.stress(quadratureParticleContainer[k]->P)+shared.externalStress)*Burgers).cross(rlgauss.col(k))
        //            return (shared.use_bvp) ? ((quadratureParticleContainer[k]->stress()+shared.externalStress+shared.bvpSolver.stress(quadratureParticleContainer[k]->P,this->source->includingSimplex()))*Burgers).cross(rlgauss.col(k))
        //			/*                   */ : ((quadratureParticleContainer[k]->stress()+shared.externalStress)*Burgers).cross(rlgauss.col(k));
        //
        //        }
        
        /**********************************************************************/
        MatrixDim stressAtQuadrature(const size_t & k)
        {/*!@param[in] k the k-th quandrature point
          *\returns the PK force at the k-th quandrature point
          */
            
            MatrixDim temp(quadratureParticleContainer[k]->stress()+shared.externalStress);
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
        {/*!\returns the matrix of PK force at the quandrature points
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
            
            /*! 1- Assembles the force vector of this segment.
             *	\f[
             *		\mathbf{K} = int_0^1 \mathbf{K}^*(u) du
             *	\f]
             */
            //                pkGauss.setZero(); // not strictly necessary
            for (int k=0;k<qOrder;++k)
            {
                stressGauss[k]=stressAtQuadrature(k);
                //                pkGauss.col(k)=pkForce(k);
                pkGauss.col(k)=(stressGauss[k]*Burgers).cross(rlgauss.col(k));
            }
            Fq.setZero();
            Quadrature<1,qOrder,QuadratureRule>::integrate(this,Fq,&LinkType::PKintegrand);
            
            /*! 2- Assembles the stiffness matrix of this segment.
             *	\f[
             *		\mathbf{K} = int_0^1 \mathbf{K}^*(u) du
             *	\f]
             */
            Kqq.setZero();
            Quadrature<1,qOrder,QuadratureRule>::integrate(this,Kqq,&LinkType::stiffness_integrand);
            
            
            
            //            //! 3-
            //            ortC.setZero();
            //
            //            Quadrature<1,qOrder,QuadratureRule>::integrate(this,ortC,&LinkType::ortC_integrand);
            
            
            
            //! 4-
            
            //            std::set<size_t> segmentDOFs;
            segmentDOFs.clear();
            
            const Eigen::VectorXi sourceDOFs(this->source->dofID());
            for(int k=0;k<sourceDOFs.rows();++k){
                segmentDOFs.insert(sourceDOFs(k));
            }
            
            const Eigen::VectorXi   sinkDOFs(this->sink->dofID());
            for(int k=0;k<sinkDOFs.rows();++k){
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
//        void addToGlobalAssembly(std::vector<Eigen::Triplet<double> >& kqqT,  Eigen::VectorXd& FQ) const
        void addToGlobalAssembly(std::deque<Eigen::Triplet<double> >& kqqT,  Eigen::VectorXd& FQ) const
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
        {/*!\returns the edge stiffness matrix
          */
            return Kqq;
        }
        
        /**********************************************************************/
        const VectorNdof& get_Fq() const
        {/*!\returns the edge force vector
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
            return -Burgers*V.cross(this->chord()).transpose();
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
        
//        /**********************************************************************/
//        std::map<double,VectorDim> boundaryCollision() const
//        {/*! The collision points with the glidePlane boundary polygon.
//          */
//            std::map<double,VectorDim> temp;
//            
//            if (shared.boundary_type){
//                VectorDim origin = this->source->get_P();
//                
//                const Eigen::Matrix<double,dim-1,Ncoeff> C1L(polynomialLocalCoeff()); // the local polynomial coefficients of this
//                PlanarSplineImplicitization<pOrder> psi(Coeff2Hermite<pOrder>::template h2c<dim-1>(C1L));
//                
//                const MatrixDim G2L(DislocationLocalReference<dim>::global2local(this->chord(),this->glidePlaneNormal));
//                
//                for (int k=0; k<pGlidePlane->segmentMeshCollisionPairContainer.size();++k)
//                {
//                    Eigen::Matrix<double,dim-1,2> H2L;
//                    H2L.col(0)=(G2L*(pGlidePlane->segmentMeshCollisionPairContainer[k].first -origin)).template segment<dim-1>(0);
//                    H2L.col(1)=(G2L*(pGlidePlane->segmentMeshCollisionPairContainer[k].second-origin)).template segment<dim-1>(0);
//                    std::set<std::pair<double,double> > intersectinParameters = psi.template intersectWith<1>(Coeff2Hermite<1>::template h2c<dim-1>(H2L),FLT_EPSILON);
//                    for (std::set<std::pair<double,double> >::const_iterator iter=intersectinParameters.begin();iter!=intersectinParameters.end();++iter)
//                    {
//                        temp.insert(std::make_pair(iter->first,this->get_r(iter->first)));
//                    }
//                }
//            }
//            
//            return temp;
//        }
        
//        /*********************************************************************/
//        CrossSlipSegment<LinkType> crossSlipSegment(const double& sinThetaCrossSlipCr,const double& crossSlipLength) const
//        {
//            return CrossSlipSegment<LinkType>(*this,sinThetaCrossSlipCr,crossSlipLength);
//        }
        
//        /*********************************************************************/
//        std::deque<const LatticePlaneBase*> conjugatePlaneNormals() const
//        {
//            return CrystalOrientation<dim>::conjugatePlaneNormal(this->flow,this->glidePlane.n);
//        }
        
        /*********************************************************************/
        VectorDim integralPK() const
        {/*!\returns The integral of the PK force over this segment.
          *!\todo FINISH HERE
          */
            return pkGauss.col(qOrder/2)*this->chord().norm();
        }
        
        /*********************************************************************/
        void addToSolidAngleJump(const VectorDim& Pf, const VectorDim& Sf, VectorDim& dispJump) const
        {
            if(is_boundarySegment())
            {
                
                const double den(Sf.dot(boundaryLoopNormal));
                if(std::fabs(den)>FLT_EPSILON) // s direction intersects plane
                {
                    const double u=(this->source->get_P()-Pf).dot(boundaryLoopNormal)/den;
                    if(u>FLT_EPSILON) // intersection in the positive sense
                    {
                        const VectorDim P=Pf+u*Sf; // intersection point on plane
                        const VectorDim C=(this->source->get_P()+this->sink->get_P())*0.5; // center of segment
                        
                        const double R=(this->source->get_P()-this->sink->get_P()).norm();
                        const double PC=(P-C).norm();
                        
                        if(PC<=R) // intersection point is inside circle
                        {
                            const VectorDim n((this->source->bndNormal()+this->sink->bndNormal()).normalized());
                            const double PCn=(P-C).dot(n);
                            if(PCn>0.0) // external side of loop
                            {
                                if(den>0.0)
                                {
                                    dispJump -= Burgers;
                                }
                                else
                                {
                                    dispJump += Burgers;
                                }
                            }
                        }
                        
                    }
                }
                
            }
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
    /*	   */ short unsigned int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule>
    const Eigen::Matrix<double,dim,dim> DislocationSegment<dim,corder,InterpolationType,qOrder,QuadratureRule>::I=Eigen::Matrix<double,dim,dim>::Identity();

    
} // namespace model
#endif


//    template <short unsigned int dim, short unsigned int corder, typename InterpolationType,
//    /*	   */ short unsigned int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule>
//    const Eigen::Matrix<double,dim,1> DislocationSegment<dim,corder,InterpolationType,qOrder,QuadratureRule>::zeroDim=Eigen::Matrix<double,dim,1>::Zero();

//    template <short unsigned int dim, short unsigned int corder, typename InterpolationType,
//    /*	   */ short unsigned int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule>
//    UniqueOutputFile<'K'> DislocationSegment<dim,corder,InterpolationType,qOrder,QuadratureRule>::k_file;




//		/**********************************************************************/
//		VectorDim displacement(const VectorDim & Rfield, const VectorDim& S) const __attribute__ ((deprecated))
//        {/*!@param[in] Rfield the field point
//          * @param[in] S the unit vector necessary to compute displacement as line integral
//          * \returns the infinitesimal dispacement field generated by this segment at Rfield
//          *
//          * The return value is calculated according to:
//          *	\f[
//          *		\mathbf{u}(\mathbf{r}_f)=\frac{1}{8\pi(1-\nu)}\int_0^1 \mathbf{u}^*(\alpha) d\alpha
//          *	\f]
//          */
//			VectorDim u_source = VectorDim::Zero();
//			Quadrature<1,qOrder,QuadratureRule>::integrate(this,u_source,&LinkType::displacement_integrand,Rfield,S);
//			return Material<Isotropic>::C4*u_source;
//		}

//		/**********************************************************************/
//		VectorDim displacement_integrand(const int & k, const VectorDim & Rfield, const VectorDim& S) const __attribute__ ((deprecated))
//        {/*!@param[in] k			the current quadrature point
//          * @param[in] Rfield	the field point
//          * @param[in] S the unit vector necessary to compute displacement as line integral
//          * \returns The infinitesimal dispacement field generated by this segment at a field point
//          *
//          * The return value is calculated according to:
//          *	\f[
//          *		\mathbf{u}^*(\alpha_k,\mathbf{r}_f) = \frac{1}{R}\left\{ \frac{2(1-\nu)}{R+\mathbf{R}\cdot\mathbf{S}}\mathbf{b} \left[\left(\mathbf{S}\times\mathbf{R}\right)\cdot\frac{d\mathbf{r}_s}{d\alpha}\right]
//          *                                          + (1-2\nu)\left( \mathbf{b}\times \frac{d\mathbf{r}_s}{d\alpha}\right)
//          *                                          +\frac{1}{R^2}\left[\left(\frac{d\mathbf{r}_s}{d\alpha}\times\mathbf{b}\right)\cdot\mathbf{R}\right]\mathbf{R} \right\}
//          *	\f]
//          *	where:
//          *  - \f$\mathbf{r}(\alpha)\f$ is the parametrized source line segment;
//          *  - \f$\mathbf{R}(\mathbf{r}_f,\alpha) = \mathbf{r}_f-\mathbf{r}_s(\alpha)\f$ is vector connecting the source point to the field point.
//          *
//          *  The calculation of the solid angle is performed transforming the surface integral into a line integral. From Stokes theorem:
//          *	\f[
//          *  \oint\frac{\hat{\mathbf{S}}\times\mathbf{R}}{R(R-\mathbf{S}\cdot\mathbf{R})}\cdot d\mathbf{l}
//          *  \underbrace{=}_{\mbox{Stokes Th.}}\int\left(\nabla\times\frac{\hat{\mathbf{S}}\times\mathbf{R}}{R(R-\mathbf{S}\cdot\mathbf{R})}\right)\cdot\hat{\mathbf{n}}dA
//          *  \underbrace{=}_{\mbox{identity}}
//          * -\int\frac{\mathbf{R}}{R^3}\cdot\mathbf{n}dA
//          *	\f]
//          * References:
//          * [1] Asvestas, J. Line integrals and physical optics. Part I. The transformation of the solid-angle surface integral to a line integral. J. Opt. Soc. Am. A, 2(6), 891–895.
//          */
//
//			const VectorDim R(Rfield-rgauss.col(k));
////            const double RaSquared(R.squaredNorm()+coreLsquared);
//            const double RaSquared(R.squaredNorm()+DislocationStress<dim>::a2);
//			const double Ra=sqrt(RaSquared);
//			return 1.0/Ra * (+ 2.0*Material<Isotropic>::C1/(Ra+R.dot(S))*Burgers*(S.cross(R)).dot(rugauss.col(k))
//                             /*            */ + Material<Isotropic>::C3*rugauss.col(k).cross(Burgers)
//                             /*            */ + 1.0/RaSquared*(rugauss.col(k).cross(Burgers)).dot(R)*R
//                             /*            */ );
//		}
//
//		/**********************************************************************/
//		MatrixDim lattice_rotation_source(const VectorDim & Rfield) const __attribute__ ((deprecated))
//        {/*! The lattice rotation tensor generated by this dislocatio segment at
//          * @param[in] Rfield the vector representing the filed point
//          */
//			MatrixDim dg(elasticDistortion(Rfield));
//			return 0.5*(dg.transpose()-dg);
//		}

//		/**********************************************************************/
//		MatrixDim elasticDistortion(const VectorDim & Rfield) const __attribute__ ((deprecated))
//        {/*! The elasticDistortion tensor generated by this dislocatio segment at
//          * @param[in] Rfield the vector representing the filed point
//          */
//			MatrixDim temp(MatrixDim::Zero());
//			Quadrature<1,qOrder,QuadratureRule>::integrate(this,temp,&LinkType::elasticDistortion_integrand,Rfield);
//			return (Material<Isotropic>::C4 * temp);
//		}
//
//		/**********************************************************************/
//		MatrixDim elasticDistortion_integrand(const int & k, const VectorDim& Rfield) const __attribute__ ((deprecated))
//        {
//			const VectorDim R(Rfield-rgauss.col(k));
//			const double RaSquared ( R.squaredNorm() + DislocationStress<dim>::a2);
//			return  ( Material<Isotropic>::C1*( 2.0 + (3.0*DislocationStress<dim>::a2/RaSquared) ) * (Burgers*(rugauss.col(k).cross(R)).transpose() + (Burgers.cross(rugauss.col(k)))*R.transpose())
//					 + (rugauss.col(k).cross(Burgers))*R.transpose() + R * (rugauss.col(k).cross(Burgers)).transpose()
//					 + R.cross(Burgers).dot(rugauss.col(k)) * (3.0/RaSquared*R*R.transpose() - I) ) / std::pow(sqrt(RaSquared),3);
//		}


//        /* stress_source ******************************************************/
//		MatrixDim stress_source(const VectorDim & Rfield) const __attribute__ ((deprecated))
//        {/*! @param[in] Rfield the vector representing the filed point
//          *  The stress tensor generated by this dislocatio segment at Rfield
//          */
//			MatrixDim temp(MatrixDim::Zero());
//			Quadrature<1,qOrder,QuadratureRule>::integrate(this,temp,&LinkType::stress_integrand,Rfield);
//			return Material<Isotropic>::C2*(temp+temp.transpose());		// extract SYMMETRIC part of stress
//		}
//
//		/* stress_integrand ***************************************************/
//		MatrixDim stress_integrand(const int & k, const VectorDim& Rfield) const __attribute__ ((deprecated))
//        {/*! Returns the asymmetric (and dimensionless) part of the stress integrand generated by the current quadrature point.
//          * @param[in] k			the current quadrature point
//          * @param[in] Rfield	the vector connecting source point (corresponding to the current quadrature point) to field point
//          *
//          * The return value is calculated according to:
//          * Cai, W., Arsenlis, A., Weinberger, C., & Bulatov, V. (2006). A non-singular continuum theory of dislocations. Journal Of The Mechanics And Physics Of Solids, 54(3), 561–587.
//          *	\f[
//          *		d\mathbf{s} = (1-\nu) \left(1+\frac{3}{2}\frac{a^2}{R_a^2}\right)\frac{\partial \mathbf{r}}{\partial u}\otimes \left(\mathbf{b}\times \mathbf{R}\right)+
//          *		\mathbf{R}\otimes\left(\frac{\partial \mathbf{r}}{\partial u}\times\mathbf{b}\right)+
//          *		\frac{1}{2} \left[ \left(\mathbf{R}\times\mathbf{b}\right)\cdot \frac{\partial \mathbf{r}}{\partial u} \right]\left[\mathbf{I}\left(1+\frac{3a^2}{R_a^2}\right)+\frac{3}{R_a^2}\mathbf{R}\otimes\mathbf{R}\right]
//          *	\f]
//          *  where \f$R_a^2=|\mathbf{R}|^2+a^2\f$ is the modified squared norm of \f$\mathbf{R}\f$.
//          *
//          *	The return value is asymmetric and dimensionless in the sense that the actual stress integrand is:
//          *	\f[
//          *		d\mathbf{\sigma}=\frac{\mu}{4\pi (1-\nu)}\left(d\mathbf{s}+d\mathbf{s}^T\right)
//          *	\f]
//          */
//			VectorDim R(Rfield-rgauss.col(k));
//			double RaSquared(R.squaredNorm() + coreLsquared);
//			return   (Material<Isotropic>::C1*(1.0+1.5*coreLsquared/RaSquared)*rugauss.col(k)*(Burgers.cross(R)).transpose()
//					  + 	R*(rugauss.col(k).cross(Burgers)).transpose()
//					  +  0.5* R.cross(Burgers).dot(rugauss.col(k)) * (I*(1.0+3.0*coreLsquared/RaSquared) + 3.0/RaSquared*R*R.transpose())
//					  )/std::pow(sqrt(RaSquared),3);
//		}

