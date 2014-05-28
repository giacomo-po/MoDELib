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
#include <boost/ptr_container/ptr_map.hpp>
#include <boost/ptr_container/ptr_vector.hpp>


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

#include <model/DislocationDynamics/CrossSlip/CrossSlipSegment.h>
#include <model/DislocationDynamics/DislocationMobility.h>


#include <model/DislocationDynamics/VirtualBoundarySlipContainer.h>


namespace model {
    
	template <short unsigned int _dim, short unsigned int corder, typename InterpolationType,
	/*	   */ short unsigned int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule>
	class DislocationSegment :
	/*	                      */ public SplineSegmentBase<DislocationSegment<_dim,corder,InterpolationType,qOrder,QuadratureRule>,
	/*                                               */ _dim, corder>,
	/*	                      */ public GlidePlaneObserver<DislocationSegment<_dim,corder,InterpolationType,qOrder,QuadratureRule> >
    {
		
        
    public:
        
        enum{dim=_dim}; // make dim available outside class
        
        
		typedef DislocationSegment<dim,corder,InterpolationType,qOrder,QuadratureRule> Derived; 		// Define "Derived" so that NetworkTypedefs.h can be used
#include <model/Network/NetworkTypedefs.h>
#include <model/Geometry/Splines/SplineEnums.h>
        
        
		typedef SplineSegmentBase<Derived,dim,corder> SegmentBaseType;
		typedef std::map<size_t,LinkType* const> AddressMapType;
		typedef typename AddressMapType::iterator AddressMapIteratorType;
		typedef Eigen::Matrix<double,dim,qOrder>	MatrixDimQorder;
		typedef Eigen::Matrix<double, 1, qOrder>	VectorQorder;
		typedef QuadPow<Ncoeff-1,qOrder,QuadratureRule> QuadPowType;
        
        typedef typename GlidePlaneObserver<LinkType>::GlidePlaneType GlidePlaneType;
		typedef typename GlidePlaneObserver<LinkType>::GlidePlaneSharedPtrType GlidePlaneSharedPtrType;
        
        //        typedef DislocationQuadratureParticle<dim> DislocationParticleType;
        
        typedef std::vector<Eigen::Matrix<double,dim,1>> vector_VectorDim;
        
        typedef DislocationParticle<dim> DislocationParticleType;
        typedef std::vector<DislocationParticleType*> QuadratureParticleContainerType;
        
        
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
		
		
		/******************************************************************/
	public: //  data members
		/******************************************************************/
		
        //	//	boost::ptr_vector<DislocationParticleType> quadratureParticleContainer;
		
        //        std::vector<std::pair<size_t,DislocationParticle<dim>* const> > quadratureParticleContainer;
        QuadratureParticleContainerType quadratureParticleContainer;
        
        
		//! The Burgers vector
		const VectorDim Burgers;
        
        //! The glide plane unit normal vector
        const VectorDim   glidePlaneNormal;
        
		const VectorDim sessilePlaneNormal;
        
		static double coreLsquared;
		
        DislocationSharedObjects<Derived> shared;
        
		
		//! A shared pointer to the GlidePlane of this segment
		const GlidePlaneSharedPtrType pGlidePlane;
        
        
//        const Eigen::Matrix<double,1,2> dH0;
        
        const DislocationMobility<dim> dm;

        
//        const double& DHs;
//        const double& DHe;
        
        //! Positions corrersponding to the quadrature points
		MatrixDimQorder rgauss;
        
        //! PK force corrersponding to the quadrature points
        MatrixDimQorder pkGauss;
        
		/******************************************************************/
	private: //  member functions
		/******************************************************************/
		
        
        
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
            //			return temp.transpose()*pkGauss.col(k)*jgauss(k);
////			return temp.transpose()*radiativeVel(pkGauss.col(k))*jgauss(k); // inverse mobility law

			return temp.transpose()*radiativeVel(pkGauss.col(k))*jgauss(k); // inverse mobility law
//            return temp.transpose()*dm.getVelocity(pkGauss.col(k),rlgauss.col(k))*jgauss(k); // inverse mobility law

            
		}
		
		/**********************************************************************/
		VectorDim radiativeVel(const VectorDim& pkF) const
        {
			const VectorDim v0(Material<Isotropic>::Binv*pkF);
			const double v0N(v0.norm());
			const double csf(0.7*Material<Isotropic>::cs);
			return (v0N>FLT_EPSILON)? csf*(1.0-std::exp(-v0N/csf))*v0/v0N : v0;
		}
		
		
#ifdef UserStressFile
#include UserStressFile
#endif
        
		/******************************************************************/
	public: // member functions
		/******************************************************************/
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
		
		/* Constructor with Nodes and FLow ************************************/
		DislocationSegment(const std::pair<NodeType*,NodeType*> nodePair, const VectorDim & Fin) :
		/* base class initialization */ SegmentBaseType::SplineSegmentBase(nodePair,Fin) ,
		/* init list       */ Burgers(this->flow * Material<Isotropic>::b),
        /* init list       */ glidePlaneNormal(CrystalOrientation<dim>::find_planeNormal(nodePair.second->get_P()-nodePair.first->get_P(),Burgers).normalized()),
        /* init list       */ sessilePlaneNormal(CrystalOrientation<dim>::get_sessileNormal(nodePair.second->get_P()-nodePair.first->get_P(),Burgers)),
		/* init list       */ pGlidePlane(this->findExistingGlidePlane(glidePlaneNormal,this->source->get_P().dot(glidePlaneNormal))), // change this
        /* init list       */ dm(glidePlaneNormal,Burgers)
        {/*! Constructor with pointers to source and sink, and flow
          *  @param[in] NodePair_in the pair of source and sink pointers
          *  @param[in] Flow_in the input flow
          */
            
            this->source->make_planeNormals();
            DislocationEnergyRules<dim>::template findEdgeConfiguration<NodeType>(*this->source); // This should not be called in edge expansion or contraction
            this->source->make_T();
            
            this->sink->make_planeNormals();
            DislocationEnergyRules<dim>::template findEdgeConfiguration<NodeType>(*this->sink); // This should not be called in edge expansion or contraction
            this->sink->make_T();
            
			assert(this->flow.squaredNorm()>0.0);
			pGlidePlane->addToGLidePlane(this);
            pkGauss.setZero(); // necessary if this is not assembled
            
            
		}
		
		/* Constructor from EdgeExpansion) ************************************/
		DislocationSegment(const std::pair<NodeType*,NodeType*> nodePair, const ExpandingEdge<LinkType>& ee) :
		/* base class initialization */ SegmentBaseType::SplineSegmentBase(nodePair,ee),
		/* init list       */ Burgers(this->flow * Material<Isotropic>::b),
        /* init list       */ glidePlaneNormal(CrystalOrientation<dim>::find_planeNormal(nodePair.second->get_P()-nodePair.first->get_P(),Burgers).normalized()),
        /* init list       */ sessilePlaneNormal(CrystalOrientation<dim>::get_sessileNormal(nodePair.second->get_P()-nodePair.first->get_P(),Burgers)),
		/* init list       */ pGlidePlane(this->findExistingGlidePlane(glidePlaneNormal,this->source->get_P().dot(glidePlaneNormal))), 			// change this
        /* init list       */ dm(glidePlaneNormal,Burgers)
        {/*! Constructor with pointers to source and sink, and ExpandingEdge
          *  @param[in] NodePair_in the pair of source and sink pointers
          *  @param[in] ee the expanding edge
          */
            
            this->source->make_planeNormals();
            DislocationEnergyRules<dim>::template findEdgeConfiguration<NodeType>(*this->source); // This should not be called in edge expansion or contraction
            this->source->make_T();
            
            this->sink->make_planeNormals();
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
			if (shared.boundary_type==softBoundary && shared.use_bvp)
            {
				if(is_boundarySegment() && this->chord().norm()>FLT_EPSILON)
                {
                    shared.vbsc.add(*this);
				}
			}
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
            
            for (unsigned int k=0;k<qOrder;++k)
            {
                MatrixNcoeff  SFCH(this->get_SFCH());
                MatrixNcoeffDim qH(this->get_qH());
                //			QuadPowType::uPow.row(k)*SFCH;
                SFgauss.row(k)=QuadPowType::uPow.row(k)*SFCH; // WHY ARE WE LOOPING TO DO THIS MATRIX MULTIPLICATION????
                rgauss.col(k)=SFgauss.row(k)*qH; // WHY ARE WE LOOPING TO DO THIS MATRIX MULTIPLICATION????
                rugauss.col(k)=QuadPowType::duPow.row(k)*SFCH.template block<Ncoeff-1,Ncoeff>(1,0)*qH; // WHY ARE WE LOOPING TO DO THIS MATRIX MULTIPLICATION????
                jgauss(k)=rugauss.col(k).norm();
                rlgauss.col(k)=rugauss.col(k)/jgauss(k);
                
                quadratureParticleContainer.push_back(particleSystem.addParticle(rgauss.col(k),rugauss.col(k),Burgers,
                                                                                 Quadrature<1,qOrder,QuadratureRule>::abscissas(k),
                                                                                 Quadrature<1,qOrder,QuadratureRule>::weights(k))
                                                      );
            
            }
            
		}
		
		/**********************************************************************/
		VectorDim pkForce(const size_t & k)
        {/*!@param[in] k the k-th quandrature point
          *\returns the PK force at the k-th quandrature point
          */
            return (shared.use_bvp) ? ((quadratureParticleContainer[k]->stress(this->source->bvpStress,this->sink->bvpStress)+shared.vbsc.stress(quadratureParticleContainer[k]->P)+shared.externalStress)*Burgers).cross(rlgauss.col(k))
			/*                   */ : ((quadratureParticleContainer[k]->stress()+shared.externalStress)*Burgers).cross(rlgauss.col(k));
            
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
            
            if (this->pSN()->nodeOrder()>=shared.minSNorderForSolve) // check that SubNetwork has at least shared.minSNorderForSolve nodes
            {                
                /*! 1- Assembles the force vector of this segment.
                 *	\f[
                 *		\mathbf{K} = int_0^1 \mathbf{K}^*(u) du
                 *	\f]
                 */
                //                pkGauss.setZero(); // not strictly necessary
                for (int k=0;k<qOrder;++k){
                    pkGauss.col(k)=pkForce(k);
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
                
                
                for (int k=0;k<Mso.cols();++k){
                    const std::set<size_t>::const_iterator f(segmentDOFs.find(sourceDOFs(k)));
                    assert(f!=segmentDOFs.end());
                    unsigned int curCol(std::distance(segmentDOFs.begin(),f));
                    Mseg.template block<Ndof/2,1>(0,curCol)=Mso.col(k);
                }
                
                
                for (int k=0;k<Msi.cols();++k){
                    const std::set<size_t>::const_iterator f(segmentDOFs.find(sinkDOFs(k)));
                    assert(f!=segmentDOFs.end());
                    unsigned int curCol(std::distance(segmentDOFs.begin(),f));
                    Mseg.template block<Ndof/2,1>(Ndof/2,curCol)=Msi.col(k);
                }
                
            }
            
            
            
		}
        
        
        //        Eigen::Matrix<double,1,Ndof> ortC_integrand(const int& k) const {
        //            return rugauss.col(k).transpose()*SFgaussEx(k)*Quadrature<1,qOrder,QuadratureRule>::abscissa(k)*(1.0-Quadrature<1,qOrder,QuadratureRule>::abscissa(k));
        //            //                                              *Quadrature<1,qOrder,QuadratureRule>::abscissa(k)*(1.0-Quadrature<1,qOrder,QuadratureRule>::abscissa(k));
        //            //return rugauss.col(k).transpose()*SFgaussEx(k);
        //        }
        
        
        
        /* addToGlobalAssembly ************************************************/
        void addToGlobalAssembly(std::vector<Eigen::Triplet<double> >& kqqT,  Eigen::VectorXd& FQ) const
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
		const VectorNdof & get_Fq() const
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
			return (this->source->meshLocation() == onMeshBoundary && this->sink->meshLocation() == onMeshBoundary);
		}
		
		
		/**********************************************************************/
		VectorDim displacement(const VectorDim & Rfield, const VectorDim& S) const
        {/*!@param[in] Rfield the field point
          * @param[in] S the unit vector necessary to compute displacement as line integral
          * \returns the infinitesimal dispacement field generated by this segment at Rfield
          *
          * The return value is calculated according to:
          *	\f[
          *		\mathbf{u}(\mathbf{r}_f)=\frac{1}{8\pi(1-\nu)}\int_0^1 \mathbf{u}^*(\alpha) d\alpha
          *	\f]
          */
			VectorDim u_source = VectorDim::Zero();
			Quadrature<1,qOrder,QuadratureRule>::integrate(this,u_source,&LinkType::displacement_integrand,Rfield,S);
			return Material<Isotropic>::C4*u_source;
		}
		
		/**********************************************************************/
		VectorDim displacement_integrand(const int & k, const VectorDim & Rfield, const VectorDim& S) const
        {/*!@param[in] k			the current quadrature point
          * @param[in] Rfield	the field point
          * @param[in] S the unit vector necessary to compute displacement as line integral
          * \returns The infinitesimal dispacement field generated by this segment at a field point
          *
          * The return value is calculated according to:
          *	\f[
          *		\mathbf{u}^*(\alpha_k,\mathbf{r}_f) = \frac{1}{R}\left\{ \frac{2(1-\nu)}{R+\mathbf{R}\cdot\mathbf{S}}\mathbf{b} \left[\left(\mathbf{S}\times\mathbf{R}\right)\cdot\frac{d\mathbf{r}_s}{d\alpha}\right]
          *                                          + (1-2\nu)\left( \mathbf{b}\times \frac{d\mathbf{r}_s}{d\alpha}\right)
          *                                          +\frac{1}{R^2}\left[\left(\frac{d\mathbf{r}_s}{d\alpha}\times\mathbf{b}\right)\cdot\mathbf{R}\right]\mathbf{R} \right\}
          *	\f]
          *	where:
          *  - \f$\mathbf{r}(\alpha)\f$ is the parametrized source line segment;
          *  - \f$\mathbf{R}(\mathbf{r}_f,\alpha) = \mathbf{r}_f-\mathbf{r}_s(\alpha)\f$ is vector connecting the source point to the field point.
          *
          *  The calculation of the solid angle is performed transforming the surface integral into a line integral. From Stokes theorem:
          *	\f[
          *  \oint\frac{\hat{\mathbf{S}}\times\mathbf{R}}{R(R-\mathbf{S}\cdot\mathbf{R})}\cdot d\mathbf{l}
          *  \underbrace{=}_{\mbox{Stokes Th.}}\int\left(\nabla\times\frac{\hat{\mathbf{S}}\times\mathbf{R}}{R(R-\mathbf{S}\cdot\mathbf{R})}\right)\cdot\hat{\mathbf{n}}dA
          *  \underbrace{=}_{\mbox{identity}}
          * -\int\frac{\mathbf{R}}{R^3}\cdot\mathbf{n}dA
          *	\f]
          * References:
          * [1] Asvestas, J. Line integrals and physical optics. Part I. The transformation of the solid-angle surface integral to a line integral. J. Opt. Soc. Am. A, 2(6), 891–895.
          */
			
			const VectorDim R(Rfield-rgauss.col(k));
            const double RaSquared(R.squaredNorm()+coreLsquared);
			const double Ra=sqrt(RaSquared);
			return 1.0/Ra * (+ 2.0*Material<Isotropic>::C1/(Ra+R.dot(S))*Burgers*(S.cross(R)).dot(rugauss.col(k))
                             /*            */ + Material<Isotropic>::C3*rugauss.col(k).cross(Burgers)
                             /*            */ + 1.0/RaSquared*(rugauss.col(k).cross(Burgers)).dot(R)*R
                             /*            */ );
		}
		
		/**********************************************************************/
		MatrixDim lattice_rotation_source(const VectorDim & Rfield) const
        {/*! The lattice rotation tensor generated by this dislocatio segment at
          * @param[in] Rfield the vector representing the filed point
          */
			MatrixDim dg(elasticDistortion(Rfield));
			return 0.5*(dg.transpose()-dg);
		}
		
		/**********************************************************************/
		MatrixDim elasticDistortion(const VectorDim & Rfield) const
        {/*! The elasticDistortion tensor generated by this dislocatio segment at
          * @param[in] Rfield the vector representing the filed point
          */
			MatrixDim temp(MatrixDim::Zero());
			Quadrature<1,qOrder,QuadratureRule>::integrate(this,temp,&LinkType::elasticDistortion_integrand,Rfield);
			return (Material<Isotropic>::C4 * temp);
		}
		
		/**********************************************************************/
		MatrixDim elasticDistortion_integrand(const int & k, const VectorDim& Rfield) const
        {
			const VectorDim R(Rfield-rgauss.col(k));
			const double RaSquared ( R.squaredNorm() + coreLsquared);
			return  ( Material<Isotropic>::C1*( 2.0 + (3.0*coreLsquared/RaSquared) ) * (Burgers*(rugauss.col(k).cross(R)).transpose() + (Burgers.cross(rugauss.col(k)))*R.transpose())
					 + (rugauss.col(k).cross(Burgers))*R.transpose() + R * (rugauss.col(k).cross(Burgers)).transpose()
					 + R.cross(Burgers).dot(rugauss.col(k)) * (3.0/RaSquared*R*R.transpose() - I) ) / std::pow(sqrt(RaSquared),3);
		}
		
		/**********************************************************************/
		Eigen::Matrix<double,dim-1,Ncoeff> hermiteLocalCoefficient() const
        {
			const MatrixDim G2L(DislocationLocalReference<dim>::global2local(this->chord(),glidePlaneNormal));
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
		std::map<double,VectorDim> boundaryCollision() const
        {/*! The collision points with the glidePlane boundary polygon.
          */
			std::map<double,VectorDim> temp;
			//	std::set<double> temp;
			
			
			if (shared.boundary_type){
				VectorDim origin = this->source->get_P();
				
				const Eigen::Matrix<double,dim-1,Ncoeff> C1L(polynomialLocalCoeff()); // the local polynomial coefficients of this
                PlanarSplineImplicitization<pOrder> psi(Coeff2Hermite<pOrder>::template h2c<dim-1>(C1L));
                
				const MatrixDim G2L(DislocationLocalReference<dim>::global2local(this->chord(),glidePlaneNormal));
				
				for (int k=0; k<pGlidePlane->segmentMeshCollisionPairContainer.size();++k)
                {
					Eigen::Matrix<double,dim-1,2> H2L;
					H2L.col(0)=(G2L*(pGlidePlane->segmentMeshCollisionPairContainer[k].first -origin)).template segment<dim-1>(0);
					H2L.col(1)=(G2L*(pGlidePlane->segmentMeshCollisionPairContainer[k].second-origin)).template segment<dim-1>(0);
                    //	PlanarSplineImplicitization<pOrder> psi(Coeff2Hermite<pOrder>::template h2c<dim-1>(C1L));
					std::set<std::pair<double,double> > intersectinParameters = psi.template intersectWith<1>(Coeff2Hermite<1>::template h2c<dim-1>(H2L),FLT_EPSILON);
					for (std::set<std::pair<double,double> >::const_iterator iter=intersectinParameters.begin();iter!=intersectinParameters.end();++iter)
                    {
						temp.insert(std::make_pair(iter->first,this->get_r(iter->first)));
					}
				}
			}
			
			return temp;
		}
		
		/*********************************************************************/
		CrossSlipSegment<LinkType> isCrossSlipSegment(const double& sinThetaCrossSlipCr,const double& crossSlipLength) const
        {
			return CrossSlipSegment<LinkType>(*this,sinThetaCrossSlipCr,crossSlipLength);
		}
        
		/*********************************************************************/
        vector_VectorDim conjugatePlaneNormal() const
        {
			return CrystalOrientation<dim>::conjugatePlaneNormal(Burgers,glidePlaneNormal);
		}
		
		/*********************************************************************/
		VectorDim integralPK() const
        {
			return pkGauss.col(qOrder/2)*this->chord().norm();
		}
        
        /* stress_source ******************************************************/
		MatrixDim stress_source(const VectorDim & Rfield) const __attribute__ ((deprecated))
        {/*! @param[in] Rfield the vector representing the filed point
          *  The stress tensor generated by this dislocatio segment at Rfield
          */
			MatrixDim temp(MatrixDim::Zero());
			Quadrature<1,qOrder,QuadratureRule>::integrate(this,temp,&LinkType::stress_integrand,Rfield);
			return Material<Isotropic>::C2*(temp+temp.transpose());		// extract SYMMETRIC part of stress
		}
        
		/* stress_integrand ***************************************************/
		MatrixDim stress_integrand(const int & k, const VectorDim& Rfield) const __attribute__ ((deprecated))
        {/*! Returns the asymmetric (and dimensionless) part of the stress integrand generated by the current quadrature point.
          * @param[in] k			the current quadrature point
          * @param[in] Rfield	the vector connecting source point (corresponding to the current quadrature point) to field point
          *
          * The return value is calculated according to:
          * Cai, W., Arsenlis, A., Weinberger, C., & Bulatov, V. (2006). A non-singular continuum theory of dislocations. Journal Of The Mechanics And Physics Of Solids, 54(3), 561–587.
          *	\f[
          *		d\mathbf{s} = (1-\nu) \left(1+\frac{3}{2}\frac{a^2}{R_a^2}\right)\frac{\partial \mathbf{r}}{\partial u}\otimes \left(\mathbf{b}\times \mathbf{R}\right)+
          *		\mathbf{R}\otimes\left(\frac{\partial \mathbf{r}}{\partial u}\times\mathbf{b}\right)+
          *		\frac{1}{2} \left[ \left(\mathbf{R}\times\mathbf{b}\right)\cdot \frac{\partial \mathbf{r}}{\partial u} \right]\left[\mathbf{I}\left(1+\frac{3a^2}{R_a^2}\right)+\frac{3}{R_a^2}\mathbf{R}\otimes\mathbf{R}\right]
          *	\f]
          *  where \f$R_a^2=|\mathbf{R}|^2+a^2\f$ is the modified squared norm of \f$\mathbf{R}\f$.
          *
          *	The return value is asymmetric and dimensionless in the sense that the actual stress integrand is:
          *	\f[
          *		d\mathbf{\sigma}=\frac{\mu}{4\pi (1-\nu)}\left(d\mathbf{s}+d\mathbf{s}^T\right)
          *	\f]
          */
			VectorDim R(Rfield-rgauss.col(k));
			double RaSquared(R.squaredNorm() + coreLsquared);
			return   (Material<Isotropic>::C1*(1.0+1.5*coreLsquared/RaSquared)*rugauss.col(k)*(Burgers.cross(R)).transpose()
					  + 	R*(rugauss.col(k).cross(Burgers)).transpose()
					  +  0.5* R.cross(Burgers).dot(rugauss.col(k)) * (I*(1.0+3.0*coreLsquared/RaSquared) + 3.0/RaSquared*R*R.transpose())
					  )/std::pow(sqrt(RaSquared),3);
		}

		
		
		/**********************************************************************/
		template <class T>
		friend T& operator << (T& os, const Derived& ds)
        {
			os  << ds.source->sID<<"\t"<< ds.sink->sID<<"\t"
			/**/<< std::setprecision(15)<<std::scientific<<ds.flow.transpose()<<"\t"
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
    
    template <short unsigned int dim, short unsigned int corder, typename InterpolationType,
    /*	   */ short unsigned int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule>
    double DislocationSegment<dim,corder,InterpolationType,qOrder,QuadratureRule>::coreLsquared=1.0;
    
    
    //////////////////////////////////////////////////////////////s
} // namespace model
#endif

            //		/* intersectWith ******************************************************/
            //		std::set<std::pair<double,double> > intersectWith( const Derived* const p_other) const {
            ////			return DislocationSegmentIntersection<dim,pOrder>(this->hermiteCoefficients(),glidePlaneNormal).intersectWith(p_other->hermiteCoefficients(),p_other->glidePlaneNormal,tol,p_other);
            //			return DislocationSegmentIntersection<dim,pOrder>(this->hermiteCoefficients(),glidePlaneNormal).intersectWith(p_other->hermiteCoefficients(),p_other->glidePlaneNormal);
            //		}

            
            
            //		/**********************************************************************/
            //		double energy() const
            //        {/*! The total elastic energy generated by this segment.
            //          *  Includes self-energy and interation energy between this segment and other segments in the network.
            //          */
            //			double temp=0.0;
            //			Quadrature<1,qOrder,QuadratureRule>::integrate(this,temp,&LinkType::energy_field);
            //			return temp;
            //		}
            //
            //
            //		//////////////////////////////////
            //		double energy_field(const int & k) const
            //        {/*!@param[in] k the k-th quadrature point
            //          *
            //          * \returns The elastic energy at the k-th quadrature
            //          * point due to other dislocation segments.
            //          */
            //
            //			double temp=0.0;
            //			for (AddressMapIteratorType aIter=this->ABbegin(); aIter!=this->ABend();++aIter){
            //				temp+=aIter->second->energy_source(rgauss.col(k), rugauss.col(k), Burgers);
            //			}
            //			return temp;
            //		}
            //
            //		//////////////////////////////////
            //		double energy_source(const VectorDim & rf, const VectorDim & ruf, const VectorDim & bf) const
            //        {/*! The elastic energy generated by this (source) segment at the input field point
            //          */
            //			double temp(0.0);
            //			Quadrature<1,qOrder,QuadratureRule>::integrate(this,temp,&LinkType::energy_integrand,rf,ruf,bf);
            //			return -Material<Isotropic>::C2*temp;
            //		}
            //
            //
            //        /**********************************************************************/
            //        double energy_integrand(const int & k,  VectorDim & rf, const VectorDim & ruf, const VectorDim & bf) const
            //        {/*!
            //          * \returns the energy_integrand
            //          */
            //			VectorDim DR= rf-rgauss.col(k);
            //			double RaSquared = DR.squaredNorm() + coreLsquared;
            //			return  (Material<Isotropic>::C1*(1.0+0.5*coreLsquared/RaSquared)*Burgers.dot(rugauss.col(k))*bf.dot(ruf)
            //					 +2.0*Material<Isotropic>::nu*(1.0+0.5*coreLsquared/RaSquared)*(bf.dot(rugauss.col(k))*Burgers.dot(ruf))
            //					 -(Burgers.dot(bf)*(1.0+coreLsquared/RaSquared)+ Burgers.dot(DR)*bf.dot(DR)*1.0/RaSquared )*ruf.dot(rugauss.col(k))
            //					 )/sqrt(RaSquared);
            //		}

            
            
//
//        /* stress_field *******************************************************/
//		MatrixDim stress_field(const size_t & k)
//        {/*! @param[in] k			the k-th quadrature point
//          *  Calculates the total stress field at the k-th gauss (field) point
//          */
//
//			//! 1- Create and set to zero a dim x dim matrix
//			MatrixDim sigma( MatrixDim::Zero() );
//
//			//! 2- Loop over other segments summing the strees field at the k-th field point
//			for (AddressMapIteratorType aIter=this->ABbegin(); aIter!=this->ABend();++aIter){
//				sigma+=aIter->second->stress_source(rgauss.col(k));
//			}
//
//			// 3- If ExternalStressFile is defined add the external stress
//#ifdef UserStressFile
//			sigma+=userStress(k);
//#endif
//
//			if (shared.use_bvp){
//				//	std::cout<<"BVP stress"<<std::endl;
//				//	std::cout<<shared.domain.stressAt(rgauss.col(k))<<std::endl;
//				sigma+=this->source->bvpStress*(1.0-Quadrature<1,qOrder,QuadratureRule>::abscissa(k))+this->sink->bvpStress*Quadrature<1,qOrder,QuadratureRule>::abscissa(k);
//			}
//
//			//return sigma;
//			//			return sigma+shared.loadController.externalStress();
//			return sigma+shared.externalStress;
//		}


