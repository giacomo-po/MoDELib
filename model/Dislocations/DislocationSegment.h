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

#ifndef VERBOSELEVEL
#define VERBOSELEVEL 3
#endif


#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#endif

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <set>
#include <boost/ptr_container/ptr_map.hpp>
#include <boost/ptr_container/ptr_vector.hpp>


#include <model/Network/Operations/EdgeExpansion.h>


#include <model/Quadrature/Quadrature.h>
#include <model/Quadrature/QuadPow.h>
#include <model/Dislocations/DislocationNetworkTraits.h>
#include <model/Dislocations/DislocationConsts.h>
#include <model/Geometry/Splines/SplineSegmentBase.h>
#include <model/Dislocations/Materials/SlipSystem.h>
#include <model/Dislocations/DislocationSharedObjects.h>
#include <model/Dislocations/GlidePlanes/GlidePlaneObserver.h>
#include <model/Dislocations/GlidePlanes/GlidePlane.h>
#include <model/Geometry/Splines/Intersection/PlanarSplineImplicitization.h>
#include <model/Geometry/Splines/Coeff2Hermite.h>
#include <model/Utilities/CompareVectorsByComponent.h>

#include <model/Dislocations/DislocationQuadratureParticle.h>
#include <model/Dislocations/DislocationLocalReference.h>
#include <model/Dislocations/Junctions/DislocationSegmentIntersection.h>

#include <model/Dislocations/CrossSlip/CrossSlipSegment.h>



namespace model {
	
	
	//	double cellSize= 300.0;
	double cellSize= 1000.0;
	
	
	template <short unsigned int dim, typename MaterialType>
	struct DislocationSegmentInitializeBeforeBase{
		
		DislocationSharedObjects<dim,MaterialType> shared;
		
		
		typedef Eigen::Matrix<double,dim,1> VectorDim;

		/********************************************/
		VectorDim find_planeNormal(const VectorDim& chord, const VectorDim& Burgers/*,const VectorDim& T1,const VectorDim& T2*/){
			enum {Nslips=MaterialType::Nslips};
			std::set<SlipSystem<dim,Nslips> > allowedSlipSystems;
			shared.material.find_slipSystem(chord,Burgers,/*T1,T2,*/allowedSlipSystems);
			return allowedSlipSystems.begin()->normal; // DON't LIKE THIS
		}

		/********************************************/
		VectorDim get_sessileNormal(const VectorDim& chord, const VectorDim& Burgers){
			assert(chord.norm()>FLT_EPSILON && "CHORD TOO SMALL");
			assert(Burgers.norm()>FLT_EPSILON && "Burgers TOO SMALL");
			VectorDim temp(chord.normalized().cross(Burgers));
			double tempNorm(temp.norm());
			if (tempNorm>FLT_EPSILON){
				temp.normalize();
			}
			else{
				temp.setZero();
			}
			return temp;
		}
		
		
		//! The glide plane normal
		const VectorDim glidePlaneNormal;
		const VectorDim sessilePlaneNormal;
		
		
		
		// Constructor with chord and Burgers
		DislocationSegmentInitializeBeforeBase(const VectorDim& chord, const VectorDim& Burgers) : glidePlaneNormal(find_planeNormal(chord,Burgers).normalized()),
		/*                                                                                      */ sessilePlaneNormal(get_sessileNormal(chord,Burgers)){}
		
		// Constructor with plane normal
		DislocationSegmentInitializeBeforeBase(const VectorDim& normal_in, const VectorDim& chord, const VectorDim& Burgers) : glidePlaneNormal(normal_in.normalized()),
		/*                                                                                                                  */ sessilePlaneNormal(get_sessileNormal(chord,Burgers)){}
		
		
		
	};
	
	
	
	template <short unsigned int dim, short unsigned int corder, typename InterpolationType,
	/*	   */ double & alpha, short unsigned int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule, 
	/*	   */ typename MaterialType>
	class DislocationSegment : public DislocationSegmentInitializeBeforeBase<dim,MaterialType>, // This must be the first base class in the inheritance structure
	/*	                    */ public SplineSegmentBase<DislocationSegment<dim,corder,InterpolationType,alpha,qOrder,QuadratureRule,MaterialType>,
	/*                                               */ dim, corder, alpha>,
	/*	                    */ public GlidePlaneObserver<dim,DislocationSegment<dim,corder,InterpolationType,alpha,qOrder,QuadratureRule,MaterialType> >{
		

public:		
		typedef DislocationSegment<dim,corder,InterpolationType,alpha,qOrder,QuadratureRule,MaterialType> Derived; 		// Define "Derived" so that NetworkTypedefs.h can be used
#include <model/Network/NetworkTypedefs.h>		
#include <model/Geometry/Splines/SplineEnums.h>
		//		typedef SplineSegmentBase<Derived,dim,corder,alpha,qOrder,QuadratureRule> SegmentBaseType;
		typedef SplineSegmentBase<Derived,dim,corder,alpha> SegmentBaseType;
		
		typedef std::map<size_t,LinkType* const> AddressMapType;
		typedef typename AddressMapType::iterator AddressMapIteratorType;
		typedef Eigen::Matrix<double,dim,qOrder>	MatrixDimQorder;
		typedef Eigen::Matrix<double, 1, qOrder>	VectorQorder;
		typedef QuadPow<Ncoeff-1,qOrder,QuadratureRule> QuadPowType;
		
		private:		
		//		typedef Eigen::Matrix<double,dim,Ndof> MatrixDimNdof;
		
		
		//! Positions corrersponding to the quadrature points
		MatrixDimQorder rgauss;
		
		//! Parametric tangents at the quadrature points
		MatrixDimQorder rugauss;
		
		//! Scalar jacobian corrersponding to the quadrature points
		VectorQorder jgauss;
		
		//! Tangents corrersponding to the quadrature points
		MatrixDimQorder rlgauss;
		
		//protected:
		
		Eigen::Matrix<double,qOrder,Ncoeff> SFgauss;
		
		
		
		//		using MaterialType::Nslips;
		
		
		enum {Nslips=MaterialType::Nslips};
		
		
		DislocationSharedObjects<dim,MaterialType> shared;
		
		
		
		
	private:		
		
		//! The static MaterialType material
		static MaterialType material;
		
		
		//! See G. Schoeck "Atomic dislocation core parameters", 2009
		
		
		//VectorDim R;
		//double Rsquared;
		double C1, C2, C3, C4;			// THIS SHOULD BE A STATIC DATA MEMBER
		
		MatrixDimQorder pkGauss;
		//VectorDim S;
		
		
		
		//! Segment Stiffness Matrix
		MatrixNdof Kqq;
		
		
		//! Segment Nodal Force Vector
		//		typedef Eigen::Matrix<double,Ndof,1> VectorNdof;
		VectorNdof Fq;
		
		static const Eigen::Matrix<double,dim,dim> I;
		
		//////////////////////////////////
		//! assemble_Kqq
		void assemble_Kqq(){
			/*! Assembles the stiffness matrix of this segment.
			 *	\f[
			 *		\mathbf{K} = int_0^1 \mathbf{K}^*(u) du
			 *	\f]
			 */
			Kqq.setZero();
			Quadrature<1,qOrder,QuadratureRule>::integrate(this,Kqq,&LinkType::stiffness_integrand);
		}
		
		MatrixDimNdof SFgaussEx(const int& k) const { // GENERALIZE THIS
			return (MatrixDimNdof()<<I*SFgauss(k,0),I*SFgauss(k,1),I*SFgauss(k,2),I*SFgauss(k,3)).finished();
		}
		
		
		//	MatrixNcoeff stiffness_integrand(const int& k) const {
		MatrixNdof stiffness_integrand(const int& k) const {
			/*! The stiffness matrix integrand evaluated at the k-th quadrature point.
			 *  @param[in] k the current quadrature point
			 *	\f[
			 *		\mathbf{K}^* = \mathbf{N}^T \mathbf{B} \mathbf{N} \frac{dl}{du}
			 *	\f]
			 */
			//			return SFgauss.row(k).transpose()*shared.material.B*SFgauss.row(k)*jgauss(k);
			//			MatrixDimNdof temp((I+rlgauss.col(k)*rlgauss.col(k).transpose())*SFgaussEx(k));
			
			//			MatrixDimNdof temp((I-rlgauss.col(k)*rlgauss.col(k).transpose())*SFgaussEx(k));
			MatrixDimNdof temp(SFgaussEx(k));
			
			return temp.transpose()*shared.material.B*temp*jgauss(k);
		}
		
		//////////////////////////////////
		//! assemble_Fq
		void assemble_Fq(){
			/*! Assembles the force vector of this segment.
			 *	\f[
			 *		\mathbf{K} = int_0^1 \mathbf{K}^*(u) du
			 *	\f]
			 */
			
			// Compute and store PK force at quadrature points
			pkGauss.setZero();
//			int k;
//#ifdef _OPENMP
//#pragma omp parallel for default(shared) private(k)
//#endif
			for (int k=0;k<qOrder;++k){
				pkGauss.col(k)=pkForce(k);	
			}
			
			//				pkGauss.col(k)=quadratureParticleVector[k].pkForce();
			
			
			// Reset Fq and integrate
			Fq.setZero();
			Quadrature<1,qOrder,QuadratureRule>::integrate(this,Fq,&LinkType::force_integrand);
		}
		
		//////////////////////////////////
		// force_integrand
		VectorNdof force_integrand(const int& k) const {
			/*! The force vector integrand evaluated at the k-th quadrature point.
			 *  @param[in] k the current quadrature point
			 */
			//			VectorNdof temp;
			//			
			//			for (int n=0;n<Ncoeff;++n){
			//				temp.template segment<dim>(n*dim)= SFgauss.row(k)(n)*pkGauss.col(k);
			//			}
			//			
			//			return temp*jgauss(k);
			
			//			MatrixDimNdof temp((I+rlgauss.col(k)*rlgauss.col(k).transpose())*SFgaussEx(k));
			
			//			MatrixDimNdof temp((I-rlgauss.col(k)*rlgauss.col(k).transpose())*SFgaussEx(k));
			MatrixDimNdof temp(SFgaussEx(k));
			
			//			for (int n=0;n<Ncoeff;++n){
			//				temp.template segment<dim>(n*dim)= SFgauss.row(k)(n)*pkGauss.col(k);
			//			}
			
			return temp.transpose()*pkGauss.col(k)*jgauss(k);
			
		}
		
		
		
		
#ifdef UserStressFile
#include UserStressFile
#endif
		
		
		
		
		/******************************************************************/		
	public: //  data members
		/******************************************************************/
		
		//		typedef DislocationQuadratureParticle<dim,cellSize,qOrder,QuadratureRule> DislocationQuadratureParticleType;
		
		typedef DislocationQuadratureParticle<dim,cellSize> DislocationQuadratureParticleType;
		boost::ptr_vector<DislocationQuadratureParticleType> quadratureParticleVector;
		
		
		//! The Burgers vector
		const VectorDim Burgers;
		const double coreL;				// THIS SHOULD BE A STATIC DATA MEMBER
		const double coreLsquared;		// THIS SHOULD BE A STATIC DATA MEMBER
		//! the unit vector normal to the glide plane !!! FINISH HERE, USE INITIALIZATION LIST IN CONSTRUCTOR !!!
		//const VectorDim planeNormal;	
		//		const MatrixDim G2L;
		
		typedef MaterialType TempMaterialType;
		
		typedef typename GlidePlaneObserver<dim,LinkType>::GlidePlaneType GlidePlaneType;
		typedef typename GlidePlaneObserver<dim,LinkType>::GlidePlaneSharedPtrType GlidePlaneSharedPtrType;
		const GlidePlaneSharedPtrType pGlidePlane;
		
		
		/******************************************************************/		
	public: // member functions
		/******************************************************************/
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
		
		/* Constructor with Nodes and FLow ****************************/
		//! Constructor with Nodes and FLow		
		DislocationSegment(const std::pair<NodeType*,NodeType*> nodePair, const VectorDim & Fin) : 
		/* base class initialization */ DislocationSegmentInitializeBeforeBase<dim,MaterialType>::DislocationSegmentInitializeBeforeBase(nodePair.second->get_P()-nodePair.first->get_P(),Fin),
		/* base class initialization */ SegmentBaseType::SplineSegmentBase(nodePair,Fin) , 
		/* initialization list       */ Burgers(this->flow * shared.material.b),
		/* initialization list       */ coreL(1.0*shared.material.b),
		/* initialization list       */ coreLsquared(std::pow(coreL,2)),
		//		/* initialization list       */ G2L(DislocationLocalReference<dim>::global2local(Burgers,this->glidePlaneNormal)), //NO, B and N may not be orthogonal
		/* initialization list       */ pGlidePlane(this->findExistingGlidePlane(this->glidePlaneNormal,this->source->get_P().dot(this->glidePlaneNormal))){			// change this		
			
			
#if VERBOSELEVEL <= 3
			std::cout<<"Creating DislocationSegment with Burgers "<<Burgers.transpose()<<" on "<<this->glidePlaneNormal.transpose()<<" plane."<<std::endl;
#endif
			
			
			assert(this->flow.squaredNorm()>0.0);
			
			C1=1.0-shared.material.nu;
			C2=shared.material.mu/(4.0*M_PI*C1);
			C3=1.0-2.0*shared.material.nu;
			C4=1.0/(8.0*M_PI*C1);
			
			pGlidePlane->addToGLidePlane(this);	
		}
		
		
		
		
		
		
		/* Constructor with Nodes and Link (Expansion) ****************************/
		//! Constructor with Nodes and FLow		
		//		DislocationSegment(const std::pair<NodeType*,NodeType*> nodePair, const LinkType* const & pL) : 
		DislocationSegment(const std::pair<NodeType*,NodeType*> nodePair, const ExpandingEdge<LinkType>& ee) : 
		//      /* base class initialization */ DislocationSegmentInitializeBeforeBase<dim,MaterialType>::DislocationSegmentInitializeBeforeBase(ee.E.glidePlaneNormal,nodePair.second->get_P()-nodePair.first->get_P(),ee.E.flow),
		/* base class initialization */ DislocationSegmentInitializeBeforeBase<dim,MaterialType>::DislocationSegmentInitializeBeforeBase(nodePair.second->get_P()-nodePair.first->get_P(),ee.E.flow),
		/* base class initialization */ SegmentBaseType::SplineSegmentBase(nodePair,ee), 
		/* initialization list       */ Burgers(this->flow * shared.material.b),
		/* initialization list       */ coreL(1.0*shared.material.b),
		/* initialization list       */ coreLsquared(std::pow(coreL,2)),
		//		/* initialization list       */ G2L(DislocationLocalReference<dim>::global2local(Burgers,this->glidePlaneNormal)), //NO, B and N may not be orthogonal
		/* initialization list       */ pGlidePlane(this->findExistingGlidePlane(this->glidePlaneNormal,this->source->get_P().dot(this->glidePlaneNormal))){			// change this		
			
			
			
#if VERBOSELEVEL <= 3
			std::cout<<"Creating DislocationSegment from Expansion with Burgers "<<Burgers.transpose()<<" on "<<this->glidePlaneNormal.transpose()<<" plane."<<std::endl;
#endif
			
			assert(this->flow.squaredNorm()>0.0);
			
			C1=1.0-shared.material.nu;
			C2=shared.material.mu/(4.0*M_PI*C1);
			C3=1.0-2.0*shared.material.nu;
			C4=1.0/(8.0*M_PI*C1);
			
			pGlidePlane->addToGLidePlane(this);
		}
		
		/* Destructor *************************************************/
		~DislocationSegment(){
			pGlidePlane->removeFromGlidePlane(this);
		}
		
		/* updateQuadGeometryKernel ***********************************/
		void updateQuadGeometryKernel(const int& k){ 
			
			MatrixNcoeff  SFCH(this->get_SFCH());
			MatrixNcoeffDim qH(this->get_qH());
			//			QuadPowType::uPow.row(k)*SFCH;
			SFgauss.row(k)=QuadPowType::uPow.row(k)*SFCH;
			rgauss.col(k)=SFgauss.row(k)*qH;
			rugauss.col(k)=QuadPowType::duPow.row(k)*SFCH.template block<Ncoeff-1,Ncoeff>(1,0)*qH;
			jgauss(k)=rugauss.col(k).norm();
			rlgauss.col(k)=rugauss.col(k)/jgauss(k);
			
			quadratureParticleVector.push_back(new DislocationQuadratureParticleType(Quadrature<1,qOrder,QuadratureRule>::abscissas(k),
			/*                                                                    */ Quadrature<1,qOrder,QuadratureRule>::weights(k),
			/*                                                                    */ rgauss.col(k),rugauss.col(k),Burgers));
		}
		
		
		
		
		/********************************************************/
		MatrixDim stress_source(const VectorDim & Rfield) const {
			/*! The stress tensor generated by this dislocatio segment at   
			 * @param[in] Rfield the vector representing the filed point
			 */
			MatrixDim sigma_source = MatrixDim::Zero();
			Quadrature<1,qOrder,QuadratureRule>::integrate(this,sigma_source,&LinkType::stress_integrand,Rfield);
			
			
			/* IF THIS SEGMENT IS AN END SEGMENT ADD HERE THE STRESS DUE TO AN INFINITE STRAIGHT SEGMENT 
			 GOING OUT OF THE BOUNDARY AND WITH LINE DIRECTION EQUAL TO THE END POINT LINE DIRECTION
			 SEE CAI.
			 */
			
			//			if ((this->sink->constraintType==boundaryNode)*0){
			//				// in this case x2 is the point at infinity and x1=this->sink
			//				// also Rfield-x2 is anti-parallel to x2-x1
			//				// therefore T(Rfield-x2)-T(Rfield-x1) becomes:
			//				sigma_source+= stress_straight_inf(this->get_ru(1.0)) - stress_straight(Rfield-this->get_r(1.0),this->get_ru(1.0));
			//			}
			//			
			//			if ((this->source->constraintType==boundaryNode)*0){
			//				// in this case x1 is the point at infinity and x2=this->source
			//				// also Rfield-x1 is parallel to x2-x1
			//				// therefore T(Rfield-x2)-T(Rfield-x1) becomes:
			//				sigma_source+= stress_straight(Rfield-this->get_r(0.0),this->get_ru(0.0)) - stress_straight_inf(this->get_ru(0.0));
			//			}
			
			return C2*(sigma_source+sigma_source.transpose());		// extract SYMMETRIC part of stress
		}
		
		/********************************************************/
		MatrixDim stress_straight(const VectorDim & R, const VectorDim & t) const {
			
			const VectorDim T=t.normalized();
			const double RdotT=R.dot(t);
			const double RdotT2=std::pow(RdotT,2);
			const double RaSquared = R.squaredNorm()+coreLsquared;
			const double Ra = std::pow(RaSquared,0.5);
			const double RaCubed = std::pow(Ra,3);
			const double A1 = - RdotT*(3.0*RaSquared-RdotT2)/std::pow(RaSquared-RdotT2,2)/RaCubed;
			const double A2 = 1.0/RaCubed-RdotT*A1;
			const double A6 = - RdotT/(RaSquared-RdotT2)/Ra;
			const double A3 = - RdotT/Ra+A6+RdotT2*A1;
			const double A4 = A6 + coreLsquared*A1;
			const double A5 = -C1*A6-coreLsquared*C1*A1*0.5;
			const double A7 = shared.material.nu/Ra - RdotT*A6 - coreLsquared*C1*A2*0.5;
			
			return R.cross(Burgers).dot(t)*(0.5*A1*R*R.transpose() + A2*t*R.transpose() + 0.5*A3*t*t.transpose() + 0.5*A4*I)
			/*  */ + A5*R.cross(Burgers)*t.transpose() + A6*t.cross(Burgers)*R.transpose()
			/*  */ + A7*t.cross(Burgers)*t.transpose();		// extract SYMMETRIC part of stress
		}
		
		/********************************************************/
		MatrixDim stress_straight_inf(const VectorDim & t) const {
			const VectorDim T=t.normalized();
			
			return MatrixDim::Zero();		// extract SYMMETRIC part of stress
		}
		
				
		/********************************************************/
		MatrixDim stress_integrand(const int & k, const VectorDim& Rfield) const {
			/*! Returns the asymmetric (and dimensionless) part of the stress integrand generated by the current quadrature point.
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
			VectorDim R=Rfield-rgauss.col(k);
			double RaSquared = R.squaredNorm() + coreLsquared;
			return   (C1*(1.0+1.5*coreLsquared/RaSquared)*rugauss.col(k)*(Burgers.cross(R)).transpose()
					  + 	R*(rugauss.col(k).cross(Burgers)).transpose() 
					  +  0.5* R.cross(Burgers).dot(rugauss.col(k)) * (I*(1.0+3.0*coreLsquared/RaSquared) + 3.0/RaSquared*R*R.transpose())
					  )/std::pow(RaSquared,1.5);
		}
		
		//////////////////////////////////
		//! stress_field. Calculates the total stress field at the k-th gauss (field) point 
		MatrixDim stress_field(const size_t & k) {
			
			//! 1- Create and set to zero a dim x dim matrix
			MatrixDim sigma = MatrixDim::Zero();
			
			//! 2- Loop over other segments summing the strees field at the k-th field point
			for (AddressMapIteratorType aIter=this->ABbegin(); aIter!=this->ABend();++aIter){
				sigma+=aIter->second->stress_source(rgauss.col(k));
			}
			
			// 3- If ExternalStressFile is defined add the external stress
#ifdef UserStressFile
			sigma+=userStress(k); 
#endif
			
			if (shared.use_bvp){
				//	std::cout<<"BVP stress"<<std::endl;
				//	std::cout<<shared.domain.stressAt(rgauss.col(k))<<std::endl;
				sigma+=this->source->bvpStress*(1.0-Quadrature<1,qOrder,QuadratureRule>::abscissa(k))+this->sink->bvpStress*Quadrature<1,qOrder,QuadratureRule>::abscissa(k);
			}
			
			//return sigma;
			//			return sigma+shared.loadController.externalStress();
			return sigma+shared.externalStress;			
		}
		
		//////////////////////////////////
		//! PK force. Calculates the PK force at the k-th field point
		VectorDim pkForce(const size_t & k){
			//			return (stress_field(k)*Burgers).cross(rlgauss.col(k));
			//			return (shared.use_bvp) ? ((quadratureParticleVector[k].stress(this->source->bvpStress,this->sink->bvpStress)+shared.loadController.externalStress())*Burgers).cross(rlgauss.col(k))
			//			/*                   */ : ((quadratureParticleVector[k].stress()+shared.loadController.externalStress())*Burgers).cross(rlgauss.col(k));			
			return (shared.use_bvp) ? ((quadratureParticleVector[k].stress(this->source->bvpStress,this->sink->bvpStress)+shared.externalStress)*Burgers).cross(rlgauss.col(k))
			/*                   */ : ((quadratureParticleVector[k].stress()+shared.externalStress)*Burgers).cross(rlgauss.col(k));			
			
		}
		
		//////////////////////////////////
		//! get_pkGauss
		MatrixDimQorder get_pkGauss() const {
			return pkGauss;
		}
		
		/********************************************************/
		double energy() const {
			/*! The total elastic energy generated by this segment. 
			 *  Includes self-energy and interation energy between this segment and other segments in the network.
			 */
			double temp=0.0;
			Quadrature<1,qOrder,QuadratureRule>::integrate(this,temp,&LinkType::energy_field);
			
			return temp;
		}
		
		
		//////////////////////////////////
		//! energy_field. Calculates the elastic energy at the k-th quadrature point due to other dislocation segments.
		double energy_field(const int & k) const {
			double temp=0.0;
			for (AddressMapIteratorType aIter=this->ABbegin(); aIter!=this->ABend();++aIter){	
				temp+=aIter->second->energy_source(rgauss.col(k), rugauss.col(k), Burgers);
			}
			return temp;			
		}
		
		//////////////////////////////////
		double energy_source(const VectorDim & rf, const VectorDim & ruf, const VectorDim & bf) const {
			/*! The elastic energy generated by this (source) segment at the input field point
			 */
			double temp(0.0);
			Quadrature<1,qOrder,QuadratureRule>::integrate(this,temp,&LinkType::energy_integrand,rf,ruf,bf);			
			return -C2*temp;	
		}
		
		//////////////////////////////////
		//! energy_integrand
		//double energy_integrand(const int & k, const VectorDim & rf, const VectorDim & ruf, const VectorDim & bf) const {
		double energy_integrand(const int & k,  VectorDim & rf, const VectorDim & ruf, const VectorDim & bf) const {
			VectorDim DR= rf-rgauss.col(k);
			double RaSquared = DR.squaredNorm() + coreLsquared;
			return  (C1*(1.0+0.5*coreLsquared/RaSquared)*Burgers.dot(rugauss.col(k))*bf.dot(ruf)
					 +2.0*shared.material.nu*(1.0+0.5*coreLsquared/RaSquared)*(bf.dot(rugauss.col(k))*Burgers.dot(ruf))
					 -(Burgers.dot(bf)*(1.0+coreLsquared/RaSquared)+ Burgers.dot(DR)*bf.dot(DR)*1.0/RaSquared )*ruf.dot(rugauss.col(k))
					 )/std::pow(RaSquared,0.5);
		}
		
		//////////////////////////////////////////////////////////////
		//! assemble
		void assemble(){
			//std::cout<<"Thread "<<omp_get_thread_num()<< " assembling DislocationSegment "<<this->sID<<std::endl;
			assemble_Fq();
			assemble_Kqq();
		}
		
		
		
		
		//////////////////////////////////////////////////////////////
		// get_Kqq
		const MatrixNdof& get_Kqq() const {
			return Kqq;
		}
		
		//////////////////////////////////////////////////////////////
		// get_Fq
		const VectorNdof & get_Fq() const {
			return Fq;
		}
		
		//////////////////////////////////////////////////////////////
		//! plasticStrainRate
		MatrixDim plasticStrainRate() const {
//			VectorDim temp(VectorDim::Zero());
//			Quadrature<1,qOrder,QuadratureRule>::integrate(this,temp,&LinkType::plasticStrainRate_integrand);
			
			VectorDim V=(this->source->get_V().template segment<dim>(0)+this->sink->get_V().template segment<dim>(0))*0.5;	
			VectorDim temp(-V.cross(this->chord()));
			return (temp*Burgers.transpose() + Burgers*temp.transpose())*0.5;
		}
		
		VectorDim plasticStrainRate_integrand(const int& k) const {
			//! !!! CHANGE THIS!!! VELOCITY SHUOLD BE CALCULATED WITH SHAPE FUNCTIONS!!
//			VectorDim V=(this->source->velocity.template segment<dim>(0)+this->sink->velocity.template segment<dim>(0))*0.5;	
			VectorDim V=(this->source->get_V().template segment<dim>(0)+this->sink->get_V().template segment<dim>(0))*0.5;	
			return -V.cross(rugauss.col(k));
		}
		
		//////////////////////////////////////////////////////////////
		//! plasticStrainRate
		bool is_boundarySegment() const {
			return (this->source->nodeMeshLocation == onMeshBoundary && this->sink->nodeMeshLocation == onMeshBoundary);
		}
		
		
		/********************************************************/
		VectorDim displacement(const VectorDim & Rfield, const VectorDim& S) const {
			/*! The infinitesimal dispacement field generated by this segment at a field point
			 * @param[in] k			the current quadrature point
			 * @param[in] Rfield	the field point
			 *
			 * The return value is calculated according to:
			 *	\f[
			 *		\mathbf{u}(\mathbf{r}_f)=\frac{1}{8\pi(1-\nu)}\int_0^1 \mathbf{u}^*(\alpha) d\alpha
			 *	\f]
			 */
			VectorDim u_source = VectorDim::Zero();
			Quadrature<1,qOrder,QuadratureRule>::integrate(this,u_source,&LinkType::displacement_integrand,Rfield,S);
			return C4*u_source;		
		}
		
		/********************************************************/
		VectorDim displacement_integrand(const int & k, const VectorDim & Rfield, const VectorDim& S) const{	
			/*! The infinitesimal dispacement field generated by this segment at a field point
			 * @param[in] k			the current quadrature point
			 * @param[in] Rfield	the field point
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
			 *  The calculation of the solid angle is performed transforming the surface integral to a line integral. From Stokes theorem:
			 *	\f[
			 *  \oint\frac{\hat{\mathbf{S}}\times\mathbf{R}}{R(R-\mathbf{S}\cdot\mathbf{R})}\cdot d\mathbf{l}
			 *  \underbrace{=}_{\mbox{Stokes Th.}}\int\left(\nabla\times\frac{\hat{\mathbf{S}}\times\mathbf{R}}{R(R-\mathbf{S}\cdot\mathbf{R})}\right)\cdot\hat{\mathbf{n}}dA
			 *  \underbrace{=}_{\mbox{identity}}
			 * -\int\frac{\mathbf{R}}{R^3}\cdot\mathbf{n}dA 
			 *	\f]
			 * References:
			 * [1] Asvestas, J. Line integrals and physical optics. Part I. The transformation of the solid-angle surface integral to a line integral. J. Opt. Soc. Am. A, 2(6), 891–895.
			 */
			
			VectorDim R=Rfield-rgauss.col(k);	
			double Ra=std::pow(R.squaredNorm()+std::pow(coreL,2.0),0.5);			
			return 1.0/Ra * (+ 2.0*C1/(Ra+R.dot(S))*Burgers*(S.cross(R)).dot(rugauss.col(k))
			/*            */ + C3*rugauss.col(k).cross(Burgers) 
			/*            */ + 1.0/std::pow(Ra,2)*(rugauss.col(k).cross(Burgers)).dot(R)*R
			/*            */ );
		}
		
		/********************************************************/
		MatrixDim lattice_rotation_source(const VectorDim & Rfield) const {
			/*! The lattice rotation tensor generated by this dislocatio segment at   
			 * @param[in] Rfield the vector representing the filed point
			 */
			
			MatrixDim dg(displacement_gradient_source(Rfield));
			return 0.5*(dg.transpose()-dg); 
		}
		
		/********************************************************/
		MatrixDim displacement_gradient_source(const VectorDim & Rfield) const {
			/*! The lattice rotation tensor generated by this dislocatio segment at   
			 * @param[in] Rfield the vector representing the filed point
			 */
			MatrixDim dis_gradient_source = MatrixDim::Zero();
			Quadrature<1,qOrder,QuadratureRule>::integrate(this,dis_gradient_source,&LinkType::displacement_gradient_integrand,Rfield);
			
			return (C4 * dis_gradient_source); 
		}
		
		
		
		/********************************************************/
		MatrixDim displacement_gradient_integrand(const int & k, const VectorDim& Rfield) const {
			
			VectorDim R=Rfield-rgauss.col(k);
			double RaSquared = R.squaredNorm() + coreLsquared;
			
			return  ( C1*( 2.0e00 + (3.0e00*coreLsquared/RaSquared) ) * (Burgers*(rugauss.col(k).cross(R)).transpose() + (Burgers.cross(rugauss.col(k)))*R.transpose())
					 + (rugauss.col(k).cross(Burgers))*R.transpose() + R * (rugauss.col(k).cross(Burgers)).transpose()
					 + R.cross(Burgers).dot(rugauss.col(k)) * (3.0/RaSquared*R*R.transpose() - I) ) / std::pow(RaSquared,1.5);
		}
		
		
		//////////////////////////////////////////////////////////////
		//HermiteCoefficient: Hermite coefficients (uniform parametrization)
		Eigen::Matrix<double,dim-1,Ncoeff> hermiteLocalCoefficient() const {
			const MatrixDim G2L(DislocationLocalReference<dim>::global2local(this->chord(),this->glidePlaneNormal));
			Eigen::Matrix<double,dim-1,Ncoeff> HrCf = Eigen::Matrix<double,dim-1,Ncoeff>::Zero();
			HrCf.col(1)= (G2L*this->sourceT()*this->chordParametricLength()).template segment<dim-1>(0);
			HrCf.col(2)= (G2L*(this->sink->get_P()-this->source->get_P())).template segment<dim-1>(0);
			HrCf.col(3)= (G2L*this->sinkT()*this->chordParametricLength()).template segment<dim-1>(0);
			return HrCf;
		}
		
		//////////////////////////////////////////////////////////////
		//hermite2Coef
		Eigen::Matrix<double,dim-1,Ncoeff> polynomialLocalCoeff() const {//change name polynomialCoeff
			return Coeff2Hermite<pOrder>::h2c<dim-1>(hermiteLocalCoefficient());
		}
		
		
		/********************************************************/
		std::map<double,VectorDim> boundaryCollision() const {
			//			std::set<double> boundaryCollision() const {
			
			/*! The collision points with the glidePlane boundary polygon.
			 */
			std::map<double,VectorDim> temp;
			//	std::set<double> temp;
			
			
			if (shared.boundary_type){
				VectorDim origin = this->source->get_P();
				
				//				const GlidePlane<dim,MaterialType>* const p_glide = &glideMap.at(origin.dot(this->glidePlaneNormal)*this->glidePlaneNormal);
				
				Eigen::Matrix<double,dim-1,Ncoeff> C1L=polynomialLocalCoeff();
				const MatrixDim G2L(DislocationLocalReference<dim>::global2local(this->chord(),this->glidePlaneNormal));
				
				
				for (int k=0; k<pGlidePlane->segmentMeshCollisionPairContainer.size();++k){
					//					for (int k=0; k<p_glide->segmentMeshCollisionPairContainer.size();++k){
					
					Eigen::Matrix<double,dim-1,2> H2L;
					H2L.col(0)=(G2L*(pGlidePlane->segmentMeshCollisionPairContainer[k].first -origin)).template segment<dim-1>(0);
					H2L.col(1)=(G2L*(pGlidePlane->segmentMeshCollisionPairContainer[k].second-origin)).template segment<dim-1>(0);
					//dislcoament
					PlanarSplineImplicitization<pOrder> psi(Coeff2Hermite<pOrder>::template h2c<dim-1>(C1L));
					std::set<std::pair<double,double> > intersectinParameters = psi.template intersectWith<1>(Coeff2Hermite<1>::template h2c<dim-1>(H2L),FLT_EPSILON);
					
					for (std::set<std::pair<double,double> >::const_iterator iter=intersectinParameters.begin();iter!=intersectinParameters.end();++iter){
						temp.insert(std::make_pair(iter->first,this->get_r(iter->first)));
						//						temp.insert(iter->first);
						
					}
				}
			}
			
			return temp;
		}
		
		

		
		/*********************************************************************/		
		CrossSlipSegment<LinkType> isCrossSlipSegment(const double& sinThetaCrossSlipCr,const double& crossSlipLength) const {
			return CrossSlipSegment<LinkType>(*this,sinThetaCrossSlipCr,crossSlipLength);
		}

		/*********************************************************************/
		VectorDim conjugatePlaneNormal() const {
			return shared.material.conjugatePlaneNormal(Burgers,this->glidePlaneNormal);		
		}
		
		/*********************************************************************/
		VectorDim integralPK() const {
			return pkGauss.col(qOrder/2)*this->chord().norm();
		}

		/* intersectWith ******************************************************/
		std::set<std::pair<double,double> > intersectWith( const Derived* const p_other  , const double& tol=FLT_EPSILON) const {
			return DislocationSegmentIntersection<dim,pOrder>(this->hermiteCoefficients(),this->glidePlaneNormal).intersectWith(p_other->hermiteCoefficients(),p_other->glidePlaneNormal,tol,p_other);
		}
		
		
		/* friend T& operator << **********************************************/
		template <class T>
		friend T& operator << (T& os, const Derived& ds){
			os  << ds.source->sID<<"\t"<< ds.sink->sID<<"\t"
			/**/<< std::setprecision(15)<<std::scientific<<ds.flow.transpose()<<"\t"
			/**/<< std::setprecision(15)<<std::scientific<<ds.glidePlaneNormal.transpose()<<"\t"
			/**/<< ds.sourceTfactor<<"\t"
			/**/<< ds.sinkTfactor<<"\t"
			/**/<< ds.pSN()->sID;
			return os;
		}
		
		
	};
	
	
	
	template <short unsigned int dim, short unsigned int corder, typename InterpolationType,
	/*	   */ double & alpha, short unsigned int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule, 
	/*	   */ typename MaterialType>
	const Eigen::Matrix<double,dim,dim> DislocationSegment<dim,corder,InterpolationType,alpha,qOrder,QuadratureRule,MaterialType>::I=Eigen::Matrix<double,dim,dim>::Identity();
	
	
	//////////////////////////////////////////////////////////////s
} // namespace model
#endif


//		/*********************************************************************/		
//		bool isCrossSlipSegment(const double& sinThetaCrossSlipCr,const double& crossSlipLength) const {
//			return this->chord().norm()>1.1*crossSlipLength
//			/*  */ && this->chord().normalized().cross(Burgers.normalized()).norm()<=sinThetaCrossSlipCr
//			/*  */ && this->source->nodeMeshLocation != onMeshBoundary && this->sink->nodeMeshLocation != onMeshBoundary // not on the boundary
//			/*  */ && (2.25*(pkGauss.col(qOrder/2)-pkGauss.col(qOrder/2).dot(this->glidePlaneNormal)*this->glidePlaneNormal).squaredNorm()<
//			/*     */ (pkGauss.col(qOrder/2)-pkGauss.col(qOrder/2).dot(conjugatePlaneNormal())*conjugatePlaneNormal()).squaredNorm());
//		}

		/*********************************************************************/
//		std::pair<VectorDim,VectorDim> crossSlipPoints(const double& crossSlipLength) const {			
//			std::pair<VectorDim,VectorDim> temp;
//			if(Burgers.dot(this->chord())>=0.0){
//				temp=std::make_pair(this->get_r(0.5)-Burgers.normalized()*crossSlipLength*0.5,this->get_r(0.5)+Burgers.normalized()*crossSlipLength*0.5);
//			}
//			else{
//				temp=std::make_pair(this->get_r(0.5)+Burgers.normalized()*crossSlipLength*0.5,this->get_r(0.5)-Burgers.normalized()*crossSlipLength*0.5);
//			}
//			return temp;
//		}
		
//		VectorDim crossSlipConjugatePoint(const double& crossSlipLength) const {
//			std::pair<VectorDim,VectorDim> Ppair(crossSlipPoints(crossSlipLength));
//			VectorDim dir(Burgers.cross(conjugatePlaneNormal()).normalized());
//			double dirDotPK(dir.dot(pkGauss.col(qOrder/2)));
//			double sgnDir((dirDotPK > 0.0) ? 1.0 : ((dirDotPK < 0.0) ? -1.0 : 0.0));
//			return 0.5*(Ppair.first+Ppair.second) + sgnDir*dir*crossSlipLength*0.5;
//		}





//		bool isCrossSlipSegment(const double& sinThetaCrossSlipCr,const double& crossSlipLength) const {
//			return this->chord().norm()> 1.1*crossSlipLength
//			/*  */ && this->chord().normalized().cross(Burgers.normalized()).norm()<=sinThetaCrossSlipCr
////			/*  */ && this->chord().norm()>3.0*crossSlipLength
//			/*  */ && (2.25*(pkGauss.col(qOrder/2)-pkGauss.col(qOrder/2).dot(this->glidePlaneNormal)*this->glidePlaneNormal).squaredNorm()<
//			/*           */ (pkGauss.col(qOrder/2)-pkGauss.col(qOrder/2).dot(conjugatePlaneNormal())*conjugatePlaneNormal()).squaredNorm());
//		}



//		//////////////////////////////////////////////////////////////
//		// get_Kqq
//		MatrixNdof get_Kqq() const {
//			MatrixNdof KqqEx=MatrixNdof::Zero();
//			for (int i=0;i<Ncoeff;++i){
//				for (int j=0;j<Ncoeff;++j){
//					KqqEx.template block<dim,dim>(dim*i,dim*j).setIdentity()*=Kqq(i,j);
//				}
//			}
//			return KqqEx;
//		}



//		typedef Eigen::Matrix<double,Ncoeff,Ncoeff> MatrixNcoeff;
//		MatrixNcoeff Kqq;

//! Expanded Segment Stiffness Matrix 
//typedef Eigen::Matrix<double,Ndof,Ndof>	MatrixNdof;
//MatrixNdof KqqEx;



//S=this->glidePlaneNormal ;
//update();

//S=this->glidePlaneNormal ;
//update();



//			quadratureParticleVector.push_back(new DislocationQuadratureParticleType(k,rgauss.col(k),rugauss.col(k),Burgers,this->source,this->sink));


//	template <short unsigned int dim, short unsigned int corder, typename InterpolationType,
//	/*	   */ double & alpha, short unsigned int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule, 
//	/*	   */ typename MaterialType>	
//	boost::ptr_map<Eigen::Matrix<double,dim,1>,GlidePlane<dim,MaterialType>,CompareVectorsByComponent<double,dim,float>  > DislocationSegment<dim,corder,InterpolationType,alpha,qOrder,QuadratureRule,MaterialType>::glideMap;


//			//shared.material.find_slipSystem(this->chord(),Burgers,allowedSlipSystems);
//			
//			//		shared.material.find_slipSystem(this->chord(),unitBurgers,allowedSlipSystems);
//			VectorDim npl=this->source->get_P().dot(this->glidePlaneNormal )*this->glidePlaneNormal ;
//			
//			
//			typename boost::ptr_map<VectorDim,GlidePlane<dim,MaterialType>,CompareVectorsByComponent<double,dim,float>  >::const_iterator glideIter=glideMap.find(npl);
//
//			//GlidePlane<dim,MaterialType>* p_glide;
//			
//			if (glideIter==glideMap.end()){ // not found
//				std::auto_ptr<GlidePlane<dim,MaterialType> > pSP (new GlidePlane<dim,MaterialType>(npl) );
//				glideMap.insert(npl,pSP);
//				//p_glide = &(glideMap.insert(npl,pSP).first->second);
//				
//			}
//			else{ //found
//				//p_glide = &(glideIter->second);
//			}


//			//shared.material.find_slipSystem(this->chord(),Burgers,allowedSlipSystems);
//			
//			//		shared.material.find_slipSystem(this->chord(),unitBurgers,allowedSlipSystems);
//			VectorDim npl=this->source->get_P().dot(this->glidePlaneNormal )*this->glidePlaneNormal;
////			if (npl.norm()<FLT_EPSILON){
////				
////			}
//			
//
//			typename boost::ptr_map<VectorDim,GlidePlane<dim,MaterialType>,CompareVectorsByComponent<double,dim,float>  >::const_iterator glideIter=glideMap.find(npl);
//			
//
//			//GlidePlane<dim,MaterialType>* p_glide;
//			
//			if (glideIter==glideMap.end()){ // not found
//				std::auto_ptr<GlidePlane<dim,MaterialType> > pSP (new GlidePlane<dim,MaterialType>(npl) );
//				glideMap.insert(npl,pSP);
//				//p_glide = &(glideMap.insert(npl,pSP).first->second);
//				
//			}
//			else{ //found
//				//p_glide = &(glideIter->second);
//			}


//#include <model/Geometry/Splines/Intersection/SplineLineIntersection.h>
//#include <model/Geometry/Splines/Intersection/PlanarSplineImplicitization.h>


//			return (quadratureParticleVector[k].stress(this->source->bvpStress,this->sink->bvpStress)*Burgers).cross(rlgauss.col(k));


//        /* normal *****************************************************/
//		VectorDim this->glidePlaneNormal  const { // MARKED FOR REMOVAL: use const planeNormal known at constructor time
//			return allowedSlipSystems.begin()->normal;		
//		}



//		/* get_Burgers ************************************************/
//		const double & get_Burgers(const size_t & i) const{ // MARKED FOR REMOVAL
//			return Burgers(i);
//		}


//			MatrixDim R = G2L();							// THIS SHOULD BE CONSTANT AND STORED BY CONSTRUCTOR
//			R.col(0)= this->flow.normalized();
//			R.col(2)= this->glidePlaneNormal ;
//			R.col(1)= R.col(2).cross(R.col(0));
//HrCf.col(0)<< this->source->get_P();




//			Eigen::Matrix<double,dim-1,Ncoeff> h2c;
//			Eigen::Matrix<double,dim-1,Ncoeff> Hermite = hermiteLocalCoefficient(R);
//			
//			h2c.col(3)=-2.0*(Hermite.col(2)-Hermite.col(0))+(    Hermite.col(1)+Hermite.col(3));
//			h2c.col(2)= 3.0*(Hermite.col(2)-Hermite.col(0))-(2.0*Hermite.col(1)+Hermite.col(3));
//			h2c.col(1)= Hermite.col(1);
//			h2c.col(0)= Hermite.col(0);				
//			return h2c;


//		/* normal *****************************************************/
//		MatrixDim global2local() const{ // MARKED FOR REMOVAL: initialize by constructor time
//			MatrixDim R;							// THIS SHOULD BE CONSTANT AND STORED BY CONSTRUCTOR
//			R.col(0)= this->flow.normalized();
//			R.col(2)= this->glidePlaneNormal ;
//			R.col(1)= R.col(2).cross(R.col(0));
//			return R;
//		}


/////////////////////////////
// Declare static data member
//	template <short unsigned int dim, short unsigned int corder, typename InterpolationType,
//	/*	   */ double & alpha, short unsigned int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule, 
//	/*	   */ typename MaterialType>
//	MaterialType DislocationSegment<dim,corder,InterpolationType,alpha,qOrder,QuadratureRule,MaterialType>::material;




//		//////////////////////////////////
//		// update
//		void update(){ // WRONG THIS IS UPDATED EVERY TIME THE NODE DOFs CHANGE!!! SHOULD ONLY BE UPDATED BEFORE ASSEMBLY!!!
//			
////			std::cout<<"Dislocation Segment "<<this->source->sID<<"->"<<this->sink->sID<<std::endl;
////			std::cout<<"I'm here 0"<<std::endl;
////			quadratureParticleVector.clear();
////			std::cout<<"I'm here 1"<<std::endl;			
////			Quadrature<1,qOrder>::execute(this,&LinkType::updateQuadGeometryKernel); // this sould be called by the network before solving (both DD and bvp)
//			
//		}




//
//			for (int k=0;k<qOrder;++k){
//				
//				
//				temp+= energy_integrand(k,rf,RaSquared,ruf,bf)*weight(k); 
//			}


//			for (int k=0;k<qOrder;++k){
//				temp=*weight(k);
//			}


//			for (int k=0;k<qOrder;++k){
//				temp+=energy_field(k)*weight(k);
//			}


//make_SF(u);
//

//			rugauss.col(k)=ru;
//			ruugauss.col(k)=ruu;
//			rlgauss.col(k)=rl;
//			rllgauss.col(k)=rll;
//			jgauss(k)=j;
//			kappagauss(k)=kappa;

//	rgauss.setZero();
//	rugauss.setZero();
//	ruugauss.setZero();

//	rlgauss.setZero();
//	rllgauss.setZero();

//	jgauss.setZero();
//	kappagauss.setZero();


//! Initialize average values
//	arcL=0.0;
//	kappam=0.0;
//	rm.setZero();
//	rlm.setZero();
//	rllm.setZero();


//double jw;

//	for(int k=0;k<qOrder;++k){
//				make_all(k); // Calculates N, r, j, rl, rll for current value of u=abscissa(k)
//				double jw=j*this->weight(k);
//				
//				arcL+=jw;
//				rm+=r*jw;
//				rlm+=rl*jw;
//				rllm+=rll*jw;
//				kappam+=kappa*jw;
//			}
//			
//			rm/=arcL;
//			rlm/=arcL;	
//			rllm/=arcL;	
//			kappam/=arcL;
