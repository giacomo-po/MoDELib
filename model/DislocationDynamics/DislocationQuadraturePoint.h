/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationQuadraturePoint_H_
#define model_DislocationQuadraturePoint_H_

#include <Eigen/Dense>
#include <model/Geometry/Splines/SplineBase.h>

namespace model
{
    template<int dim,int corder>
    struct DislocationQuadraturePoint
    {
        
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        
        static constexpr int Ncoeff= SplineBase<dim,corder>::Ncoeff;
        static constexpr int pOrder= SplineBase<dim,corder>::pOrder;
        static constexpr int Ndof= SplineBase<dim,corder>::Ndof;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        typedef Eigen::Matrix<double,dim,Ndof> MatrixDimNdof;
        typedef Eigen::Matrix<double,Ncoeff,Ncoeff> MatrixNcoeff;
        typedef Eigen::Matrix<double,Ncoeff,dim> MatrixNcoeffDim;
        typedef   QuadratureDynamic<1,UniformOpen,1,2,3,4,5,6,7,8,16,32,64,128,256,512,1024> QuadratureDynamicType;
        typedef QuadPowDynamic<pOrder,UniformOpen,1,2,3,4,5,6,7,8,16,32,64,128,256,512,1024> QuadPowDynamicType;
        
        const size_t sourceID;
        const size_t sinkID;
        const int qID;
        const Eigen::Matrix<double,1,Ncoeff> SF;    // Spline shape-functions at this quadrature point
        const VectorDim r;                          // position
        const VectorDim ru;                         // parametric tangent dr/du with u in [0:1]
        const double j;                             // jacobian dl/du
        const VectorDim rl;                         // unit tangent dr/dl
        const double dL;                            // line length corresponding to this quadrature point
        
        MatrixDim stress;
        VectorDim pkForce;
        VectorDim glideVelocity;
        
        /**********************************************************************/
        MatrixDimNdof SFgaussEx() const
        { /*! The MatrixDimNdof matrix of shape functions at the k-th quadrature point
           */
            MatrixDimNdof temp(MatrixDimNdof::Zero());
            for (size_t n=0;n<Ncoeff;++n)
            {
                temp.template block<dim,dim>(0,n*dim)=MatrixDim::Identity()*SF(n);
            }
            return temp;
        }
        
        /**********************************************************************/
        template<typename LinkType>
        static VectorDim  getGlideVelocity(const LinkType& seg,
                                           const VectorDim& fPK,
                                           const MatrixDim& S,
                                           const VectorDim& rl,
                                           const double& dL
                                           )
        {
            VectorDim n(seg.glidePlaneNormal()); // plane normal
            VectorDim b(seg.burgers()); // Burgers vector
            VectorDim t(rl);            // tangent vector
            
            // Select right-handed normal whenever possible
            if(seg.loopLinks().size()==1)
            {// pick right-handed normal for n
                const typename LinkType::LoopLinkType& loopLink(**seg.loopLinks().begin());
                if(std::fabs(loopLink.loop()->slippedArea())>FLT_EPSILON)
                {
                    if(seg.source->sID!=loopLink.source()->sID)
                    {// NetworkLink and LoopLink are oriented in opposite direction
                        b*=-1.0;
                        t*=-1.0;
                    }
                }
            }
            
            VectorDim glideForce = fPK-fPK.dot(n)*n;
            double glideForceNorm(glideForce.norm());
            
            if(glideForceNorm<FLT_EPSILON && seg.network().use_stochasticForce)
            {
                glideForce=seg.chord().cross(n);
                glideForceNorm=glideForce.norm();
                if(glideForceNorm>FLT_EPSILON)
                {
                    glideForce/=glideForceNorm;
                }
            }
            
            VectorDim vv=VectorDim::Zero();
            if(glideForceNorm>FLT_EPSILON)
            {
                double v =seg.network().poly.mobility->velocity(S,
                                                                b,
                                                                t,
                                                                n,
                                                                seg.network().poly.T,
                                                                dL,
                                                                seg.network().get_dt(),
                                                                seg.network().use_stochasticForce);
                
                
                assert((seg.network().use_stochasticForce || v>= 0.0) && "Velocity must be a positive scalar");
                const bool useNonLinearVelocity=true;
                if(useNonLinearVelocity && v>FLT_EPSILON)
                {
                    v= 1.0-std::exp(-v);
                }
                
                vv= v * glideForce/glideForceNorm;
            }
            return vv;
        }
        
        /**********************************************************************/
        template<typename LinkType>
        DislocationQuadraturePoint(const LinkType& seg,
                                   const int& q,const int& qOrder,
                                   const MatrixNcoeff& SFCH,
                                   const MatrixNcoeffDim& qH) :
        /* init */ sourceID(seg.source->sID)
        /* init */,sinkID(seg.sink->sID)
        /* init */,qID(q)
        /* init */,SF(QuadPowDynamicType::uPow(qOrder).row(qID)*SFCH)
        /* init */,r(SF*qH)
        /* init */,ru(QuadPowDynamicType::duPow(qOrder).row(qID)*SFCH.template block<Ncoeff-1,Ncoeff>(1,0)*qH)
        /* init */,j(ru.norm())
        /* init */,rl(ru/j)
        /* init */,dL(j*QuadratureDynamicType::weight(qOrder,qID))
        /* init */,stress(MatrixDim::Zero())
        /* init */,pkForce(VectorDim::Zero())
        /* init */,glideVelocity(VectorDim::Zero())
        {
            
        }
        
        /**********************************************************************/
        template<typename LinkType>
        void updateForcesAndVelocities(const LinkType& seg)
        {
            pkForce=(stress*seg.burgers()).cross(rl);
            glideVelocity=getGlideVelocity(seg,pkForce,stress,rl,dL);
        }
        
        
    };
    
    template<int dim,int corder>
    class DislocationQuadraturePointContainer : public std::deque<DislocationQuadraturePoint<dim,corder>>
    
    {
        static constexpr int Ncoeff= SplineBase<dim,corder>::Ncoeff;
        static constexpr int Ndof= SplineBase<dim,corder>::Ndof;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef Eigen::Matrix<double,Ndof,1> VectorNdof;
        typedef Eigen::Matrix<double,Ndof,Ndof> MatrixNdof;
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        typedef Eigen::Matrix<double,dim,Ndof> MatrixDimNdof;
        typedef DislocationQuadraturePointContainer<dim,corder> DislocationQuadraturePointContainerType;
        typedef DislocationQuadraturePoint<dim,corder> DislocationQuadraturePointType;
        typedef std::deque<DislocationQuadraturePointType> BaseContainerType;
        typedef Eigen::Matrix<double,Ncoeff,Ncoeff> MatrixNcoeff;
        typedef Eigen::Matrix<double,Ncoeff,dim> MatrixNcoeffDim;
        typedef typename DislocationQuadraturePointType::QuadratureDynamicType QuadratureDynamicType;
        typedef typename DislocationQuadraturePointType::QuadPowDynamicType QuadPowDynamicType;
        
        /**********************************************************************/
        VectorNdof nodalVelocityLinearKernel(const int& k) const
        { /*!@param[in] k the current quadrature point
           *\returns The kernel N^T(k)*v(k)*j(k)
           */
            return quadraturePoint(k).SFgaussEx().transpose()*quadraturePoint(k).glideVelocity*quadraturePoint(k).j;
        }
        
        /**********************************************************************/
        MatrixNdof nodalVelocityBilinearKernel(const int& k) const
        { /*! @param[in] k the current quadrature point
           *  The stiffness matrix integrand evaluated at the k-th quadrature point.
           *	\f[
           *		\mathbf{K}^* = \mathbf{N}^T \mathbf{B} \mathbf{N} \frac{dl}{du}
           *	\f]
           */
            const MatrixDimNdof SFEx(quadraturePoint(k).SFgaussEx());
            return SFEx.transpose()*SFEx*quadraturePoint(k).j;
        }
        
        /**********************************************************************/
        VectorDim pkKernel(const int& k) const
        { /*!@param[in] k the current quadrature point
           *\returns dF_PK/du=dF_PK/dL*dL/du at quadrature point k, where
           * u in [0,1] is the spline parametrization
           */
            return quadraturePoint(k).pkForce*quadraturePoint(k).j;
        }
        
        /**********************************************************************/
        VectorDim glideVelocityKernel(const int& k) const
        {
            return quadraturePoint(k).glideVelocity*this->quadraturePoint(k).j;
        }
        
    public:
        
        /**********************************************************************/
        template<typename LinkType>
        void updateGeometry(const LinkType& seg,
                            const double& quadPerLength)
        {
            
            this->clear();
            
            if(    !seg.hasZeroBurgers()
               &&  !seg.isBoundarySegment()
               &&  !seg.isSessile())
            {
                const int order=QuadPowDynamicType::lowerOrder(quadPerLength*seg.chord().norm());
                const MatrixNcoeff  SFCH(seg.sfCoeffs());
                const MatrixNcoeffDim qH(seg.hermiteDofs());
                for(int q=0;q<order;++q)
                {
                    this->emplace_back(seg,q,order,SFCH,qH);
                }
            }
        }
        
        /**********************************************************************/
        const DislocationQuadraturePointContainerType& quadraturePoints() const
        {
            return *this;
        }
        
        /**********************************************************************/
        DislocationQuadraturePointContainerType& quadraturePoints()
        {
            return *this;
        }
        
        /**********************************************************************/
        const DislocationQuadraturePointType& quadraturePoint(const int& k) const
        {
            return this->operator[](k);
        }
        
        /**********************************************************************/
        template<typename LinkType>
        void updateForcesAndVelocities(const LinkType& seg,
                                       const double& quadPerLength,
                                       const std::deque<StressStraight<dim>,Eigen::aligned_allocator<StressStraight<dim>>>& straightSegmentsDeq)
        {
            updateGeometry(seg,quadPerLength);
            
            if(this->size())
            {
                // Compute stress of straightSegmentsDeq
                const double L0(seg.chord().norm());
                const VectorDim c(0.5*(seg.source->get_P()+seg.sink->get_P()));
                for(const auto& ss : straightSegmentsDeq)
                {
                    SegmentSegmentDistance<dim> ssd(ss.P0,ss.P1,
                                                    seg.source->get_P(),seg.sink->get_P());
                    
                    //                        const double dr(ssd.dMin/(L0+ss.length));
                    const double dr(ssd.dMin/(L0));
                    
                    if(dr<10.0)
                    {// full interaction
                        for (auto& qPoint : quadraturePoints())
                        {
                            qPoint.stress += ss.stress(qPoint.r);
                        }
                    }
                    else if(dr<100.0)
                    {// 2pt interpolation
                        const MatrixDim stressSource(ss.stress(seg.source->get_P()));
                        const MatrixDim stressSink(ss.stress(seg.sink->get_P()));
                        for (auto& qPoint : quadraturePoints())
                        {
                            const double u(QuadratureDynamicType::abscissa(this->size(),qPoint.qID));
                            qPoint.stress += (1.0-u)*stressSource+u*stressSink;
                        }
                    }
                    else
                    {// 1pt interpolation
                        const MatrixDim stressC(ss.stress(c));
                        for (auto& qPoint : quadraturePoints())
                        {
                            qPoint.stress += stressC;
                        }
                    }
                }
                
                // Add other stress contributions, and compute PK force
                for (auto& qPoint : quadraturePoints())
                {
                    // Add stress of ExternalLoadController
                    if(seg.network().use_externalStress)
                    {
                        qPoint.stress += seg.network().extStressController.externalStress(qPoint.r);
                    }
                    
                    // Add BVP stress
                    if(seg.network().use_bvp)
                    {
                        qPoint.stress += seg.network().bvpSolver.stress(qPoint.r,seg.source->includingSimplex());
                    }
                    
                    // Add GB stress
                    for(const auto& sStraight : seg.network().poly.grainBoundaryDislocations() )
                    {
                        qPoint.stress += sStraight.stress(qPoint.r);
                    }
                    
                    // Add EshelbyInclusions stress
                    for(const auto& inclusion : seg.network().eshelbyInclusions() )
                    {
                        for(const auto& shift : seg.network().periodicShifts)
                        {
                            qPoint.stress += inclusion.second.stress(qPoint.r-shift);
                        }
                    }
                    
                    qPoint.updateForcesAndVelocities(seg);
                }
            }
        }
        
        /**********************************************************************/
        VectorNdof nodalVelocityVector() const
        { /*\returns The segment-integrated nodal velocity vector int N^T*v dL
           */
            VectorNdof Fq(VectorNdof::Zero());
            QuadratureDynamicType::integrate(this->size(),this,Fq,&DislocationQuadraturePointContainerType::nodalVelocityLinearKernel);
            return Fq;
        }
        
        /**********************************************************************/
        template<typename LinkType>
        MatrixNdof nodalVelocityMatrix(const LinkType& seg) const
        {
            MatrixNdof Kqq(MatrixNdof::Zero());
            if(corder==0)
            {// Analytical integral can be performed
                const double L(seg.chord().norm());
                Kqq<<L/3.0,  0.0,  0.0, L/6.0,  0.0,  0.0,
                /**/   0.0,L/3.0,  0.0,   0.0,L/6.0,  0.0,
                /**/   0.0,  0.0,L/3.0,   0.0,  0.0,L/6.0,
                /**/ L/6.0,  0.0,  0.0, L/3.0,  0.0,  0.0,
                /**/   0.0,L/6.0,  0.0,   0.0,L/3.0,  0.0,
                /**/   0.0,  0.0,L/6.0,   0.0,  0.0,L/3.0;
            }
            else
            {// Numerical integral must be performed
                assert(0 && "WE NEED TO INCREASE qOrder FOR THE FOLLOWING INTEGRATION, SINCE EVEN FOR LINEAR SEGMENTS Kqq IS NOT INTEGRATED CORRECLTY FOR SMALL qOrder");
                assert(0 && "THIS MUST RETURN A NON_ZERO VALUE EVEN FOR SEGMENTS WHICH DONT HAAVE QUADRATURE POINTS TO CORRECTLY SOVE K*V=F");
                QuadratureDynamicType::integrate(this->size(),this,Kqq,&DislocationQuadraturePointContainerType::nodalVelocityBilinearKernel);
            }
            return Kqq;
        }
        
        /*************************************************************/
        VectorDim glideVelocityIntegral() const
        {/*!\returns The integral of the glide velocity over the segment.
          */
            VectorDim V(VectorDim::Zero());
            QuadratureDynamicType::integrate(this->quadraturePoints().size(),this,V,&DislocationQuadraturePointContainerType::glideVelocityKernel);
            return V;
        }
        
        /**********************************************************************/
        VectorDim pkIntegral() const
        {/*!\returns The integral of the PK force over the segment.
          */
            VectorDim F(VectorDim::Zero());
            QuadratureDynamicType::integrate(this->quadraturePoints().size(),this,F,&DislocationQuadraturePointContainerType::pkKernel);
            return F;
        }
        
        //        /**********************************************************************/
        //        double arcLength() const
        //        {
        //            return SplineSegmentType::template arcLength<16,UniformOpen>();
        //        }
        
        
    };
    
}
#endif
