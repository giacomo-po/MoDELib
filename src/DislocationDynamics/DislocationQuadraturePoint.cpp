/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

//Implements non-uniform open quadrature rule

#ifndef model_DislocationQuadraturePoint_cpp_
#define model_DislocationQuadraturePoint_cpp_

#include <Eigen/Dense>
#include <SplineBase.h>
#include <QuadratureDynamic.h>
#include <QuadPowDynamic.h>
#include <StressStraight.h>
#include <SegmentSegmentDistance.h>
#include <StraightDislocationSegment.h>
#include <EshelbyInclusionBase.h>
#include <DefectiveCrystalParameters.h>
#include <Polygon2D.h>
#include <CatmullRomSplineSegment.h>
#include <DislocationSegment.h>
#include <DislocationQuadraturePoint.h>

namespace model
{

template<int dim,int corder>
DislocationQuadraturePoint<dim,corder>::DislocationQuadraturePoint(const LinkType& parentSegment,
                                                                   const int& q,const int& qOrder,
                                                                   const MatrixNcoeff& SFCH,
                                                                   const MatrixNcoeffDim& qH) :
/* init */ sourceID(parentSegment.source->sID)
/* init */,sinkID(parentSegment.sink->sID)
/* init */,qID(q)
/* init */,SF(QuadPowDynamicType::uPow(qOrder).row(qID)*SFCH)
/* init */,r(SF*qH)
/* init */,ru(QuadPowDynamicType::duPow(qOrder).row(qID)*SFCH.template block<Ncoeff-1,Ncoeff>(1,0)*qH)
/* init */,j(ru.norm())
/* init */,rl(j > FLT_EPSILON ? (ru/j).eval() : VectorDim::Zero())
/* init */,dL(j*QuadratureDynamicType::weight(qOrder,qID))
/* init */,stress(MatrixDim::Zero())
/* init */,pkForce(VectorDim::Zero())
/* init */,stackingFaultForce(VectorDim::Zero())
/* init */,lineTensionForce(VectorDim::Zero())
/* init */,glideVelocity(VectorDim::Zero())
/* init */,elasticEnergyPerLength(0.0)
/* init */,coreEnergyPerLength(0.0)
/* init */,inclusionID(-1)
{
    for(const auto& inclusion : parentSegment.network().eshelbyInclusions() )
    {// Add EshelbyInclusions stress
        if(inclusion.second->contains(r))
        {
            inclusionID=inclusion.second->sID;
            break;
        }
    }
}

template<int dim,int corder>
DislocationQuadraturePoint<dim,corder>::DislocationQuadraturePoint(const size_t sourceID_in,
                                                                   const size_t sinkID_in,
                                                                   const int& q,const int& qOrder,
                                                                   const MatrixNcoeff& SFCH,
                                                                   const MatrixNcoeffDim& qH) :
/* init */ sourceID(sourceID_in)
/* init */,sinkID(sinkID_in)
/* init */,qID(q)
/* init */,SF(QuadPowDynamicType::uPow(qOrder).row(qID)*SFCH)
/* init */,r(SF*qH)
/* init */,ru(QuadPowDynamicType::duPow(qOrder).row(qID)*SFCH.template block<Ncoeff-1,Ncoeff>(1,0)*qH)
/* init */,j(ru.norm())
/* init */,rl(j > FLT_EPSILON ? (ru/j).eval() : VectorDim::Zero())
/* init */,dL(j*QuadratureDynamicType::weight(qOrder,qID))
/* init */,stress(MatrixDim::Zero())
/* init */,pkForce(VectorDim::Zero())
/* init */,stackingFaultForce(VectorDim::Zero())
/* init */,lineTensionForce(VectorDim::Zero())
/* init */,glideVelocity(VectorDim::Zero())
/* init */,elasticEnergyPerLength(0.0)
/* init */,coreEnergyPerLength(0.0)
/* init */,inclusionID(-1)
{
    
}


template<int dim,int corder>
DislocationQuadraturePoint<dim,corder>::DislocationQuadraturePoint() :
/* init */ sourceID(0)
/* init */,sinkID(0)
/* init */,qID(0)
/* init */,SF(Eigen::Matrix<double,1,Ncoeff>::Zero())
/* init */,r(VectorDim::Zero())
/* init */,ru(VectorDim::Zero())
/* init */,j(0.0)
/* init */,rl(VectorDim::Zero())
/* init */,dL(0.0)
/* init */,stress(MatrixDim::Zero())
/* init */,pkForce(VectorDim::Zero())
/* init */,stackingFaultForce(VectorDim::Zero())
/* init */,lineTensionForce(VectorDim::Zero())
/* init */,glideVelocity(VectorDim::Zero())
/* init */,elasticEnergyPerLength(0.0)
/* init */,coreEnergyPerLength(0.0)
/* init */,inclusionID(-1)
{
    
}


//
//        template <class T>
//        friend T& operator << (T& os, const DislocationQuadraturePoint<dim,corder>& dqp)
//        {
//            os  << dqp.sourceID<<" "<<dqp.sinkID<<" "<<dqp.qID<<" "
//            /**/<< std::setprecision(15)<<std::scientific<< dqp.SF<<" "
//            /**/<< dqp.r.transpose()<<" "
//            /**/<< dqp.ru.transpose()<<" "
//            /**/<< dqp.j<<" "
//            /**/<< dqp.pkForce.transpose()<<" "
//            /**/<< dqp.elasticEnergyPerLength;
//            return os;
//        }


template<int dim,int corder>
typename DislocationQuadraturePoint<dim,corder>::MatrixDimNdof DislocationQuadraturePoint<dim,corder>::SFgaussEx() const
{ /*! The MatrixDimNdof matrix of shape functions at the k-th quadrature point
   */
    MatrixDimNdof temp(MatrixDimNdof::Zero());
    for (size_t n=0;n<Ncoeff;++n)
    {
        temp.template block<dim,dim>(0,n*dim)=MatrixDim::Identity()*SF(n);
    }
    return temp;
}


//        template<typename LinkType>
template<int dim,int corder>
typename DislocationQuadraturePoint<dim,corder>::VectorDim  DislocationQuadraturePoint<dim,corder>::getGlideVelocity(const LinkType& parentSegment,
                                                                                                                     const VectorDim& r,
                                                                                                                     const VectorDim& fPK,
                                                                                                                     const MatrixDim& S,
                                                                                                                     const VectorDim& rl,
                                                                                                                     const double& dL,
                                                                                                                     const int& inclusionID
                                                                                                                     )
{
    // std::cout<<"parentSegment "<<parentSegment.tag()<<" has slipSystem Compare to nullPtr"<<(parentSegment.slipSystem()==nullptr)<<std::endl;
    
    if(parentSegment.slipSystem())
    {
        VectorDim n(parentSegment.glidePlaneNormal()); // plane normal
        VectorDim b(parentSegment.burgers()); // Burgers vector
        VectorDim t(rl);            // tangent vector
        
        // Select right-handed normal whenever possible
        if(parentSegment.loopLinks().size()==1)
        {// pick right-handed normal for n
            const typename LinkType::LoopLinkType& loopLink(**parentSegment.loopLinks().begin());
            if(std::fabs(loopLink.loop->slippedArea())>FLT_EPSILON)
            {
                if(parentSegment.source->sID!=loopLink.source->sID)
                {// NetworkLink and LoopLink are oriented in opposite direction
                    b*=-1.0;
                    t*=-1.0;
                }
            }
        }
        
        VectorDim glideForce = fPK-fPK.dot(n)*n;
        double glideForceNorm(glideForce.norm());
        
        if(glideForceNorm<FLT_EPSILON && parentSegment.network().stochasticForceGenerator)
        {
            glideForce=parentSegment.chord().cross(n);
            glideForceNorm=glideForce.norm();
            if(glideForceNorm>FLT_EPSILON)
            {
                glideForce/=glideForceNorm;
            }
        }
        
        VectorDim vv=VectorDim::Zero();
        if(glideForceNorm>FLT_EPSILON)
        {
            
            //                    double v =parentSegment.network().poly.mobility->velocity(S,b,t,n,
            double v =parentSegment.slipSystem()->mobility->velocity(S,b,t,n,
                                                                     parentSegment.network().poly.T,
                                                                     dL,parentSegment.network().simulationParameters.dt,parentSegment.network().stochasticForceGenerator);
            // std::cout<<"v is "<<v<<std::endl;
            if(v<0.0 && v>=-FLT_EPSILON)
            {
                v=0.0; // kill roundoff errors for small negative velocities
            }
            
            if(inclusionID>=0)
            {
                v*=parentSegment.network().eshelbyInclusions().at(inclusionID)->mobilityReduction;
            }
            
            assert((parentSegment.network().stochasticForceGenerator || v>= 0.0) && "Velocity must be a positive scalar");
            const bool useNonLinearVelocity=true;
            if(useNonLinearVelocity && v>FLT_EPSILON)
            {
                v= 1.0-std::exp(-v);
            }
            
            //                    for(const auto& inclusion : parentSegment.network().eshelbyInclusions() )
            //                    {// Add EshelbyInclusions stress
            //                        if(inclusion.second.contains(r))
            //                        {
            //                            v*=inclusion.second.mobilityReduction;
            //                        }
            //                    }
            
            
            vv= v * glideForce/glideForceNorm;
            // std::cout<<"vv is "<<vv.transpose()<<std::endl;
            
        }
        return vv;
        
    }
    else
    {
        return VectorDim::Zero();
    }
    
}

//
//        std::set<const EshelbyInclusion<dim>*> getInclusions(const std::map<size_t,EshelbyInclusion<dim>>& inclusionContainer,
//                                                            const VectorDim& x)
//        {
//            std::set<const EshelbyInclusion<dim>*> temp;
//            for(const auto& pair : inclusionContainer)
//            {
//                if(pair.second.cointains(x))
//                {
//                    temp.insert(&pair.second);
//                }
//            }
//            return temp;
//        }



//        template<typename LinkType>
template<int dim,int corder>
void DislocationQuadraturePoint<dim,corder>::updateForcesAndVelocities(const LinkType& parentSegment)
{
    pkForce=(stress*parentSegment.burgers()).cross(rl);
    //            glideVelocity=getGlideVelocity(parentSegment,r,pkForce+stackingFaultForce+lineTensionForce,stress,rl,dL,inclusionID);
    glideVelocity=getGlideVelocity(parentSegment,r,pkForce+stackingFaultForce+lineTensionForce,stress+forceToStress(stackingFaultForce+lineTensionForce,parentSegment),rl,dL,inclusionID);
}

//        template<typename LinkType>
template<int dim,int corder>
typename DislocationQuadraturePoint<dim,corder>::MatrixDim DislocationQuadraturePoint<dim,corder>::forceToStress(const VectorDim& force,const LinkType& parentSegment) const
{
    // std::cout<<" Curvature norm "<<cmSeg.get_rll(paramU).norm()<<std::endl;
    const double burgersNorm(parentSegment.burgers().norm());
    if(burgersNorm>FLT_EPSILON)
    {
        const double resolvedtensionStress (force.norm()/burgersNorm);
        const VectorDim unitBurgers(parentSegment.burgers()/burgersNorm);
        const MatrixDim tensionStress(resolvedtensionStress*(unitBurgers*parentSegment.glidePlaneNormal().transpose() + parentSegment.glidePlaneNormal()*unitBurgers.transpose()));
        const VectorDim eqPKForce ((tensionStress*parentSegment.burgers()).cross(rl));
        const double eqForceNorm(eqPKForce.norm());
        if (eqForceNorm>FLT_EPSILON)
        {
            if (eqPKForce.dot(force) > 0.0)
            {
                return tensionStress;
            }
            else
            {
                return -tensionStress;
            }
        }
    }
    return MatrixDim::Zero();
    
    //            else
    //            {
    //                return MatrixDim::Zero();
    //            }
}

template struct DislocationQuadraturePoint<3,0>;





template<int dim,int corder>
typename DislocationQuadraturePointContainer<dim,corder>::VectorNdof DislocationQuadraturePointContainer<dim,corder>::nodalVelocityLinearKernel(const int& k) const
{ /*!@param[in] k the current quadrature point
   *\returns The kernel N^T(k)*v(k)*j(k)
   */
    return quadraturePoint(k).SFgaussEx().transpose()*quadraturePoint(k).glideVelocity*quadraturePoint(k).j;
}


template<int dim,int corder>
typename DislocationQuadraturePointContainer<dim,corder>::MatrixNdof DislocationQuadraturePointContainer<dim,corder>::nodalVelocityBilinearKernel(const int& k) const
{ /*! @param[in] k the current quadrature point
   *  The stiffness matrix integrand evaluated at the k-th quadrature point.
   *	\f[
   *		\mathbf{K}^* = \mathbf{N}^T \mathbf{B} \mathbf{N} \frac{dl}{du}
   *	\f]
   */
    const MatrixDimNdof SFEx(quadraturePoint(k).SFgaussEx());
    return SFEx.transpose()*SFEx*quadraturePoint(k).j;
}


template<int dim,int corder>
typename DislocationQuadraturePointContainer<dim,corder>::VectorDim DislocationQuadraturePointContainer<dim,corder>::pkKernel(const int& k) const
{ /*!@param[in] k the current quadrature point
   *\returns dF_PK/du=dF_PK/dL*dL/du at quadrature point k, where
   * u in [0,1] is the spline parametrization
   */
    return quadraturePoint(k).pkForce*quadraturePoint(k).j;
}


//Added by Yash
template<int dim,int corder>
typename DislocationQuadraturePointContainer<dim,corder>::MatrixDim DislocationQuadraturePointContainer<dim,corder>::stressKernel(const int& k) const
{ /*!@param[in] k the current quadrature point
   *\returns d_sigma/du=d_sigma/dL*dL/du at quadrature point k, where
   * u in [0,1] is the spline parametrization
   */
    return quadraturePoint(k).stress*quadraturePoint(k).j;
}


template<int dim,int corder>
typename DislocationQuadraturePointContainer<dim,corder>::VectorDim DislocationQuadraturePointContainer<dim,corder>::glideVelocityKernel(const int& k) const
{
    return quadraturePoint(k).glideVelocity*this->quadraturePoint(k).j;
}

//        template<typename LinkType>
//        void computeMatrixStackingFaultForces(const LinkType& parentSegment)
//        {
//            const double eps=1.0e-2;
//            MatrixDim temp(MatrixDim::Zero());
//            if(parentSegment.isGlissile())
//            {// slipSystem must exist
//                assert(parentSegment.glidePlanes().size()==1);
//                const auto& glidePlane(*parentSegment.glidePlanes().begin());
//                for(const auto& loopLink : parentSegment.loopLinks())
//                {
//                    
//                    const auto& slipSystem();
//                    
//                    if(loopLink->loop->slipSystem())
//                    {// gamma surface must exist to perform force calculation
//                        if(loopLink->loop->slipSystem()->gammaSurface)
//                        {// gamma surface must exist to perform force calculation
//                            VectorDim outDir((loopLink->sink->get_P() - loopLink->source->get_P()).cross(loopLink->loop->rightHandedUnitNormal()));
//                            const double outDirNorm(outDir.norm());
//                            if(outDirNorm>FLT_EPSILON)
//                            {
//                                outDir/=outDirNorm;
//                                std::vector<std::pair<VectorDim,VectorDim>> qPointSlip(quadraturePoints().size(),std::make_pair(VectorDim::Zero(),VectorDim::Zero())); // accumulated b1 and b2 for each qPoint
//                                for(const auto& otherLoop: parentSegment.network().loops())
//                                {
//                    
//                                    const auto& otherSlipSystem(otherLoop.second.lock()->slipSystem());
//                                    const bool otherSlipSystemIsPartial(otherSlipSystem? otherSlipSystem->isPartial() : false);
//                                    
//                                    
//                                    if(otherLoop.second.lock()->slipSystem() && glidePlane==otherLoop.second.lock()->glidePlane.get())
//                                    {
//                                        if(otherLoop.second.lock()->slipSystem()->isPartial())
//                                        {// only partial dislocations will contribute to a change in gamma surface
//                                            const double nRdotnR(loopLink->loop->rightHandedUnitNormal().dot(otherLoop.second.lock()->rightHandedUnitNormal()));
//                                            if(std::fabs(nRdotnR)>FLT_EPSILON)
//                                            {
//                                                std::vector<Eigen::Matrix<double,dim-1,1>> otherLocalNodes; // local position of other loop on parentSegment's loop
//                                                for(const auto& otherLoopLink : otherLoop.second.lock()->linkSequence())
//                                                {
//                                                    otherLocalNodes.push_back((loopLink->loop->slipSystem()->gammaSurface->G2L*(otherLoopLink->source->get_P()-glidePlane->P)).template segment<dim-1>(0));
//                                                }
//                                                
//                                                for(size_t q=0;q<quadraturePoints().size();++q)
//                                                {
//                                                    const auto& qPoint(quadraturePoints()[q]);
//                                                    const Eigen::Matrix<double,dim-1,1> x1((loopLink->loop->slipSystem()->gammaSurface->G2L*(qPoint.r + eps*outDir-glidePlane->P)).template segment<dim-1>(0));
//                                                    const Eigen::Matrix<double,dim-1,1> x2((loopLink->loop->slipSystem()->gammaSurface->G2L*(qPoint.r - eps*outDir-glidePlane->P)).template segment<dim-1>(0));
//                                                    
//                                                    const int wn1(Polygon2D::windingNumber(x1,otherLocalNodes));
//                                                    const int wn2(Polygon2D::windingNumber(x2,otherLocalNodes));
//                                                    qPointSlip[q].first -=wn1*otherLoop.second.lock()->burgers(); // slip vector is negative the burgers vector
//                                                    qPointSlip[q].second-=wn2*otherLoop.second.lock()->burgers(); // slip vector is negative the burgers vector
//                                                    
//                                                    //                                                        if(nRdotnR>FLT_EPSILON)
//                                                    //                                                        {// same rightHandedUnitNormal
//                                                    //                                                            qPointSlip[q].first +=wn1*otherLoop.second.lock()->burgers();
//                                                    //                                                            qPointSlip[q].second+=wn2*otherLoop.second.lock()->burgers();
//                                                    //                                                        }
//                                                    //                                                        else
//                                                    //                                                        {// opposite rightHandedUnitNormal
//                                                    //                                                            qPointSlip[q].first -=wn1*otherLoop.second.lock()->burgers();
//                                                    //                                                            qPointSlip[q].second-=wn2*otherLoop.second.lock()->burgers();
//                                                    //                                                        }
//                                                }
//                                            }
//                                        }
//                                    }
//                                }
//                                
//                                
//                                for(size_t q=0;q<quadraturePoints().size();++q)
//                                {
//                                    if((qPointSlip[q].first-qPointSlip[q].second).squaredNorm()>FLT_EPSILON)
//                                    {
//                                        const double gamma1(parentSegment.slipSystem()->misfitEnergy(qPointSlip[q].first));  // outer point
//                                        const double gamma2(parentSegment.slipSystem()->misfitEnergy(qPointSlip[q].second)); // inner point
//                                        quadraturePoints()[q].stackingFaultForce+= -(gamma2-gamma1)*outDir; // * fact
//                                    }
//                                }
//                                
//                            }
//                            
//                        }
//                    }
//                }
//            }
//        }




//        template<typename LinkType>
template<int dim,int corder>
void DislocationQuadraturePointContainer<dim,corder>::updateForcesAndVelocities(const LinkType& parentSegment,
                                                                                const double&,
                                                                                const double& )
{
    // updateQuadraturePoints(parentSegment,quadPerLength,isClimbStep);
    
    
    if(this->size())
    {
        if(parentSegment.network().computeDDinteractions)
        {
            const double L0(parentSegment.chord().norm());
            const VectorDim c(0.5*(parentSegment.source->get_P()+parentSegment.sink->get_P()));
            for(const auto& link : parentSegment.network().networkLinks())
            {
                if(   !link.second.lock()->hasZeroBurgers()
                   && !(link.second.lock()->isBoundarySegment() && parentSegment.network().simulationParameters.simulationType==DefectiveCrystalParameters::FINITE_NO_FEM) // exclude boundary segments even if they are non-zero Burgers
                   //                           && !(link.second.lock()->isVirtualBoundarySegment() && parentSegment.network().simulationParameters.simulationType==DefectiveCrystalParameters::PERIODIC_IMAGES)
                   )
                {
                    
                    const StraightDislocationSegment<dim>& ss(link.second.lock()->straight);
                    
                    
                    for(const auto& shift : parentSegment.network().periodicShifts)
                    {
                        SegmentSegmentDistance<dim> ssd(ss.P0,ss.P1,
                                                        parentSegment.source->get_P()+shift,parentSegment.sink->get_P()+shift);
                        
                        //                        const double dr(ssd.dMin/(L0+ss.length));
                        const double dr(ssd.dMin/(L0));
                        
                        if(dr<10.0)
                        {// full interaction
                            for (auto& qPoint : quadraturePoints())
                            {
                                qPoint.stress += ss.stress(qPoint.r+shift);
                            }
                        }
                        else if(dr<100.0)
                        {// 2pt interpolation
                            const MatrixDim stressSource(ss.stress(parentSegment.source->get_P()+shift));
                            const MatrixDim stressSink(ss.stress(parentSegment.sink->get_P()+shift));
                            for (auto& qPoint : quadraturePoints())
                            {
                                const double u(QuadratureDynamicType::abscissa(this->size(),qPoint.qID));
                                qPoint.stress += (1.0-u)*stressSource+u*stressSink;
                            }
                        }
                        else
                        {// 1pt interpolation
                            const MatrixDim stressC(ss.stress(c+shift));
                            for (auto& qPoint : quadraturePoints())
                            {
                                qPoint.stress += stressC;
                            }
                        }
                    }
                    
                    if(parentSegment.network().computeElasticEnergyPerLength)
                    {
                        for(const auto& shift : parentSegment.network().periodicShifts)
                        {
                            for (auto& qPoint : quadraturePoints())
                            {
                                qPoint.elasticEnergyPerLength += ss.elasticInteractionEnergy(qPoint.r+shift,qPoint.rl,parentSegment.burgers());
                            }
                        }
                    }
                    
                }
            }
        }
        
        
        // Stacking fault contribution in the matrix
        //                computeMatrixStackingFaultStress(parentSegment);
        //                computeMatrixStackingFaultForces(parentSegment);
        
        
        if(parentSegment.slipSystem() && parentSegment.glidePlanes().size()==1)
        {
            const auto& glidePlane(**parentSegment.glidePlanes().begin());
            const auto& slipSystem(*parentSegment.slipSystem());
            if(slipSystem.planeNoise)
            {
                for (auto& qPoint : quadraturePoints())
                {
                    const auto noiseVal(slipSystem.gridInterp(qPoint.r-glidePlane.P));
                    qPoint.stress += std::get<0>(noiseVal);
                }
            }
        }
        
        // Add other stress contributions, and compute PK force
        for (auto& qPoint : quadraturePoints())
        {
                        
            if(parentSegment.network().externalLoadController)
            {// Add stress of externalLoadController
                qPoint.stress += parentSegment.network().externalLoadController->stress(qPoint.r);
            }
            
            if(parentSegment.network().bvpSolver)
            {// Add BVP stress
                qPoint.stress += parentSegment.network().bvpSolver->stress(qPoint.r,parentSegment.source->includingSimplex());
            }
            
            //                    for(const auto& sStraight : parentSegment.network().poly.grainBoundaryDislocations() )
            //                    {// Add GB stress
            //                        qPoint.stress += sStraight.stress(qPoint.r);
            //                    }
            
            for(const auto& inclusion : parentSegment.network().eshelbyInclusions() )
            {// Add EshelbyInclusions stress
                for(const auto& shift : parentSegment.network().periodicShifts)
                {
                    qPoint.stress += inclusion.second->stress(qPoint.r+shift);
                }
            }
            
            //Add line tension contribution
            if (parentSegment.network().useLineTension && !parentSegment.hasZeroBurgers() && parentSegment.chordLength()>FLT_EPSILON)
            {
                //Add the line tension contribution due to only non -zero segments
                for (const auto &ll : parentSegment.loopLinks())
                {
                    // To discuss what would happen for links connected to boundary??
                    //  parameter may be same or opposite to the segment
                    
                    assert(ll->source->periodicPrev()!=nullptr && "PeriodicPrev must exist from network link for line tension contribution");
                    assert(ll->sink->periodicNext()!=nullptr && "PeriodicNext must exist from network link for line tension contribution");
                    
                    const VectorDim prevNodePos(ll->prev->twin()? ll->prev->twin()->source->periodicPrev()->get_P()
                                                -ll->source->periodicPlanePatch()->shift
                                                +ll->prev->twin()->source->periodicPrev()->periodicPlanePatch()->shift
                                                : ll->source->periodicPrev()->get_P());
                    const VectorDim nextNodePos(ll->next->twin()? ll->next->twin()->  sink->periodicNext()->get_P()
                                                -ll->  sink->periodicPlanePatch()->shift
                                                +ll->next->twin()->  sink->periodicNext()->periodicPlanePatch()->shift
                                                : ll->  sink->periodicNext()->get_P());
                    CatmullRomSplineSegment<dim> cmSeg(prevNodePos,ll->source->get_P(),ll->sink->get_P(),nextNodePos);
                    const double alpha (parentSegment.network().alphaLineTension);
                    const double paramUTemp ((qPoint.r-parentSegment.source->get_P()).norm()/parentSegment.chordLength());
                    const double paramU (ll->source->networkNode==parentSegment.source ? paramUTemp : 1.0-paramUTemp);
                    const VectorDim llunitTangent(cmSeg.get_rl(paramU));
                    const double qPointEnergyDensity (alpha * parentSegment.network().poly.C2 * (ll->loop->burgers().squaredNorm()-parentSegment.network().poly.nu*(std::pow(llunitTangent.dot(ll->loop->burgers()),2))));
                    const VectorDim qPointForceVector (qPointEnergyDensity * cmSeg.get_rll(paramU));
                    qPoint.coreEnergyPerLength+=qPointEnergyDensity;
                    qPoint.lineTensionForce+=qPointForceVector;
                }
            }
            
            qPoint.updateForcesAndVelocities(parentSegment);
        }
        
    }
}



//        template<typename LinkType>
template<int dim,int corder>
void DislocationQuadraturePointContainer<dim,corder>::updateQuadraturePoints(const LinkType& parentSegment,
                                                                             const double& quadPerLength,
                                                                             const bool& isClimbStep)
{
    
    this->clear();
    
    if(    !parentSegment.hasZeroBurgers()
       &&  !parentSegment.isBoundarySegment()
       &&  (!parentSegment.isSessile() || isClimbStep || parentSegment.network().computeElasticEnergyPerLength)
       &&  !parentSegment.isVirtualBoundarySegment())
    {
        const int order=QuadPowDynamicType::lowerOrder(quadPerLength*parentSegment.chord().norm());
        const MatrixNcoeff  SFCH(parentSegment.sfCoeffs());
        const MatrixNcoeffDim qH(parentSegment.hermiteDofs());
        for(int q=0;q<order;++q)
        {
            this->emplace_back(parentSegment,q,order,SFCH,qH);
        }
    }
}


template<int dim,int corder>
const typename DislocationQuadraturePointContainer<dim,corder>::DislocationQuadraturePointContainerType& DislocationQuadraturePointContainer<dim,corder>::quadraturePoints() const
{
    return *this;
}


template<int dim,int corder>
typename DislocationQuadraturePointContainer<dim,corder>::DislocationQuadraturePointContainerType& DislocationQuadraturePointContainer<dim,corder>::quadraturePoints()
{
    return *this;
}


template<int dim,int corder>
const typename DislocationQuadraturePointContainer<dim,corder>::DislocationQuadraturePointType& DislocationQuadraturePointContainer<dim,corder>::quadraturePoint(const int& k) const
{
    return this->operator[](k);
}




// VectorNdof nodalVelocityVector() const
// { /*\returns The segment-integrated nodal velocity vector int N^T*v dL
//    */
//     VectorNdof Fq(VectorNdof::Zero());
//     QuadratureDynamicType::integrate(this->size(),this,Fq,&DislocationQuadraturePointContainerType::nodalVelocityLinearKernel);
//     return Fq;
// }

//        template<typename LinkType>
template<int dim,int corder>
typename DislocationQuadraturePointContainer<dim,corder>::VectorNdof DislocationQuadraturePointContainer<dim,corder>::nodalVelocityVector(const LinkType& parentSegment) const
{ /*\returns The segment-integrated nodal velocity vector int N^T*v dL
   */
    VectorNdof Fq(VectorNdof::Zero());
    if(parentSegment.chordLength()>FLT_EPSILON)
    {
        QuadratureDynamicType::integrate(this->size(),this,Fq,&DislocationQuadraturePointContainerType::nodalVelocityLinearKernel);
    }
    return Fq;
}


//        template<typename LinkType>
template<int dim,int corder>
typename DislocationQuadraturePointContainer<dim,corder>::MatrixNdof DislocationQuadraturePointContainer<dim,corder>::nodalVelocityMatrix(const LinkType& parentSegment) const
{
    MatrixNdof Kqq(MatrixNdof::Zero());
    if(corder==0)
    {// Analytical integral can be performed
        const double L(parentSegment.chordLength());
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
template<int dim,int corder>
typename DislocationQuadraturePointContainer<dim,corder>::VectorDim DislocationQuadraturePointContainer<dim,corder>::glideVelocityIntegral() const
{/*!\returns The integral of the glide velocity over the segment.
  */
    VectorDim V(VectorDim::Zero());
    QuadratureDynamicType::integrate(this->quadraturePoints().size(),this,V,&DislocationQuadraturePointContainerType::glideVelocityKernel);
    return V;
}


template<int dim,int corder>
typename DislocationQuadraturePointContainer<dim,corder>::VectorDim DislocationQuadraturePointContainer<dim,corder>::pkIntegral() const
{/*!\returns The integral of the PK force over the segment.
  */
    VectorDim F(VectorDim::Zero());
    QuadratureDynamicType::integrate(this->quadraturePoints().size(),this,F,&DislocationQuadraturePointContainerType::pkKernel);
    return F;
}


//Added by Yash
template<int dim,int corder>
typename DislocationQuadraturePointContainer<dim,corder>::MatrixDim DislocationQuadraturePointContainer<dim,corder>::stressIntegral() const
{/*!\returns The integral of the stress over the segment.
  */
    MatrixDim stressInt(MatrixDim::Zero());
    QuadratureDynamicType::integrate(this->quadraturePoints().size(),this,stressInt,&DislocationQuadraturePointContainerType::stressKernel);
    return stressInt;
}

template class DislocationQuadraturePointContainer<3,0>;


}
#endif
