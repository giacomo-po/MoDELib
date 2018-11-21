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

#include <PlanarDislocationSegment.h>

//#include <DislocationQuadraturePoint.h>


#ifndef NDEBUG
#define VerboseDislocationSegment(N,x) if(verboseDislocationSegment>=N){model::cout<<x;}
#else
#define VerboseDislocationSegment(N,x)
#endif


namespace model
{
    
    /**************************************************************************/
    /**************************************************************************/
    template <int dim, short unsigned int corder, typename InterpolationType>
    struct DislocationSegment : public PlanarDislocationSegment<DislocationSegment<dim,corder,InterpolationType>>
    {
        
        
        typedef typename TypeTraits<DislocationSegment<dim,corder,InterpolationType>>::NodeType NodeType;

        
        typedef PlanarDislocationSegment<DislocationSegment<dim,corder,InterpolationType>> PlanarDislocationSegmentType;
        
        static constexpr int Ncoeff=PlanarDislocationSegmentType::Ncoeff;
        static constexpr int pOrder=PlanarDislocationSegmentType::pOrder;

        typedef typename PlanarDislocationSegmentType::VectorDim VectorDim;
        typedef typename PlanarDislocationSegmentType::MatrixDim MatrixDim;
        typedef DislocationParticle<dim> DislocationParticleType;


        /******************************************************************/
        DislocationSegment(const std::shared_ptr<NodeType>& nI,
                           const std::shared_ptr<NodeType>& nJ) :
        /* init */ PlanarDislocationSegmentType(nI,nJ)
        {

        }
        

        
//        /**********************************************************************/
//        Eigen::Matrix<double,dim-1,Ncoeff> hermiteLocalCoefficient() const __attribute__ ((deprecated))
//        {
//            const MatrixDim G2L(DislocationLocalReference<dim>::global2local(this->chord(),this->glidePlaneNormal()));
//            Eigen::Matrix<double,dim-1,Ncoeff> HrCf = Eigen::Matrix<double,dim-1,Ncoeff>::Zero();
//            HrCf.col(1)= (G2L*this->sourceT()*this->chordParametricLength()).template segment<dim-1>(0);
//            HrCf.col(2)= (G2L*(this->sink->get_P()-this->source->get_P())).template segment<dim-1>(0);
//            HrCf.col(3)= (G2L*this->sinkT()*this->chordParametricLength()).template segment<dim-1>(0);
//            return HrCf;
//        }
//        
//        /**********************************************************************/
//        Eigen::Matrix<double,dim-1,Ncoeff> polynomialLocalCoeff() const __attribute__ ((deprecated)) //change name polynomialCoeff
//        {
//            return Coeff2Hermite<pOrder>::template h2c<dim-1>(hermiteLocalCoefficient());
//        }
        

        
        /**********************************************************************/
        const MatrixDim& midPointStress() const __attribute__ ((deprecated))
        {/*!\returns The stress matrix for the centre point over this segment.*/
            //            return stressGauss[this->quadraturePoints().size()/2];
            
            return this->quadraturePoints().size()? quadraturePoint(this->quadraturePoints().size()/2).stress : MatrixDim::Zero();
            
        }
        
    };

}
#endif


//        /**********************************************************************/
//        void assemble() __attribute__ ((deprecated))
//        {/*!Computes the following:
//          * - edge stiffness matrix Kqq
//          *	\f[
//          *		\mathbf{K} = int_0^1 \mathbf{K}^*(u) du
//          *	\f]
//          * - nodal force vector Fq
//          *	\f[
//          *		\mathbf{K} = int_0^1 \mathbf{K}^*(u) du
//          *	\f]
//          * - edge-to-component matrix Mseg
//          */
//
//
//            if(!this->hasZeroBurgers())
//            {
//                assert(0 && "THIS FUNCTION MUST BE REWORKED");
//
////                Fq.setZero();
////                Kqq.setZero();
//
//                //                //! 1- Compute and store stress and PK-force at quadrature points
//                //                stressGauss.clear();
//                //                if(!this->network().use_bvp && this->isBoundarySegment())
//                //                {
//                //                    pkGauss.setZero(dim,this->quadraturePoints().size());
//                //                }
//                //                else
//                //                {
//                //                    for (unsigned int k=0;k<this->quadraturePoints().size();++k)
//                //                    {
//                //                        MatrixDim temp(quadratureParticleContainer[k]->stress());
//                //
//                //                        if(this->network().use_externalStress)
//                //                        {
//                //                            temp+=this->network().extStressController.externalStress(this->quadraturePoint(k).r);
//                //                        }
//                //
//                //                        if(this->network().use_bvp)
//                //                        {
//                //                            temp += this->network().bvpSolver.stress(this->quadraturePoint(k).r,this->source->includingSimplex());
//                //                        }
//                //
//                //                        for(const auto& sStraight : this->network().poly.grainBoundaryDislocations() )
//                //                        {
//                //                            temp+=sStraight.stress(this->quadraturePoint(k).r);
//                //                        }
//                //
//                //                        stressGauss.push_back(temp);
//                //                        pkGauss.col(k)=(stressGauss[k]*Burgers).cross(this->quadraturePoint(k).rl);
//                //                        vGauss.col(k)=getGlideVelocity(k);
//                //                    }
//                //                }
//                //
//                //
//                //                QuadratureDynamicType::integrate(this->quadraturePoints().size(),this,Fq,&LinkType::velocityIntegrand);
//                //
//                //
//                //                if(corder==0)
//                //                {
//                //                    const double L=this->chord().norm();
//                //                    Kqq<<L/3.0,    0.0,    0.0, L/6.0,    0.0,    0.0,
//                //                    /**/ 0.0,    L/3.0,    0.0,     0.0,L/6.0,    0.0,
//                //                    /**/ 0.0,    0.0,L/3.0,     0.0,    0.0,    L/6.0,
//                //                    /**/ L/6.0,    0.0,    0.0, L/3.0,    0.0,    0.0,
//                //                    /**/ 0.0,L/6.0,    0.0,     0.0,    L/3.0,    0.0,
//                //                    /**/ 0.0,    0.0,L/6.0,     0.0,    0.0,    L/3.0;
//                ////                    Kqq*=this->chord().norm();
//                //                }
//                //                else
//                //                {
//                //                    assert(0 && "WE NEED TO INCREASE qOrder FOR THE FOLLOWING INTEGRATION, SINCE EVEN FOR LINEAR SEGMENTS Kqq IS NOT INTEGRATED CORRECLTY FOR SMALL qOrder");
//                //                    QuadratureDynamicType::integrate(this->quadraturePoints().size(),this,Kqq,&LinkType::stiffness_integrand);
//                //                }
//                //
//                //                h2posMap=this->hermite2posMap();
//                //
//                //                Mseg.setZero(Ncoeff*dim,h2posMap.size()*dim);
//                //
//                //                size_t c=0;
//                //                for(const auto& pair : h2posMap)
//                //                {
//                //                    for(int r=0;r<Ncoeff;++r)
//                //                    {
//                //                        Mseg.template block<dim,dim>(r*dim,c*dim)=pair.second.first(r)*MatrixDim::Identity();
//                //                    }
//                //                    c++;
//                //                }
//
//            }
//
//        }


//        /**********************************************************************/
//        void updateQuadraturePoints(ParticleSystem<DislocationParticleType>& particleSystem) __attribute__ ((deprecated))
//        {/*! @param[in] particleSystem the ParticleSystem of DislocationParticle
//          *  Computes all geometric properties at the k-th quadrature point
//          */
//
//            assert(0 && "THIS FUNCTION MUST BE REWORKED");
//
////            this->quadraturePoints().updateQuadraturePoints(*this,this->quadPerLength);
////
////            quadratureParticleContainer.clear();
////            quadratureParticleContainer.reserve(this->quadraturePoints().size());
////
////            if(!this->hasZeroBurgers())
////            {
////
////                if(!this->isBoundarySegment())
////                {
////                    const bool enableField=isGlissile();
////                    for (unsigned int k=0;k<this->quadraturePoints().size();++k)
////                    {
////                        quadratureParticleContainer.push_back(particleSystem.addParticle(this->quadraturePoint(k).r,
////                                                                                         this->source->sID,this->sink->sID,k,
////                                                                                         this->quadraturePoint(k).ru,
////                                                                                         Burgers,
////                                                                                         QuadratureDynamicType::abscissa(this->quadraturePoints().size(),k),
////                                                                                         QuadratureDynamicType::weight(this->quadraturePoints().size(),k),
////                                                                                         true,enableField,  // stressSource enabled, stressField enabled,
////                                                                                         true,enableField,  //   dispSource enabled,   dispField enabled,
////                                                                                         true,enableField));//   enrgSource enabled,   enrgField enabled,
////                    }
////
////                    if(!this->network().use_bvp) // not using FEM correction
////                    {
////                        if(this->network().useVirtualExternalLoops)
////                        {
////
////                            // Place Quadrature-particles on P1->P2
////                            if(this->source->isBoundaryNode())
////                            {
////                                const VectorDim& P1(this->source->get_P());
////                                const VectorDim d21=-(this->source->bndNormal()-this->source->bndNormal().dot(this->glidePlaneNormal())*this->glidePlaneNormal()).normalized();
////
////                                const VectorDim P2(P1-d21*this->virtualSegmentDistance);
////                                //const VectorDim d21=(P1-P2).normalized();
////                                const size_t qOrder12=QuadPowDynamicType::lowerOrder(this->quadPerLength*this->virtualSegmentDistance);
////
////                                //                                const VectorDim d21=(P1-P2).normalized();
////
////                                for (unsigned int k=0;k<qOrder12;++k)
////                                {
////
////                                    particleSystem.addParticle(P2+d21*this->virtualSegmentDistance*QuadratureDynamicType::abscissa(qOrder12,k),
////                                                               this->source->sID,this->sink->sID,this->quadraturePoints().size()+k,
////                                                               d21*this->virtualSegmentDistance,
////                                                               Burgers,
////                                                               QuadratureDynamicType::abscissa(qOrder12,k),
////                                                               QuadratureDynamicType::weight(qOrder12,k),
////                                                               true,false,  // stressSource disabled, stressField enabled,
////                                                               true,false,   //   dispSource  enabled,   dispField enabled,
////                                                               true,false);
////                                }
////
////                            }
////
////                            // Place Quadrature-particles on P1->P2
////                            if(this->sink->isBoundaryNode())
////                            {
////                                const VectorDim& P1(this->sink->get_P());
////                                const VectorDim d12=(this->sink->bndNormal()-this->sink->bndNormal().dot(this->glidePlaneNormal())*this->glidePlaneNormal()).normalized();
////
////                                const VectorDim P2(P1+d12*this->virtualSegmentDistance);
////                                //const VectorDim d12=(P2-P1).normalized();
////                                const size_t qOrder12=QuadPowDynamicType::lowerOrder(this->quadPerLength*this->virtualSegmentDistance);
////
////                                //                                const VectorDim d21=(P1-P2).normalized();
////
////                                for (unsigned int k=0;k<qOrder12;++k)
////                                {
////
////                                    particleSystem.addParticle(P1+d12*this->virtualSegmentDistance*QuadratureDynamicType::abscissa(qOrder12,k),
////                                                               this->source->sID,this->sink->sID,this->quadraturePoints().size()+k,
////                                                               d12*this->virtualSegmentDistance,
////                                                               Burgers,
////                                                               QuadratureDynamicType::abscissa(qOrder12,k),
////                                                               QuadratureDynamicType::weight(qOrder12,k),
////                                                               true,false,  // stressSource disabled, stressField enabled,
////                                                               true,false,   //   dispSource  enabled,   dispField enabled,
////                                                               true,false);
////                                }
////
////                            }
////
////                        }
////                    }
////
////
////                }
////                else // boundary segment
////                {
////                    if(this->network().use_bvp) // using FEM correction
////                    {
////                        if(this->network().useVirtualExternalLoops)
////                        {
////                            for (unsigned int k=0;k<this->quadraturePoints().size();++k)
////                            {
////                                quadratureParticleContainer.push_back(particleSystem.addParticle(this->quadraturePoint(k).r,
////                                                                                                 this->source->sID,this->sink->sID,k,
////                                                                                                 this->quadraturePoint(k).ru,Burgers,
////                                                                                                 QuadratureDynamicType::abscissa(this->quadraturePoints().size(),k),
////                                                                                                 QuadratureDynamicType::weight(this->quadraturePoints().size(),k),
////                                                                                                 false,true,  // stressSource enabled, stressField enabled,
////                                                                                                 false,true,  //   dispSource enabled,   dispField enabled,
////                                                                                                 false,true));//   enrgSource enabled,   enrgField enabled,
////                            }
////
////                            const VectorDim P1(this->source->get_P());
////                            const VectorDim P2(P1+this->source->bndNormal()*this->virtualSegmentDistance);
////                            const VectorDim P3(this->sink->get_P());
////                            const VectorDim P4(P3+this->sink->bndNormal()*this->virtualSegmentDistance);
////
////                            const size_t qOrder12=QuadPowDynamicType::lowerOrder(this->quadPerLength*this->virtualSegmentDistance);
////
////                            // Place Quadrature-particles on P1->P2
////                            if(!this->source->isPureBoundaryNode())
////                            {
////                                const VectorDim d12=(P2-P1).normalized();
////
////                                for (unsigned int k=0;k<qOrder12;++k)
////                                {
////
////                                    particleSystem.addParticle(P1+d12*this->virtualSegmentDistance*QuadratureDynamicType::abscissa(qOrder12,k),
////                                                               this->source->sID,this->sink->sID,this->quadraturePoints().size()+k,
////                                                               d12*this->virtualSegmentDistance,
////                                                               Burgers,
////                                                               QuadratureDynamicType::abscissa(qOrder12,k),
////                                                               QuadratureDynamicType::weight(qOrder12,k),
////                                                               true,false,  // stressSource disabled, stressField enabled,
////                                                               true,false,   //   dispSource  enabled,   dispField enabled,
////                                                               true,false);
////                                }
////
////                            }
////
////                            // Place Quadrature-particles on P2->P4
////                            const double L24=(P4-P2).norm();
////                            const size_t qOrder24=QuadPowDynamicType::lowerOrder(this->quadPerLength*L24);
////                            const VectorDim d24=(P4-P2)/L24;
////                            for (unsigned int k=0;k<qOrder24;++k)
////                            {
////                                particleSystem.addParticle(P2+d24*L24*QuadratureDynamicType::abscissa(qOrder24,k),
////                                                           this->source->sID,this->sink->sID,this->quadraturePoints().size()+qOrder12+k,
////                                                           d24*L24,
////                                                           Burgers,
////                                                           QuadratureDynamicType::abscissa(qOrder24,k),
////                                                           QuadratureDynamicType::weight(qOrder24,k),
////                                                           true,false,  // stressSource disabled, stressField enabled,
////                                                           true,false,   //   dispSource  enabled,   dispField enabled,
////                                                           true,false);
////                            }
////
////                            // Place Quadrature-particles on P4->P3
////                            if(!this->sink->isPureBoundaryNode())
////                            {
////                                const VectorDim d43=(P3-P4).normalized();
////                                for (unsigned int k=0;k<qOrder12;++k)
////                                {
////
////                                    particleSystem.addParticle(P4+d43*this->virtualSegmentDistance*QuadratureDynamicType::abscissa(qOrder12,k),
////                                                               this->source->sID,this->sink->sID,this->quadraturePoints().size()+qOrder12+qOrder24+k,
////                                                               d43*this->virtualSegmentDistance,
////                                                               Burgers,
////                                                               QuadratureDynamicType::abscissa(qOrder12,k),
////                                                               QuadratureDynamicType::weight(qOrder12,k),
////                                                               true,false,  // stressSource disabled, stressField enabled,
////                                                               true,false,   //   dispSource  enabled,   dispField enabled,
////                                                               true,false);
////                                }
////
////                            }
////
////                        }
////                        else // bnd segment, with bvp, without virtual segments
////                        {
////                            for (unsigned int k=0;k<this->quadraturePoints().size();++k)
////                            {
////                                quadratureParticleContainer.push_back(particleSystem.addParticle(this->quadraturePoint(k).r,
////                                                                                                 this->source->sID,this->sink->sID,k,
////                                                                                                 this->quadraturePoint(k).ru,Burgers,
////                                                                                                 QuadratureDynamicType::abscissa(this->quadraturePoints().size(),k),
////                                                                                                 QuadratureDynamicType::weight(this->quadraturePoints().size(),k),
////                                                                                                 true,true,  // stressSource enabled, stressField enabled,
////                                                                                                 true,true,  //   dispSource enabled,   dispField enabled,
////                                                                                                 true,true));//   enrgSource enabled,   enrgField enabled,
////                            }
////                        }
////                    }
////                    else // bonudary segment without bvp, do not place quadrature particles
////                    {
////
////                    }
////                }
////            }
//        }



//        /**********************************************************************/
//        void addToSolidAngleJump(const VectorDim& Pf, const VectorDim& Sf, VectorDim& dispJump) const __attribute__ ((deprecated))
//        {
//            if(this->isBoundarySegment() && this->network().useVirtualExternalLoops)
//            {
//                // first triangle is P1->P2->P3, second triangle is P2->P4->P3
//                const VectorDim P1(this->source->get_P());
//                const VectorDim P2(P1+this->source->bndNormal()*this->virtualSegmentDistance);
//                const VectorDim P3(this->sink->get_P());
//                const VectorDim P4(P3+this->sink->bndNormal()*this->virtualSegmentDistance);
//
//                dispJump += this->burgers()*LineSimplexIntersection<dim>::lineTriangleIntersection(Pf,Sf,P1,P2,P3);
//                dispJump += this->burgers()*LineSimplexIntersection<dim>::lineTriangleIntersection(Pf,Sf,P2,P4,P3);
//            }
//        }

