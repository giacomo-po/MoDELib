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

#include <memory>
#include <set>
#include <algorithm>    // std::set_intersection, std::sort


#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Sparse>

#include <model/Quadrature/Quadrature.h>
#include <model/Quadrature/QuadPow.h>
#include <model/DislocationDynamics/DislocationNetworkTraits.h>
#include <model/Geometry/Splines/SplineSegment.h>
#include <model/DislocationDynamics/Materials/Material.h>
#include <model/DislocationDynamics/GlidePlanes/GlidePlaneObserver.h>
#include <model/DislocationDynamics/GlidePlanes/GlidePlane.h>
#include <model/Geometry/Splines/Coeff2Hermite.h>
#include <model/DislocationDynamics/ElasticFields/DislocationParticle.h>
#include <model/ParticleInteraction/ParticleSystem.h>
#include <model/DislocationDynamics/DislocationLocalReference.h>
#include <model/IO/UniqueOutputFile.h>
#include <model/DislocationDynamics/Polycrystals/GrainBoundary.h>
#include <model/Geometry/LineSimplexIntersection.h>
#include <model/DislocationDynamics/BoundingLineSegments.h>
#include <model/DislocationDynamics/BoundingLineSegments.h>
#include <model/Mesh/MeshPlane.h>
#include <model/DislocationDynamics/DislocationQuadraturePoint.h>

namespace model
{
    
    /**************************************************************************/
    /**************************************************************************/
    template <int _dim, short unsigned int _corder, typename InterpolationType>
    class DislocationSegment : public SplineSegment<DislocationSegment<_dim,_corder,InterpolationType>,_dim,_corder>,
    /*                      */ private std::set<const MeshPlane<_dim>*>,
    /*                      */ private std::set<const GrainBoundary<_dim>*>,
    /*                      */ private std::set<const Grain<_dim>*>,
    /*                      */ private BoundingLineSegments<_dim>,
    /*                      */ public DislocationQuadraturePointContainer<_dim,_corder>
    {
        
        
    public:
        
        static constexpr int dim=_dim; // make dim available outside class
        static constexpr int corder=_corder; // make dim available outside class
        typedef SplineSegmentBase<dim,corder> SplineSegmentBaseType;
        typedef DislocationSegment<dim,corder,InterpolationType> LinkType;
        typedef typename TypeTraits<LinkType>::LoopNetworkType NetworkType;
        typedef LoopLink<LinkType> LoopLinkType;
        typedef DislocationNode<dim,corder,InterpolationType> NodeType;
        typedef SplineSegment<LinkType,dim,corder> SplineSegmentType;
        static constexpr int Ncoeff=SplineSegmentBaseType::Ncoeff;
        static constexpr int Ndof=SplineSegmentType::Ndof;
        typedef typename SplineSegmentType::MatrixNdof MatrixNdof;
        typedef typename SplineSegmentType::VectorNdof VectorNdof;
        typedef typename SplineSegmentType::VectorDim VectorDim;
        typedef typename SplineSegmentType::MatrixDim MatrixDim;
        typedef typename SplineSegmentType::MatrixDimNdof MatrixDimNdof;
        typedef typename SplineSegmentType::MatrixNcoeff MatrixNcoeff;
        typedef typename SplineSegmentType::MatrixNcoeffDim MatrixNcoeffDim;
        static constexpr int pOrder=SplineSegmentType::pOrder;
        typedef Eigen::Matrix<double,Ncoeff,1>     VectorNcoeff;
        typedef DislocationQuadraturePoint<_dim,_corder> DislocationQuadraturePointType;
        typedef typename DislocationQuadraturePointType::QuadPowDynamicType QuadPowDynamicType;
        typedef typename DislocationQuadraturePointType::QuadratureDynamicType QuadratureDynamicType;
        typedef DislocationParticle<dim> DislocationParticleType;
        typedef std::vector<DislocationParticleType*> QuadratureParticleContainerType;
        typedef LatticeVector<dim> LatticeVectorType;
        typedef ReciprocalLatticeDirection<dim> ReciprocalLatticeDirectionType;
        typedef std::set<const GrainBoundary<dim>*> GrainBoundaryContainerType;
        typedef std::set<const Grain<dim>*> GrainContainerType;
        typedef GlidePlane<dim> GlidePlaneType;
        typedef MeshPlane<dim> MeshPlaneType;
        typedef std::set<const MeshPlaneType*> MeshPlaneContainerType;
        typedef typename TypeTraits<LinkType>::MeshLocation MeshLocation;
        
        
    private:
        
        std::map<size_t,
        /*    */ std::pair<VectorNcoeff,VectorDim>,
        /*    */ std::less<size_t>,
        /*    */ Eigen::aligned_allocator<std::pair<size_t, std::pair<VectorNcoeff,VectorDim>> >
        /*    */ > h2posMap;
        Eigen::Matrix<double, Ndof, Eigen::Dynamic> Mseg;
        MatrixNdof Kqq; //! Segment Stiffness Matrix
        VectorNdof Fq; //! Segment Nodal Force Vector
        VectorDim Burgers; //! The Burgers vector
        
    public:
        
        
        static const Eigen::Matrix<double,_dim,_dim> I;
        static const Eigen::Matrix<double,_dim,1> zeroVector;
        static double quadPerLength;
        static double virtualSegmentDistance;
        QuadratureParticleContainerType quadratureParticleContainer;
        
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        
        /******************************************************************/
        static void initFromFile(const std::string& fileName)
        {
            
            LinkType::alpha=TextFileParser(fileName).readScalar<double>("parametrizationExponent",true);
            assert((LinkType::alpha)>=0.0 && "parametrizationExponent MUST BE >= 0.0");
            assert((LinkType::alpha)<=1.0 && "parametrizationExponent MUST BE <= 1.0");
            
            quadPerLength=TextFileParser(fileName).readScalar<double>("quadPerLength",true);
            assert((LinkType::quadPerLength)>=0.0 && "quadPerLength MUST BE >= 0.0");
            
            virtualSegmentDistance=TextFileParser(fileName).readScalar<double>("virtualSegmentDistance",true);
            
        }
        
        /******************************************************************/
        DislocationSegment(const std::shared_ptr<NodeType>& nI,
                           const std::shared_ptr<NodeType>& nJ) :
        /* init */ SplineSegmentType(nI,nJ),
        /* init */ Burgers(VectorDim::Zero())
        //        /* init */ qOrder(QuadPowDynamicType::lowerOrder(quadPerLength*this->chord().norm()))
        {/*! Constructor with pointers to source and sink, and flow
          *  @param[in] NodePair_in the pair of source and sink pointers
          *  @param[in] Flow_in the input flow
          */
            
            //            pkGauss.setZero(dim,this->quadraturePoints().size()); // necessary if this is not assembled
        }
        
        /**********************************************************************/
        const MeshPlaneContainerType& meshPlanes() const
        {
            return *this;
        }
        
        /**********************************************************************/
        MeshPlaneContainerType& meshPlanes()
        {
            return *this;
        }
        
        /**********************************************************************/
        const GrainContainerType& grains() const
        {
            return *this;
        }
        
        /**********************************************************************/
        GrainContainerType& grains()
        {
            return *this;
        }
        
        /**********************************************************************/
        GrainBoundaryContainerType& grainBoundaries()
        {
            return *this;
        }
        
        /**********************************************************************/
        const GrainBoundaryContainerType& grainBoundaries() const
        {
            return *this;
        }
        
        /**********************************************************************/
        const BoundingLineSegments<dim>& boundingBoxSegments() const
        {
            return *this;
        }
        
        /**********************************************************************/
        BoundingLineSegments<dim>& boundingBoxSegments()
        {
            return *this;
        }
        
        /**********************************************************************/
        bool addMeshPlane(const MeshPlaneType& gp)
        {
            const bool success=meshPlanes().insert(&gp).second;
            if(success)
            {
                const bool sourceContained(gp.contains(this->source->get_P()));
                const bool   sinkContained(gp.contains(this->  sink->get_P()));
                if(!(sourceContained && sinkContained))
                {
                    model::cout<<"DislocationSegment "<<this->source->sID<<"->"<<this->sink->sID<<std::endl;
                    model::cout<<"sourceContained "<<sourceContained<<std::endl;
                    model::cout<<"  sinkContained "<<sinkContained<<std::endl;
                    assert(false &&  "Glide Plane does not contain source or sink");
                }
                boundingBoxSegments().updateWithMeshPlane(gp); // Update _boundingBoxSegments. This must be called before updateGlidePlaneIntersections
                grains().insert(&this->network().poly.grain(gp.regionIDs.first));    // Insert new grain in grainSet
                grains().insert(&this->network().poly.grain(gp.regionIDs.second));   // Insert new grain in grainSet
            }
            return success;
        }
        
        /**********************************************************************/
        void addLink(LoopLinkType* const pL)
        {
            SplineSegmentType::addLink(pL); // forward to base class
            
            // Modify Burgers vector
            if(pL->source()->sID==this->source->sID)
            {
                Burgers+=pL->flow().cartesian();
            }
            else
            {
                Burgers-=pL->flow().cartesian();
            }
            
            addMeshPlane(pL->loop()->glidePlane);
        }
        
        /**********************************************************************/
        void removeLink(LoopLinkType* const pL)
        {
            SplineSegmentType::removeLink(pL);  // forward to base class
            
            // Modify Burgers vector
            if(pL->source()->sID==this->source->sID)
            {
                Burgers-=pL->flow().cartesian();
            }
            else
            {
                Burgers+=pL->flow().cartesian();
            }
            
            
            meshPlanes().clear();
            boundingBoxSegments().clear();
            for(const auto& loopLink : this->loopLinks())
            {
                addMeshPlane(loopLink->loop()->glidePlane);
            }
            
            addGrainBoundaryPlanes();
            
        }
        
        /**********************************************************************/
        size_t addGrainBoundaryPlanes()
        {
            size_t addedGp=0;
            for(const auto& gb : this->source->grainBoundaries())
            {
                if(this->sink->grainBoundaries().find(gb)!=this->sink->grainBoundaries().end())
                {
                    grainBoundaries().insert(gb);
                    //                    const auto& gp(gb->glidePlanes().begin()->second);// HERE BEGIN IS TEMPORARY, UNTIL WE STORE THE GLIDE PLANE OF THE CSL AND DSCL
                    //                    addedGp+=addMeshPlane(*gp.get());
                    addedGp+=addMeshPlane(*gb);
                }
            }
            
            return addedGp;
        }
        
        /**********************************************************************/
        void addToStressStraight(std::deque<StressStraight<dim>,Eigen::aligned_allocator<StressStraight<dim>>>& straightSegmentsDeq) const
        {
            if(!hasZeroBurgers())
            {
                
                if(this->network().use_bvp) // using FEM correction
                {
                    assert(0 && "FINISH HERE");
                }
                else
                {// Don't add stress
                    if(!isBoundarySegment())
                    {
                        straightSegmentsDeq.emplace_back(this->source->get_P(),
                                                         this->sink->get_P(),
                                                         burgers());
                    }
                    else
                    {
                        
                    }
                }
            }
        }
        
        /**********************************************************************/
        void updateQuadraturePoints(ParticleSystem<DislocationParticleType>& particleSystem)
        {/*! @param[in] particleSystem the ParticleSystem of DislocationParticle
          *  Computes all geometric properties at the k-th quadrature point
          */
            
            this->quadraturePoints().updateGeometry(*this,quadPerLength);
            
            quadratureParticleContainer.clear();
            quadratureParticleContainer.reserve(this->quadraturePoints().size());
            
            if(!hasZeroBurgers())
            {
                
                if(!isBoundarySegment())
                {
                    const bool enableField=isGlissile();
                    for (unsigned int k=0;k<this->quadraturePoints().size();++k)
                    {
                        quadratureParticleContainer.push_back(particleSystem.addParticle(this->quadraturePoint(k).r,
                                                                                         this->source->sID,this->sink->sID,k,
                                                                                         this->quadraturePoint(k).ru,
                                                                                         Burgers,
                                                                                         QuadratureDynamicType::abscissa(this->quadraturePoints().size(),k),
                                                                                         QuadratureDynamicType::weight(this->quadraturePoints().size(),k),
                                                                                         true,enableField,  // stressSource enabled, stressField enabled,
                                                                                         true,enableField,  //   dispSource enabled,   dispField enabled,
                                                                                         true,enableField));//   enrgSource enabled,   enrgField enabled,
                    }
                    
                    if(!this->network().use_bvp) // not using FEM correction
                    {
                        if(this->network().use_virtualSegments)
                        {
                            
                            // Place Quadrature-particles on P1->P2
                            if(this->source->isBoundaryNode())
                            {
                                const VectorDim& P1(this->source->get_P());
                                const VectorDim d21=-(this->source->bndNormal()-this->source->bndNormal().dot(glidePlaneNormal())*glidePlaneNormal()).normalized();
                                
                                const VectorDim P2(P1-d21*virtualSegmentDistance);
                                //const VectorDim d21=(P1-P2).normalized();
                                const size_t qOrder12=QuadPowDynamicType::lowerOrder(quadPerLength*virtualSegmentDistance);
                                
                                //                                const VectorDim d21=(P1-P2).normalized();
                                
                                for (unsigned int k=0;k<qOrder12;++k)
                                {
                                    
                                    particleSystem.addParticle(P2+d21*virtualSegmentDistance*QuadratureDynamicType::abscissa(qOrder12,k),
                                                               this->source->sID,this->sink->sID,this->quadraturePoints().size()+k,
                                                               d21*virtualSegmentDistance,
                                                               Burgers,
                                                               QuadratureDynamicType::abscissa(qOrder12,k),
                                                               QuadratureDynamicType::weight(qOrder12,k),
                                                               true,false,  // stressSource disabled, stressField enabled,
                                                               true,false,   //   dispSource  enabled,   dispField enabled,
                                                               true,false);
                                }
                                
                            }
                            
                            // Place Quadrature-particles on P1->P2
                            if(this->sink->isBoundaryNode())
                            {
                                const VectorDim& P1(this->sink->get_P());
                                const VectorDim d12=(this->sink->bndNormal()-this->sink->bndNormal().dot(glidePlaneNormal())*glidePlaneNormal()).normalized();
                                
                                const VectorDim P2(P1+d12*virtualSegmentDistance);
                                //const VectorDim d12=(P2-P1).normalized();
                                const size_t qOrder12=QuadPowDynamicType::lowerOrder(quadPerLength*virtualSegmentDistance);
                                
                                //                                const VectorDim d21=(P1-P2).normalized();
                                
                                for (unsigned int k=0;k<qOrder12;++k)
                                {
                                    
                                    particleSystem.addParticle(P1+d12*virtualSegmentDistance*QuadratureDynamicType::abscissa(qOrder12,k),
                                                               this->source->sID,this->sink->sID,this->quadraturePoints().size()+k,
                                                               d12*virtualSegmentDistance,
                                                               Burgers,
                                                               QuadratureDynamicType::abscissa(qOrder12,k),
                                                               QuadratureDynamicType::weight(qOrder12,k),
                                                               true,false,  // stressSource disabled, stressField enabled,
                                                               true,false,   //   dispSource  enabled,   dispField enabled,
                                                               true,false);
                                }
                                
                            }
                            
                        }
                    }
                    
                    
                }
                else // boundary segment
                {
                    if(this->network().use_bvp) // using FEM correction
                    {
                        if(this->network().use_virtualSegments)
                        {
                            for (unsigned int k=0;k<this->quadraturePoints().size();++k)
                            {
                                quadratureParticleContainer.push_back(particleSystem.addParticle(this->quadraturePoint(k).r,
                                                                                                 this->source->sID,this->sink->sID,k,
                                                                                                 this->quadraturePoint(k).ru,Burgers,
                                                                                                 QuadratureDynamicType::abscissa(this->quadraturePoints().size(),k),
                                                                                                 QuadratureDynamicType::weight(this->quadraturePoints().size(),k),
                                                                                                 false,true,  // stressSource enabled, stressField enabled,
                                                                                                 false,true,  //   dispSource enabled,   dispField enabled,
                                                                                                 false,true));//   enrgSource enabled,   enrgField enabled,
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
                                                               this->source->sID,this->sink->sID,this->quadraturePoints().size()+k,
                                                               d12*virtualSegmentDistance,
                                                               Burgers,
                                                               QuadratureDynamicType::abscissa(qOrder12,k),
                                                               QuadratureDynamicType::weight(qOrder12,k),
                                                               true,false,  // stressSource disabled, stressField enabled,
                                                               true,false,   //   dispSource  enabled,   dispField enabled,
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
                                                           this->source->sID,this->sink->sID,this->quadraturePoints().size()+qOrder12+k,
                                                           d24*L24,
                                                           Burgers,
                                                           QuadratureDynamicType::abscissa(qOrder24,k),
                                                           QuadratureDynamicType::weight(qOrder24,k),
                                                           true,false,  // stressSource disabled, stressField enabled,
                                                           true,false,   //   dispSource  enabled,   dispField enabled,
                                                           true,false);
                            }
                            
                            // Place Quadrature-particles on P4->P3
                            if(!this->sink->isPureBoundaryNode())
                            {
                                const VectorDim d43=(P3-P4).normalized();
                                for (unsigned int k=0;k<qOrder12;++k)
                                {
                                    
                                    particleSystem.addParticle(P4+d43*virtualSegmentDistance*QuadratureDynamicType::abscissa(qOrder12,k),
                                                               this->source->sID,this->sink->sID,this->quadraturePoints().size()+qOrder12+qOrder24+k,
                                                               d43*virtualSegmentDistance,
                                                               Burgers,
                                                               QuadratureDynamicType::abscissa(qOrder12,k),
                                                               QuadratureDynamicType::weight(qOrder12,k),
                                                               true,false,  // stressSource disabled, stressField enabled,
                                                               true,false,   //   dispSource  enabled,   dispField enabled,
                                                               true,false);
                                }
                                
                            }
                            
                        }
                        else // bnd segment, with bvp, without virtual segments
                        {
                            for (unsigned int k=0;k<this->quadraturePoints().size();++k)
                            {
                                quadratureParticleContainer.push_back(particleSystem.addParticle(this->quadraturePoint(k).r,
                                                                                                 this->source->sID,this->sink->sID,k,
                                                                                                 this->quadraturePoint(k).ru,Burgers,
                                                                                                 QuadratureDynamicType::abscissa(this->quadraturePoints().size(),k),
                                                                                                 QuadratureDynamicType::weight(this->quadraturePoints().size(),k),
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
            }
        }
        
        /**********************************************************************/
        void assemble(const std::deque<StressStraight<dim>,Eigen::aligned_allocator<StressStraight<dim>>>& straightSegmentsDeq)
        {
            this->updateForcesAndVelocities(*this,quadPerLength,straightSegmentsDeq);
            quadratureParticleContainer.clear();
            quadratureParticleContainer.reserve(this->quadraturePoints().size());
            
            Fq= this->quadraturePoints().size()? this->nodalVelocityVector() : VectorNdof::Zero();
            Kqq=this->nodalVelocityMatrix(*this);
            h2posMap=this->hermite2posMap();
            
            Mseg.setZero(Ncoeff*dim,h2posMap.size()*dim);
            size_t c=0;
            for(const auto& pair : h2posMap)
            {
                for(int r=0;r<Ncoeff;++r)
                {
                    Mseg.template block<dim,dim>(r*dim,c*dim)=pair.second.first(r)*MatrixDim::Identity();
                }
                c++;
            }
        }
        
        
        /**********************************************************************/
        void assemble()
        {/*!Computes the following:
          * - edge stiffness matrix Kqq
          *	\f[
          *		\mathbf{K} = int_0^1 \mathbf{K}^*(u) du
          *	\f]
          * - nodal force vector Fq
          *	\f[
          *		\mathbf{K} = int_0^1 \mathbf{K}^*(u) du
          *	\f]
          * - edge-to-component matrix Mseg
          */
            
            Fq.setZero();
            Kqq.setZero();
            
            if(!hasZeroBurgers())
            {
                assert(0 && "THIS FUNCTION MUST BE REWORKED");
                
                //                //! 1- Compute and store stress and PK-force at quadrature points
                //                stressGauss.clear();
                //                if(!this->network().use_bvp && isBoundarySegment())
                //                {
                //                    pkGauss.setZero(dim,this->quadraturePoints().size());
                //                }
                //                else
                //                {
                //                    for (unsigned int k=0;k<this->quadraturePoints().size();++k)
                //                    {
                //                        MatrixDim temp(quadratureParticleContainer[k]->stress());
                //
                //                        if(this->network().use_externalStress)
                //                        {
                //                            temp+=this->network().extStressController.externalStress(this->quadraturePoint(k).r);
                //                        }
                //
                //                        if(this->network().use_bvp)
                //                        {
                //                            temp += this->network().bvpSolver.stress(this->quadraturePoint(k).r,this->source->includingSimplex());
                //                        }
                //
                //                        for(const auto& sStraight : this->network().poly.grainBoundaryDislocations() )
                //                        {
                //                            temp+=sStraight.stress(this->quadraturePoint(k).r);
                //                        }
                //
                //                        stressGauss.push_back(temp);
                //                        pkGauss.col(k)=(stressGauss[k]*Burgers).cross(this->quadraturePoint(k).rl);
                //                        vGauss.col(k)=getGlideVelocity(k);
                //                    }
                //                }
                //
                //
                //                QuadratureDynamicType::integrate(this->quadraturePoints().size(),this,Fq,&LinkType::velocityIntegrand);
                //
                //
                //                if(corder==0)
                //                {
                //                    const double L=this->chord().norm();
                //                    Kqq<<L/3.0,    0.0,    0.0, L/6.0,    0.0,    0.0,
                //                    /**/ 0.0,    L/3.0,    0.0,     0.0,L/6.0,    0.0,
                //                    /**/ 0.0,    0.0,L/3.0,     0.0,    0.0,    L/6.0,
                //                    /**/ L/6.0,    0.0,    0.0, L/3.0,    0.0,    0.0,
                //                    /**/ 0.0,L/6.0,    0.0,     0.0,    L/3.0,    0.0,
                //                    /**/ 0.0,    0.0,L/6.0,     0.0,    0.0,    L/3.0;
                ////                    Kqq*=this->chord().norm();
                //                }
                //                else
                //                {
                //                    assert(0 && "WE NEED TO INCREASE qOrder FOR THE FOLLOWING INTEGRATION, SINCE EVEN FOR LINEAR SEGMENTS Kqq IS NOT INTEGRATED CORRECLTY FOR SMALL qOrder");
                //                    QuadratureDynamicType::integrate(this->quadraturePoints().size(),this,Kqq,&LinkType::stiffness_integrand);
                //                }
                //
                //                h2posMap=this->hermite2posMap();
                //
                //                Mseg.setZero(Ncoeff*dim,h2posMap.size()*dim);
                //
                //                size_t c=0;
                //                for(const auto& pair : h2posMap)
                //                {
                //                    for(int r=0;r<Ncoeff;++r)
                //                    {
                //                        Mseg.template block<dim,dim>(r*dim,c*dim)=pair.second.first(r)*MatrixDim::Identity();
                //                    }
                //                    c++;
                //                }
                
            }
            
        }
        
        /**********************************************************************/
        void addToGlobalAssembly(std::deque<Eigen::Triplet<double> >& kqqT,
                                 Eigen::VectorXd& FQ) const
        {/*!\param[in] kqqT the stiffness matrix of the network component
          * \param[in] FQ the force vector of the network component
          */
            
            if(!hasZeroBurgers())
            {
                const Eigen::MatrixXd tempKqq(Mseg.transpose()*Kqq*Mseg); // Create the temporaty stiffness matrix and push into triplets
                
                
                size_t localI=0;
                for(const auto& pairI : h2posMap)
                {
                    for(int dI=0;dI<dim;++dI)
                    {
                        const size_t globalI=pairI.first*dim+dI;
                        
                        size_t localJ=0;
                        for(const auto& pairJ : h2posMap)
                        {
                            for(int dJ=0;dJ<dim;++dJ)
                            {
                                const size_t globalJ=pairJ.first*dim+dJ;
                                
                                if (std::fabs(tempKqq(localI,localJ))>FLT_EPSILON)
                                {
                                    //                                    kqqT.push_back(Eigen::Triplet<double>(globalI,globalJ,tempKqq(localI,localJ)));
                                    kqqT.emplace_back(globalI,globalJ,tempKqq(localI,localJ));
                                }
                                
                                localJ++;
                            }
                        }
                        
                        localI++;
                    }
                }
                
                const Eigen::VectorXd tempFq(Mseg.transpose()*Fq); // Create temporary force vector and add to global FQ
                
                localI=0;
                for(const auto& pairI : h2posMap)
                {
                    for(int dI=0;dI<dim;++dI)
                    {
                        const size_t globalI=pairI.first*dim+dI;
                        
                        FQ(globalI)+=tempFq(localI);
                        
                        localI++;
                    }
                }
                
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
            
            //\todo this integral should be calculated using shape functions
            
            const VectorDim V((this->source->get_V().template segment<dim>(0)+this->sink->get_V().template segment<dim>(0))*0.5);
            return -Burgers*V.cross(this->chord()).transpose()*(!isBoundarySegment())/this->network().mesh.volume();
        }
        
        /**********************************************************************/
        MatrixDim plasticStrainRate() const
        {
            const VectorDim temp(plasticDistortionRate());
            return (temp+temp.transpose())*0.5;
        }
        
        /**********************************************************************/
        bool isGrainBoundarySegment() const
        {
            return grainBoundaries().size();
        }
        
        /**********************************************************************/
        Eigen::Matrix<double,dim-1,Ncoeff> hermiteLocalCoefficient() const
        {
            const MatrixDim G2L(DislocationLocalReference<dim>::global2local(this->chord(),glidePlaneNormal()));
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
        void addToSolidAngleJump(const VectorDim& Pf, const VectorDim& Sf, VectorDim& dispJump) const
        {
            if(isBoundarySegment() && this->network().use_virtualSegments)
            {
                // first triangle is P1->P2->P3, second triangle is P2->P4->P3
                const VectorDim P1(this->source->get_P());
                const VectorDim P2(P1+this->source->bndNormal()*virtualSegmentDistance);
                const VectorDim P3(this->sink->get_P());
                const VectorDim P4(P3+this->sink->bndNormal()*virtualSegmentDistance);
                
                dispJump += Burgers*LineSimplexIntersection<dim>::lineTriangleIntersection(Pf,Sf,P1,P2,P3);
                dispJump += Burgers*LineSimplexIntersection<dim>::lineTriangleIntersection(Pf,Sf,P2,P4,P3);
            }
        }
        
        /**********************************************************************/
        const VectorDim& burgers() const
        {
            return Burgers;
        }
        
        /**********************************************************************/
        const VectorDim& glidePlaneNormal() const
        {
            return meshPlanes().size()==1? (*meshPlanes().begin())->unitNormal : zeroVector;
        }
        
        /**********************************************************************/
        bool isSessile() const
        {
            return    !isGlissile()
            /*  */ && !isBoundarySegment()
            /*  */ && !isGrainBoundarySegment()
            /*  */ && !hasZeroBurgers();
        }
        
        /**********************************************************************/
        bool isGlissile() const
        {/*\returns true if ALL the following conditions are met
          * - the segment is confined by only one plane
          * - its Burgers vector is non-zero
          * - all loops containing this segment are glissile
          */
            bool temp(meshPlanes().size()==1 && !hasZeroBurgers());
            if(temp)
            {
                for(const auto& loopLink : this->loopLinks())
                {
                    temp*=loopLink->loop()->isGlissile;
                }
            }
            return temp;
        }
        
        /**********************************************************************/
        bool hasZeroBurgers() const
        {
            return Burgers.squaredNorm()<FLT_EPSILON;
        }
        
        //        /**********************************************************************/
        //        double arcLength() const
        //        {
        //            return SplineSegmentType::template arcLength<16,UniformOpen>();
        //        }
        
        //        /**********************************************************************/
        //        VectorDim velocity(const double& u) const
        //        {
        //            return this->source->get_V().template segment<dim>(0)*(1.0-u)+this->sink->get_V().template segment<dim>(0)*u;
        //        }
        
        /**********************************************************************/
        const MatrixDim& midPointStress() const
        {/*!\returns The stress matrix for the centre point over this segment.*/
            //            return stressGauss[this->quadraturePoints().size()/2];
            
            return this->quadraturePoints().size()? quadraturePoint(this->quadraturePoints().size()/2).stress : MatrixDim::Zero();
            
        }
        
        /**********************************************************************/
        bool isBoundarySegment() const // THIS IS CALLED MANY TIMES< CONSIDER STORING
        {/*!\returns true if both nodes are boundary nodes, and the midpoint is
          * on the boundary.
          */
            return this->source->isBoundaryNode() &&
            /*  */ this->sink->isBoundaryNode() &&
            /*  */ boundingBoxSegments().contains(0.5*(this->source->get_P()+this->sink->get_P())).first;
        }
        
        /**********************************************************************/
        MeshLocation meshLocation() const
        {/*!\returns the position of *this relative to the bonudary:
          * 1 = inside mesh
          * 2 = on mesh boundary
          */
            
            MeshLocation temp = MeshLocation::outsideMesh;
            
            
            if(isBoundarySegment())
            {
                temp=MeshLocation::onMeshBoundary;
            }
            else
            {
                if(isGrainBoundarySegment())
                {
                    temp=MeshLocation::onRegionBoundary;
                }
                else
                {
                    temp=MeshLocation::insideMesh;
                }
            }
            
            return temp;
        }
        
        /**********************************************************************/
        template <class T>
        friend T& operator << (T& os, const LinkType& ds)
        {
            os  << ds.source->sID<<"\t"<< ds.sink->sID<<"\t"
            /**/<< std::setprecision(15)<<std::scientific<<ds.Burgers.transpose()<<"\t"
            /**/<< std::setprecision(15)<<std::scientific<<ds.glidePlaneNormal().transpose()<<"\t"
            /**/<<SplineSegmentBase<dim,corder>::sourceT(ds).transpose()<<"\t"
            /**/<<SplineSegmentBase<dim,corder>::sinkT(ds).transpose()<<"\t"
            /**/<<ds.meshLocation();
            return os;
        }
        
    };
    
    // Static Data
    template <int dim, short unsigned int corder, typename InterpolationType>
    const Eigen::Matrix<double,dim,dim> DislocationSegment<dim,corder,InterpolationType>::I=Eigen::Matrix<double,dim,dim>::Identity();
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    const Eigen::Matrix<double,dim,1> DislocationSegment<dim,corder,InterpolationType>::zeroVector=Eigen::Matrix<double,dim,1>::Zero();
    
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    double DislocationSegment<dim,corder,InterpolationType>::quadPerLength=0.2;
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    double DislocationSegment<dim,corder,InterpolationType>::virtualSegmentDistance=200.0;
    
}
#endif

