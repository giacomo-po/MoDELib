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


#ifndef model_PlanarDislocationSegment_H
#define model_PlanarDislocationSegment_H

#include <memory>
#include <set>
#include <algorithm>    // std::set_intersection, std::sort


#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Sparse>

#include <Quadrature.h>
#include <QuadPow.h>
#include <DislocationNetworkTraits.h>
#include <SplineSegment.h>
//#include <Material.h>
#include <GlidePlane.h>
#include <Coeff2Hermite.h>
#include <DislocationParticle.h>
#include <UniqueOutputFile.h>
#include <GrainBoundary.h>
#include <LineSimplexIntersection.h>
//#include <BoundingLineSegments.h>
//#include <BoundingLineSegments.h>
#include <MeshPlane.h>
#include <DislocationQuadraturePoint.h>
#include <StraightDislocationSegment.h>
#include <ConfinedDislocationObject.h>
#include <DislocationLoopIO.h>

#ifndef NDEBUG
#define VerbosePlanarDislocationSegment(N,x) if(verbosePlanarDislocationSegment>=N){model::cout<<x;}
#else
#define VerbosePlanarDislocationSegment(N,x)
#endif


namespace model
{
    
    /**************************************************************************/
    /**************************************************************************/
    template <typename Derived>
    class PlanarDislocationSegment : public SplineSegment<Derived,TypeTraits<Derived>::dim,TypeTraits<Derived>::corder>
    /*                            */,public ConfinedDislocationObject<TypeTraits<Derived>::dim>
    //    /*                            */,private std::set<const GlidePlane<TypeTraits<Derived>::dim>*>
    ////    /*                      */,private std::set<const GrainBoundary<TypeTraits<Derived>::dim>*>
    //    /*                            */,private std::set<const Grain<TypeTraits<Derived>::dim>*>
    //    /*                            */,private std::set<const PlanarMeshFace<TypeTraits<Derived>::dim>*>
    //    /*                            */,private BoundingLineSegments<TypeTraits<Derived>::dim>
    /*                            */,public DislocationQuadraturePointContainer<TypeTraits<Derived>::dim,TypeTraits<Derived>::corder>
    {
        
    public:
        
        static constexpr int dim=TypeTraits<Derived>::dim; // make dim available outside class
        static constexpr int corder=TypeTraits<Derived>::corder; // make dim available outside class
        typedef SplineSegmentBase<dim,corder> SplineSegmentBaseType;
        typedef ConfinedDislocationObject<dim> ConfinedDislocationObjectType;
        typedef Derived LinkType;
        typedef StraightDislocationSegment<dim> StraightDislocationSegmentType;
        typedef typename TypeTraits<LinkType>::LoopType LoopType;
        typedef typename TypeTraits<LinkType>::LoopNetworkType NetworkType;
        typedef LoopLink<LinkType> LoopLinkType;
        typedef typename TypeTraits<LinkType>::NodeType NodeType;
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
        typedef DislocationQuadraturePoint<dim,corder> DislocationQuadraturePointType;
        typedef DislocationParticle<dim> DislocationParticleType;
        typedef LatticeVector<dim> LatticeVectorType;
        typedef ReciprocalLatticeDirection<dim> ReciprocalLatticeDirectionType;
        //        typedef std::set<const GrainBoundary<dim>*> GrainBoundaryContainerType;
        //        typedef std::set<const Grain<dim>*> GrainContainerType;
        //        typedef GlidePlane<dim> GlidePlaneType;
        //        typedef MeshPlane<dim> MeshPlaneType;
        //        typedef std::set<const GlidePlaneType*> GlidePlaneContainerType;
        typedef typename TypeTraits<LinkType>::MeshLocation MeshLocation;
        //        typename std::set<const PlanarMeshFace<TypeTraits<Derived>::dim>*> PlanarMeshFaceContainerType;
        
    private:
        
        std::map<size_t,
        /*    */ std::pair<VectorNcoeff,VectorDim>,
        /*    */ std::less<size_t>
        /*    */ > h2posMap;
        Eigen::Matrix<double, Ndof, Eigen::Dynamic> Mseg;
        MatrixNdof Kqq; //! Segment Stiffness Matrix
        VectorNdof Fq; //! Segment Nodal Force Vector
        VectorDim Burgers; //! The Burgers vector
        double BurgersNorm;
        //        bool _isBoundarySegment;
        //        bool _isGrainBoundarySegment;
        
    public:
        
        static const Eigen::Matrix<double,TypeTraits<Derived>::dim,TypeTraits<Derived>::dim> I;
        static const Eigen::Matrix<double,TypeTraits<Derived>::dim,1> zeroVector;
        static double quadPerLength;
        static int verbosePlanarDislocationSegment;
        StraightDislocationSegment<dim> straight;
        
    public:
        
        /******************************************************************/
        static void initFromFile(const std::string& fileName)
        {
            LinkType::alpha=TextFileParser(fileName).readScalar<double>("parametrizationExponent",true);
            assert((LinkType::alpha)>=0.0 && "parametrizationExponent MUST BE >= 0.0");
            assert((LinkType::alpha)<=1.0 && "parametrizationExponent MUST BE <= 1.0");
            quadPerLength=TextFileParser(fileName).readScalar<double>("quadPerLength",true);
            assert((LinkType::quadPerLength)>=0.0 && "quadPerLength MUST BE >= 0.0");
            verbosePlanarDislocationSegment=TextFileParser("inputFiles/DD.txt").readScalar<int>("verbosePlanarDislocationSegment",true);
        }
        
        /******************************************************************/
        PlanarDislocationSegment(const std::shared_ptr<NodeType>& nI,
                                 const std::shared_ptr<NodeType>& nJ) :
//        /* init */ ConfinedDislocationObjectType(nI->network())
        /* init */ SplineSegmentType(nI,nJ)
        /* init */,ConfinedDislocationObjectType(this->source->get_P(),this->sink->get_P())
        /* init */,Burgers(VectorDim::Zero())
        /* init */,BurgersNorm(Burgers.norm())
        //        /* init */,_isBoundarySegment(false)
        //        /* init */,_isGrainBoundarySegment(false)
        /* init */,straight(this->source->get_P(),this->sink->get_P(),Burgers,this->chordLength(),this->unitDirection())
        {/*! Constructor with pointers to source and sink, and flow
          *  @param[in] NodePair_in the pair of source and sink pointers
          *  @param[in] Flow_in the input flow
          */
            VerbosePlanarDislocationSegment(1,"Creating PlanarDislocationSegment "<<this->tag()<<std::endl;);
            //            VerbosePlanarDislocationSegment(3,"source->isBoundaryNode() "<<this->source->isBoundaryNode()<<std::endl;);
            //            VerbosePlanarDislocationSegment(3,"sink->isBoundaryNode() "<<this->sink->isBoundaryNode()<<std::endl;);
            //            VerbosePlanarDislocationSegment(2,"_isBoundarySegment "<<_isBoundarySegment<<std::endl;);
            //            VerbosePlanarDislocationSegment(3,"midpoint is boundary "<<boundingBoxSegments().contains(0.5*(this->source->get_P()+this->sink->get_P()))<<std::endl;);
        }
        
        /**********************************************************************/
        ~PlanarDislocationSegment()
        {
            VerbosePlanarDislocationSegment(1,"Destroying PlanarDislocationSegment "<<this->tag()<<std::endl;);
        }
        
        /**********************************************************************/
        void updateGeometry()
        {
            SplineSegmentType::updateGeometry();
            ConfinedDislocationObjectType::updateGeometry(this->source->get_P(),this->sink->get_P());
            straight.updateGeometry();
            //            addMeshFaces();
            //            _isBoundarySegment=this->source->isBoundaryNode() && this->sink->isBoundaryNode() && boundingBoxSegments().contains(0.5*(this->source->get_P()+this->sink->get_P()));
            //            _isBoundarySegment=boundingBoxSegments().boundaryNormal(0.5*(this->source->get_P()+this->sink->get_P())).squaredNorm()>FLT_EPSILON;
        }
        
        
        
        /**********************************************************************/
        void addLoopLink(LoopLinkType* const pL)
        {
            VerbosePlanarDislocationSegment(2,"PlanarDislocationSegment "<<this->tag()<<", adding LoopLink "<<pL->tag()<<std::endl;);
            SplineSegmentType::addLoopLink(pL); // forward to base class
            if(pL->source()->sID==this->source->sID)
            {// Update Burgers vector
                Burgers+=pL->flow().cartesian();
            }
            else
            {
                Burgers-=pL->flow().cartesian();
            }
            BurgersNorm=Burgers.norm();
            straight.updateGeometry();
            
            VerbosePlanarDislocationSegment(3,"adding GlidePlane with bounding box:\n"<<pL->loop()->glidePlane->meshIntersections<<std::endl;);
            ConfinedDislocationObjectType::addGlidePlane(pL->loop()->glidePlane.get());
        }
        

        
        /**********************************************************************/
        void removeLoopLink(LoopLinkType* const pL)
        {
            VerbosePlanarDislocationSegment(2,"PlanarDislocationSegment "<<this->tag()<<", removing LoopLink "<<pL->tag()<<std::endl;);
            SplineSegmentType::removeLoopLink(pL);  // forward to base class
            
            if(pL->source()->sID==this->source->sID)
            {// Modify Burgers vector
                Burgers-=pL->flow().cartesian();
            }
            else
            {
                Burgers+=pL->flow().cartesian();
            }
            BurgersNorm=Burgers.norm();
            straight.updateGeometry(); // update b x t
            
            ConfinedDislocationObjectType::clear();
            for(const auto& loopLink : this->loopLinks())
            {
                ConfinedDislocationObjectType::addGlidePlane(loopLink->loop()->glidePlane.get());
            }
            
            
            //            if(!pL->loop()->isVirtualBoundaryLoop())
            //            {
            //                glidePlanes().clear();
            //                grains().clear();
            //                grainBoundaries.clear();
            //                boundingBoxSegments().clear();
            //                for(const auto& loopLink : this->loopLinks())
            //                {
            //                    if(loopLink->loop()->glidePlane)
            //                    {
            //                        addMeshPlane(*loopLink->loop()->glidePlane.get());
            //                    }
            //                }
            ////                if(glidePlanes().size())
            ////                {
            ////                    addMeshFaces();
            ////                }
            //                addMeshFaces();
            //
            ////                _isBoundarySegment=this->source->isBoundaryNode() && this->sink->isBoundaryNode() && boundingBoxSegments().contains(0.5*(this->source->get_P()+this->sink->get_P()));
            ////                _isBoundarySegment=boundingBoxSegments().boundaryNormal(0.5*(this->source->get_P()+this->sink->get_P())).squaredNorm()>FLT_EPSILON;
            //            }
        }
        

        
        /**********************************************************************/
        void assembleGlide()
        {
            this->updateForcesAndVelocities(*this,quadPerLength,false);
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
        const VectorDim& burgers() const
        {
            return Burgers;
        }
        
        /**********************************************************************/
        const VectorDim& glidePlaneNormal() const
        {
            return this->glidePlanes().size()==1? (*this->glidePlanes().begin())->unitNormal : zeroVector;
        }
        
        /**********************************************************************/
        std::set<const LoopType*> virtualLoops() const
        {//!\returns a set of pointers to the virtualBoundaryLoops of this segment
            std::set<const LoopType*> temp;
            for(const auto& loopLink : this->loopLinks())
            {
                if(loopLink->loop()->isVirtualBoundaryLoop())
                {
                    temp.insert(loopLink->loop().get());
                }
            }
            return temp;
        }
        
        
        bool hasPeriodicLoop() const
        {
            bool temp(false);
            for(const auto& loopLink : this->loopLinks())
            {
                if(loopLink->loop()->loopType==DislocationLoopIO<dim>::PERIODICLOOP)
                {
                    temp=true;
                }
            }
            return temp;
        }
        
        /**********************************************************************/
        bool hasZeroBurgers() const
        {
            return BurgersNorm<FLT_EPSILON;
        }
        
        /**********************************************************************/
        bool isBoundarySegment() const 
        {/*!\returns true if both nodes are boundary nodes, and the midpoint is
          * on the boundary.
          */
            return this->isOnExternalBoundary();
        }
        
        /**********************************************************************/
        bool isGrainBoundarySegment() const
        {
            return this->isOnInternalBoundary();
        }
        
        /**********************************************************************/
        bool isVirtualBoundarySegment() const
        {//!\returns true if all loops of this segment are virtualBoundaryLoops
            bool temp(true);
            for(const auto& loopLink : this->loopLinks())
            {
                temp*=loopLink->loop()->isVirtualBoundaryLoop();
                if(!temp)
                {
                    break;
                }
            }
            return temp;
        }
        
        /**********************************************************************/
        bool isGlissile() const
        {/*\returns true if ALL the following conditions are met
          * - the segment is confined by only one plane
          * - its Burgers vector is non-zero
          * - all loops containing this segment are glissile
          */
//            bool temp(this->glidePlanes().size()==1 && !hasZeroBurgers() && !isVirtualBoundarySegment());
            bool temp(this->glidePlanes().size()==1 && !isVirtualBoundarySegment());
            if(temp)
            {
                for(const auto& loopLink : this->loopLinks())
                {
                    temp*=loopLink->loop()->loopType==DislocationLoopIO<dim>::GLISSILELOOP;
                }
            }
            return temp;
        }
        
        /**********************************************************************/
        bool isSessile() const
        {
            return !isVirtualBoundarySegment() && !isGlissile();
        }
        
        /**********************************************************************/
        bool isPureSessile() const
        {
            bool temp(true);
            for(const auto& loopLink : this->loopLinks())
            {
                temp*=loopLink->loop()->loopType==DislocationLoopIO<dim>::SESSILELOOP;
            }
            return temp;
        }
        
//        /**********************************************************************/
//        std::set<size_t> faceIDs() const
//        {
//            std::set<size_t> key;
//            for(const auto& face : this->meshFaces())
//            {
//                key.insert(face->sID);
//            }
//            return key;
//        }
//
//        /**********************************************************************/
//        std::set<size_t> imageFaceIDs() const
//        {
//            assert(this->network().mesh.regions().size()==1);
//            const auto& region(*this->network().mesh.regions().begin()->second);
//
//            std::set<size_t> key;
//            for(const auto& face : this->meshFaces())
//            {
//                key.insert(region.parallelFaces().at(face->sID));
//            }
//            return key;
//        }
//
//        /**********************************************************************/
//        std::set<const PlanarMeshFace<dim>*> imageFaces() const
//        {
//            assert(this->network().mesh.regions().size()==1);
//            const auto& region(*this->network().mesh.regions().begin()->second);
//
//            std::set<const PlanarMeshFace<dim>*> temp;
//            for(const auto& face : this->meshFaces())
//            {
//                temp.insert(region.faces().at(region.parallelFaces().at(face->sID)).get());
//            }
//            return temp;
//        }
        
        /**********************************************************************/
        MeshLocation meshLocation() const
        {/*!\returns the position of *this relative to the bonudary:
          * 1 = inside mesh
          * 2 = on mesh boundary
          */
            if(isBoundarySegment())
            {
                return MeshLocation::onMeshBoundary;
            }
            else if(isGrainBoundarySegment())
            {
                return MeshLocation::onRegionBoundary;
            }
            else if(isVirtualBoundarySegment())
            {
                return MeshLocation::outsideMesh;
            }
            else
            {
                return MeshLocation::insideMesh;
            }
        }
    };
    
    // Static Data
    template <typename Derived>
    const Eigen::Matrix<double,TypeTraits<Derived>::dim,TypeTraits<Derived>::dim> PlanarDislocationSegment<Derived>::I=Eigen::Matrix<double,TypeTraits<Derived>::dim,TypeTraits<Derived>::dim>::Identity();
    
    template <typename Derived>
    const Eigen::Matrix<double,TypeTraits<Derived>::dim,1> PlanarDislocationSegment<Derived>::zeroVector=Eigen::Matrix<double,TypeTraits<Derived>::dim,1>::Zero();
    
    template <typename Derived>
    double PlanarDislocationSegment<Derived>::quadPerLength=0.2;
    
    template <typename Derived>
    int PlanarDislocationSegment<Derived>::verbosePlanarDislocationSegment=0;
}
#endif

