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

#ifndef model_DislocationSegment_cpp_
#define model_DislocationSegment_cpp_

#include <DislocationSegment.h>



namespace model
{
    template <int dim, short unsigned int corder>
    DislocationSegment<dim,corder>::DislocationSegment(LoopNetworkType* const net,
                                                                         const std::shared_ptr<NetworkNodeType>& nI,
                                                                         const std::shared_ptr<NetworkNodeType>& nJ) :
    /* init */ NetworkLink<DislocationSegment>(net,nI,nJ)
    /* init */,SplineSegment<dim,corder>(this->source->get_P(),this->sink->get_P())
    // /* init */,ConfinedDislocationObjectType(this->source->get_P(),this->sink->get_P())
    /* init */,Burgers(VectorDim::Zero())
    /* init */,BurgersNorm(Burgers.norm())
    /* init */,straight(this->source->get_P(),this->sink->get_P(),Burgers,this->chordLength(),this->unitDirection())
    /* init */,_slipSystem(nullptr)
    {
        VerboseDislocationSegment(1,"Constructing DislocationSegment "<<this->tag()<<std::endl);
    }
    
    // template <int dim, short unsigned int corder>
    // void DislocationSegment<dim,corder>::updateSlipSystem()
    // {
    //     std::set<std::shared_ptr<SlipSystem>> ssSet;
    //     for(const auto& loopLink : this->loopLinks())
    //     {
    //         if(loopLink->loop->slipSystem())
    //         {
    //             ssSet.insert(loopLink->loop->slipSystem());
    //         }
    //     }
        
    //     if(ssSet.size()==1)
    //     {// a unique slip system found. TO DO. This fails for planar glissile junctions, since the two loop links have different slip systems but the resultant is glissile on a third slip system.
    //         _slipSystem=*ssSet.begin();
    //     }
    //     else
    //     {
    //         _slipSystem=nullptr;
    //     }
    //     if(_slipSystem)
    //     {
    //         VerboseDislocationSegment(3,"_slipSystem= "<<_slipSystem->s.cartesian().transpose()<<std::endl;);
    //         VerboseDislocationSegment(3,"_slipSystem= "<<_slipSystem->unitNormal.transpose()<<std::endl;);
    //     }
    // }

    // template <int dim, short unsigned int corder>
    // void DislocationSegment<dim, corder>::updateSlipSystem()
    // {
    //     if (this->grains().size() == 0)
    //     {
    //          _slipSystem=nullptr;
    //     }
    //     else if (this->grains().size() == 1)
    //     {
    //         std::set<std::shared_ptr<SlipSystem>> ssSet;

    //         for (const auto &loopLink : this->loopLinks())
    //         {
    //             if (loopLink->loop->slipSystem())
    //             {
    //                 ssSet.insert(loopLink->loop->slipSystem());
    //             }
    //         }

    //         if (ssSet.size() == 1)
    //         { // a unique slip system found. TO DO. This fails for planar glissile junctions, since the two loop links have different slip systems but the resultant is glissile on a third slip system.
    //             _slipSystem = *ssSet.begin();
    //         }
    //         else if (ssSet.size() > 1)
    //         {
    //             const auto firstSlipSystem(*ssSet.begin());
    //             //            const auto& firstN(firstSlipSystem.n);
    //             std::shared_ptr<RationalLatticeDirection<dim>> s(new RationalLatticeDirection<dim> (firstSlipSystem->s * 0));
    //             for (const auto &ss : ssSet)
    //             {
    //                 _slipSystem = nullptr;
    //                 if (ss->n.cross(firstSlipSystem->n).squaredNorm() == 0)
    //                 {
    //                     if (ss->n.cartesian().dot(firstSlipSystem->n.cartesian()) > 0.0)
    //                     { // aligned normals
    //                         // s = s + ss->s;
    //                         s.reset(new RationalLatticeDirection<dim>(*s + ss->s));
    //                     }
    //                     else
    //                     { // opposite normals
    //                         // s = s - ss->s;
    //                         s.reset(new RationalLatticeDirection<dim>(*s - ss->s));

    //                     }
    //                 }
    //                 else
    //                 {
    //                     s.reset(new RationalLatticeDirection<dim>(firstSlipSystem->s * 0));
    //                     // s = firstSlipSystem->s * 0;
    //                     break;
    //                 }
    //             }
    //             for (const auto &ss : (*this->grains().begin())->slipSystems())
    //             {
    //                 if (ss->isSameAs(*s, firstSlipSystem->n))
    //                 {
    //                     _slipSystem = ss;
    //                     break;
    //                 }
    //             }
    //         }
    //         else
    //         { // ssSet.size()==0
    //             _slipSystem = nullptr;
    //         }
    //     }
    //     else
    //     {
    //         assert(false && "FINISH THIS FOR MULTIPLE GRAINS");
    //         _slipSystem = nullptr;
    //     }

    //     if (_slipSystem)
    //     {
    //         VerboseDislocationSegment(3, "_slipSystem= " << _slipSystem->s.cartesian().transpose() << std::endl;);
    //         VerboseDislocationSegment(3, "_slipSystem= " << _slipSystem->unitNormal.transpose() << std::endl;);
    //     }
    // }
       template <int dim, short unsigned int corder>
    void DislocationSegment<dim, corder>::updateSlipSystem()
    {
        if (this->grains().size() == 0)
        {
             _slipSystem=nullptr;
        }
        else if (this->grains().size() == 1)
        {
            std::set<std::shared_ptr<SlipSystem>> ssSet;

            for (const auto &loopLink : this->loopLinks())
            {
                if (loopLink->loop->slipSystem())
                {
                    ssSet.insert(loopLink->loop->slipSystem());
                }
            }

            if (ssSet.size() == 1)
            { // a unique slip system found. TO DO. This fails for planar glissile junctions, since the two loop links have different slip systems but the resultant is glissile on a third slip system.
                _slipSystem = *ssSet.begin();
            }
            else
            { // ssSet.size()==0
                _slipSystem = nullptr;
            }
        }
        else
        {
            assert(false && "FINISH THIS FOR MULTIPLE GRAINS");
            _slipSystem = nullptr;
        }

        if (_slipSystem)
        {
            VerboseDislocationSegment(3, "_slipSystem= " << _slipSystem->s.cartesian().transpose() << std::endl;);
            VerboseDislocationSegment(3, "_slipSystem= " << _slipSystem->unitNormal.transpose() << std::endl;);
        }
    }

    /**********************************************************************/
    template <int dim, short unsigned int corder>
    void DislocationSegment<dim,corder>::addLoopLink(LoopLinkType* const pL)
    {
        VerboseDislocationSegment(2,"DislocationSegment "<<this->tag()<<", adding LoopLink "<<pL->tag()<<std::endl;);
        NetworkLink<DislocationSegment>::addLoopLink(pL); // forward to base class
        if(pL->source->networkNode==this->source)
        {// Update Burgers vector
            Burgers+=pL->flow().cartesian();
        }
        else
        {
            Burgers-=pL->flow().cartesian();
        }
        BurgersNorm=Burgers.norm();
        straight.updateGeometry();
        
        //assert(0 && "Following addGlidePlane is wrong, must consider shifts");
        
        if(this->network().simulationParameters.isPeriodicSimulation())
        {
            const auto periodicPlanePatch(pL->periodicPlanePatch());
            if(periodicPlanePatch)
            {
                // this->confinedObject().addGlidePlane(periodicPlanePatch->glidePlane.get());
               // /* Check for the consistency of the glide planes
                //assert if the confined object glide planes are the same as that of the each node glide planes
                bool glidePlanesConsistency(false);
                for (const auto &Seggps : glidePlanes())
                {
                    glidePlanesConsistency = false;
                    for (const auto &sourceGPs : this->source->glidePlanes())
                    {
                        if (Seggps == sourceGPs)
                        {
                            glidePlanesConsistency = true;
                            break;
                        }
                    }
                }
                // std::cout<<" A " <<glidePlanesConsistency<<std::endl;

                if (glidePlanesConsistency)
                {
                    for (const auto &Seggps : glidePlanes())
                    {
                        glidePlanesConsistency=false;
                        for (const auto &sinkGPs : this->sink->glidePlanes())
                        {
                            if (Seggps==sinkGPs)
                            {
                                glidePlanesConsistency=true;
                                break;
                            }
                        }
                    }
                }
                // std::cout<<" B " <<glidePlanesConsistency<<std::endl;

                assert(glidePlanesConsistency && "Not a consistent definition of glide planes for segments and segment nodes");
              // Check for the consistency of glide planes end here  */
            }
        }
        // else
        // {
        //     this->confinedObject().addGlidePlane(pL->loop->glidePlane.get());
        // }


        updateSlipSystem();
    }
    
    
    
    /**********************************************************************/
    template <int dim, short unsigned int corder>
    void DislocationSegment<dim,corder>::removeLoopLink(LoopLinkType* const pL)
    {
       VerboseDislocationSegment(2,"DislocationSegment "<<this->tag()<<", removing LoopLink "<<pL->tag()<<std::endl;);
        NetworkLink<DislocationSegment>::removeLoopLink(pL);  // forward to base class
        
        if(pL->source->networkNode==this->source)
        {// Modify Burgers vector
            Burgers-=pL->flow().cartesian();
        }
        else
        {
            Burgers+=pL->flow().cartesian();
        }
        BurgersNorm=Burgers.norm();
        straight.updateGeometry(); // update b x t
        
        // this->confinedObject().clear();
        // for(const auto& loopLink : this->loopLinks())
        // {
        //     // std::cout<<"loopLink"<<loopLink->tag()<<" Adding glidePlane "<<std::endl;
        //     if(this->network().simulationParameters.isPeriodicSimulation())
        //     {
        //         const auto periodicPlanePatch(loopLink->periodicPlanePatch());
        //         // if(periodicPlanePatch)
        //         // {
        //         //     this->confinedObject().addGlidePlane(periodicPlanePatch->glidePlane.get());
        //         // }
        //     }
        //     else
        //     {
        //         this->confinedObject().addGlidePlane(loopLink->loop->glidePlane.get());
        //     }
        // }
        
        updateSlipSystem();
        
    }
    
    template <int dim, short unsigned int corder>
    void DislocationSegment<dim,corder>::updateGeometry()
    {
        VerboseDislocationSegment(2,"DislocationSegment "<<this->tag()<<", updateGeometry "<<std::endl;);
        //        if(this->network().simulationParameters.isPeriodicSimulation())
        //        {
        //            bool periodicPlanePatchesReady(true);
        //            for(const auto& loopLink : this->loopLinks())
        //            {
        //                VerboseDislocationSegment(3,"loopLink "<<loopLink->tag()<<", (loopLink->periodicPlanePatch()!=nullptr)= "<<(loopLink->periodicPlanePatch()!=nullptr)<<std::endl;);
        //                periodicPlanePatchesReady*=(loopLink->periodicPlanePatch()!=nullptr);
        //            }
        //            if(periodicPlanePatchesReady)
        //            {
        //                SplineSegmentType::updateGeometry();
        //                this->confinedObject().clear();
        //                this->updateConfinement(this->source->get_P(),this->sink->get_P());
        //                for(const auto& loopLink : this->loopLinks())
        //                {
        //                    VerboseDislocationSegment(3,"adding glidePlane for loopLink "<<loopLink->tag()<<std::endl;);
        //                        this->confinedObject().addGlidePlane(loopLink->periodicPlanePatch()->glidePlane.get());
        //                }
        //                straight.updateGeometry();
        //                updateSlipSystem();
        //            }
        //        }
        //        else
        //        {
        //            SplineSegmentType::updateGeometry();
        //            this->updateConfinement(this->source->get_P(),this->sink->get_P());
        //            straight.updateGeometry();
        //        }
        SplineSegmentType::updateGeometry();
        // this->updateConfinement(this->source->get_P(),this->sink->get_P());
        straight.updateGeometry();
        VerboseDislocationSegment(2,"DislocationSegment "<<this->tag()<<", updateGeometry DONE"<<std::endl;);
    }
    
    template <int dim, short unsigned int corder>
    bool DislocationSegment<dim,corder>::isGlissile() const
    {/*\returns true if ALL the following conditions are met
      * - the segment is confined by only one plane
      * - all loops containing this segment are glissile
      */
        //            bool temp(this->glidePlanes().size()==1 && !hasZeroBurgers() && !isVirtualBoundarySegment());
        bool temp(this->glidePlanes().size()==1 && !isVirtualBoundarySegment());
        if(temp)
        {
            for(const auto& loopLink : this->loopLinks())
            {
                temp*=loopLink->loop->loopType==DislocationLoopIO<dim>::GLISSILELOOP;
            }
        }
        return temp;
    }
    
    template <int dim, short unsigned int corder>
    bool DislocationSegment<dim,corder>::isSessile() const
    {
        return !isVirtualBoundarySegment() && !isGlissile();
    }
    
    /******************************************************************/
    template <int dim, short unsigned int corder>
    void DislocationSegment<dim,corder>::initFromFile(const std::string& fileName)
    {
        //        LinkType::alpha=TextFileParser(fileName).readScalar<double>("parametrizationExponent",true);
        //        assert((LinkType::alpha)>=0.0 && "parametrizationExponent MUST BE >= 0.0");
        //        assert((LinkType::alpha)<=1.0 && "parametrizationExponent MUST BE <= 1.0");
        quadPerLength=TextFileParser(fileName).readScalar<double>("quadPerLength",true);
        //            assembleWithTangentProjection=TextFileParser(fileName).readScalar<int>("assembleWithTangentProjection",true);
        assert((NetworkLinkType::quadPerLength)>=0.0 && "quadPerLength MUST BE >= 0.0");
        verboseDislocationSegment=TextFileParser("inputFiles/DD.txt").readScalar<int>("verboseDislocationSegment",true);
    }
    
    template <int dim, short unsigned int corder>
    typename DislocationSegment<dim,corder>::MeshLocation DislocationSegment<dim,corder>::meshLocation() const
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
    
    template <int dim, short unsigned int corder>
    bool DislocationSegment<dim,corder>::isVirtualBoundarySegment() const
    {//!\returns true if all loops of this segment are virtualBoundaryLoops
        bool temp(true);
        for(const auto& loopLink : this->loopLinks())
        {
            temp*=loopLink->loop->isVirtualBoundaryLoop();
            if(!temp)
            {
                break;
            }
        }
        return temp;
    }
    
    template <int dim, short unsigned int corder>
    bool DislocationSegment<dim,corder>::hasZeroBurgers() const
    {
        return BurgersNorm<FLT_EPSILON;
    }
    
    template <int dim, short unsigned int corder>
    const std::shared_ptr<SlipSystem>&  DislocationSegment<dim,corder>::slipSystem() const
    {
        return _slipSystem;
    }
    
    template <int dim, short unsigned int corder>
    bool DislocationSegment<dim,corder>::isBoundarySegment() const
    {/*!\returns true if both nodes are boundary nodes, and the midpoint is
      * on the boundary.
      */
        return this->isOnExternalBoundary();
    }
    
    template <int dim, short unsigned int corder>
    bool DislocationSegment<dim,corder>::isGrainBoundarySegment() const
    {
        return this->isOnInternalBoundary();
    }
    
    template <int dim, short unsigned int corder>
    const typename DislocationSegment<dim,corder>::VectorDim& DislocationSegment<dim,corder>::burgers() const
    {
        return Burgers;
    }
    
    template <int dim, short unsigned int corder>
    void DislocationSegment<dim,corder>::addToGlobalAssembly(std::deque<Eigen::Triplet<double> >& kqqT,
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
    
    template <int dim, short unsigned int corder>
    void DislocationSegment<dim,corder>::updateQuadraturePointsSeg()
    {
        this->updateQuadraturePoints(*this,quadPerLength,false);
    }

    template <int dim, short unsigned int corder>
    void DislocationSegment<dim,corder>::assembleGlide(const bool& computeForcesanVelocities)
    {
        VerboseDislocationSegment(2,"DislocationSegment "<<this->tag()<<", assembleGlide"<<std::endl;);
        if (computeForcesanVelocities)
        {
            this->updateForcesAndVelocities(*this,quadPerLength,false);
        }
        Fq= this->quadraturePoints().size()? this->nodalVelocityVector(*this) : VectorNdof::Zero();
        Kqq=this->nodalVelocityMatrix(*this);
        h2posMap.clear();
        switch (corder)
        {
            case 0:
            {
                h2posMap.emplace(this->source->networkID(),std::make_pair((VectorNcoeff()<<1.0,0.0).finished(),this->source->get_P()));
                h2posMap.emplace(this->  sink->networkID(),std::make_pair((VectorNcoeff()<<0.0,1.0).finished(),this->  sink->get_P()));
                break;
            }
                
            default:
            {
                assert(0 && "IMPLEMENT THIS CASE FOR CURVED SEGMENTS");
                break;
            }
        }
        //        h2posMap=this->hermite2posMap();
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

    template <int dim, short unsigned int corder>
    int DislocationSegment<dim, corder>::velocityGroup(const double &maxVelocity, const std::set<int> &subcyclingBins) const
    // This represents how many number of steps taken for updating the velocity
    {
        if (this->quadraturePoints().size() == 0)
        {
            return *subcyclingBins.begin();
        }
        else
        {
            if (maxVelocity > FLT_EPSILON)
            {
                const double avgV((0.5 * (this->source->get_V() + this->sink->get_V())).norm());

                if (avgV >= FLT_EPSILON)
                {
                    const double velRat(maxVelocity / avgV);
                    std::map<double, int> temp;
                    for (const auto &sbin : subcyclingBins)
                    {
                        temp.emplace(fabs(velRat - sbin), sbin);
                    }
                    return temp.begin()->second;
                }
                else
                {
                    return *subcyclingBins.rbegin();
                }

                // if (relV >= 0.1)
                // {
                //     return 1;
                // }
                // else if (relV >= 0.01)
                // {
                //     return 10;
                // }
                // else if (relV >= 0.001)
                // {
                //     return 100;
                // }
                // else
                // {
                //     return 1000;
                // }
            }
            else
            {
                return *subcyclingBins.begin();
            }
        }

        // return DislocationSegmentIO<dim>::velocityGroup(this->source->get_V(),this->sink->get_V(),maxVelocity);
    }

    template <int dim, short unsigned int corder>
    const typename DislocationSegment<dim,corder>::VectorDim& DislocationSegment<dim,corder>::glidePlaneNormal() const
    {
        return this->glidePlanes().size()==1? (*this->glidePlanes().begin())->unitNormal : zeroVector;
    }
    
    template <int dim, short unsigned int corder>
    typename DislocationSegment<dim,corder>::GlidePlaneContainerType DislocationSegment<dim,corder>::glidePlanes() const
    {
        GlidePlaneContainerType temp;
        for (const auto &ln : this->loopLinks())
        {
            if (ln->periodicPlanePatch())
            {
                temp.emplace(ln->periodicPlanePatch()->glidePlane.get());
            }
        }
        return temp;
    }

    template <int dim, short unsigned int corder>
    typename DislocationSegment<dim,corder>::PlanarMeshFaceContainerType DislocationSegment<dim,corder>::meshFaces() const
    {
        PlanarMeshFaceContainerType temp;
        const PlanarMeshFaceContainerType sourceMeshFaces(this->source->meshFaces());
        const PlanarMeshFaceContainerType sinkMeshFaces(this->sink->meshFaces());

        for (const auto &sourceMeshFace : sourceMeshFaces)
        {
            if (sinkMeshFaces.find(sourceMeshFace) != sinkMeshFaces.end())
            {
                temp.emplace(sourceMeshFace);
            }
        }
        return temp;
    }
    template <int dim, short unsigned int corder>
    bool DislocationSegment<dim,corder>::isOnExternalBoundary() const
    { /*!\returns _isOnExternalBoundarySegment.
       */
        bool _isOnExternalBoundary(false);
        for (const auto &face : meshFaces())
        {
            if (face->regionIDs.first == face->regionIDs.second)
            {
                _isOnExternalBoundary = true;
            }
        }

        return _isOnExternalBoundary;
    }

    template <int dim, short unsigned int corder>
    bool DislocationSegment<dim,corder>::isOnInternalBoundary() const
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

    template <int dim, short unsigned int corder>
    bool DislocationSegment<dim,corder>::isOnBoundary() const
    {
        return isOnExternalBoundary() || isOnInternalBoundary();
    }
    template <int dim, short unsigned int corder>
    typename DislocationSegment<dim,corder>::GrainContainerType DislocationSegment<dim,corder>::grains() const
    {
        GrainContainerType temp;
        for (const auto &glidePlane : glidePlanes())
        {
            temp.insert(&glidePlane->grain);
        }
        return temp;
    }

    template <int dim, short unsigned int corder>
    typename DislocationSegment<dim,corder>::VectorDim DislocationSegment<dim,corder>::bndNormal() const
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

    template <int dim, short unsigned int corder>
    std::vector<std::pair<const GlidePlane<dim> *const, const GlidePlane<dim> *const>> DislocationSegment<dim,corder>::parallelAndCoincidentGlidePlanes(const GlidePlaneContainerType &other) const
    {
        std::vector<std::pair<const GlidePlane<dim> *const, const GlidePlane<dim> *const>> pp;

        for (const auto &plane : glidePlanes())
        {
            for (const auto &otherPlane : other)
            {
                //                    if(plane!=otherPlane && plane->n.cross(otherPlane->n).squaredNorm()==0)
                if (plane->n.cross(otherPlane->n).squaredNorm() == 0)
                { // parallel planes
                    pp.emplace_back(plane, otherPlane);
                }
            }
        }
        return pp;
    }

    template <int dim, short unsigned int corder>
    typename DislocationSegment<dim, corder>::VectorDim DislocationSegment<dim, corder>::snapToGlidePlanes(const VectorDim &P)
    {
        GlidePlaneContainerType gps(glidePlanes());
        switch (gps.size())
        {
        case 0:
        {
            assert(false && "Glide plane size must be larger than 0");
            return P;
            break;
        }
        case 1:
        {
            return (*gps.begin())->snapToPlane(P);
            break;
        }
        case 2:
        {
            const PlanePlaneIntersection<dim> ppi(**gps.begin(), **gps.rbegin());
            assert(ppi.type == PlanePlaneIntersection<dim>::INCIDENT && "Intersection must be incident");
            return ppi.P + (P - ppi.P).dot(ppi.d) * ppi.d;
            break;
        }
        case 3:
        {
            const PlanePlaneIntersection<dim> ppi(**gps.begin(), **gps.rbegin());
            assert(ppi.type == PlanePlaneIntersection<dim>::INCIDENT && "Intersection must be incident");
            const auto iterP(++gps.begin());
            const PlaneLineIntersection<dim> pli((*iterP)->P, (*iterP)->unitNormal, ppi.P, ppi.d);
            assert(pli.type == PlaneLineIntersection<dim>::INCIDENT && "Plane line intersection must be incident");
            return pli.P;
            break;
        }
        default:
        {
            const auto iterPlane1(gps.begin());
            const auto iterPlane2(++gps.begin());
            auto iterPlane3(++(++gps.begin()));

            const PlanePlaneIntersection<dim> ppi(**iterPlane1, **iterPlane2);
            assert(ppi.type == PlanePlaneIntersection<dim>::INCIDENT && "Intersection must be incident");
            const PlaneLineIntersection<dim> pli((*iterPlane3)->P, (*iterPlane3)->unitNormal, ppi.P, ppi.d);
            const VectorDim snappedPos(pli.P);
            assert(pli.type == PlaneLineIntersection<dim>::INCIDENT && "Plane line intersection must be incident");
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
            while (++iterPlane3 != gps.end())
            {
                if (!(*iterPlane3)->contains(snappedPos))
                {
                    std::cout << " Glide Plane " << (*iterPlane3)->P.transpose() << " " << (*iterPlane3)->unitNormal.transpose() << std::endl;
                    assert(false && "Plane must contain the position for glide plane size >=3");
                }
            }
            return snappedPos;
            break;
        }
        }
    }

    template <int dim, short unsigned int corder>
    const typename DislocationSegment<dim,corder>::MatrixDim DislocationSegment<dim,corder>::I=DislocationSegment<dim,corder>::MatrixDim::Identity();
    
    template <int dim, short unsigned int corder>
    const typename DislocationSegment<dim,corder>::VectorDim DislocationSegment<dim,corder>::zeroVector=DislocationSegment<dim,corder>::VectorDim::Zero();
    
    template <int dim, short unsigned int corder>
    double DislocationSegment<dim,corder>::quadPerLength=0.2;
    
    //    template <int dim, short unsigned int corder>
    //    bool DislocationSegment<Derived>::assembleWithTangentProjection=false;
    
    template <int dim, short unsigned int corder>
    int DislocationSegment<dim,corder>::verboseDislocationSegment=0;
        
        
    template class DislocationSegment<3,0>;
    
}
#endif
