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
    template <int dim, short unsigned int corder, typename InterpolationType>
    DislocationSegment<dim,corder,InterpolationType>::DislocationSegment(LoopNetworkType* const net,
                                                                         const std::shared_ptr<NetworkNodeType>& nI,
                                                                         const std::shared_ptr<NetworkNodeType>& nJ) :
    /* init */ NetworkLink<DislocationSegment>(net,nI,nJ)
    /* init */,SplineSegment<dim,corder>(this->source->get_P(),this->sink->get_P())
    /* init */,ConfinedDislocationObjectType(this->source->get_P(),this->sink->get_P())
    /* init */,Burgers(VectorDim::Zero())
    /* init */,BurgersNorm(Burgers.norm())
    /* init */,straight(this->source->get_P(),this->sink->get_P(),Burgers,this->chordLength(),this->unitDirection())
    /* init */,_slipSystem(nullptr)
    {
        VerboseDislocationSegment(1,"Constructing DislocationSegment "<<this->tag()<<std::endl);
    }
    
    // template <int dim, short unsigned int corder, typename InterpolationType>
    // void DislocationSegment<dim,corder,InterpolationType>::updateSlipSystem()
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

    // template <int dim, short unsigned int corder, typename InterpolationType>
    // void DislocationSegment<dim, corder, InterpolationType>::updateSlipSystem()
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
       template <int dim, short unsigned int corder, typename InterpolationType>
    void DislocationSegment<dim, corder, InterpolationType>::updateSlipSystem()
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
    template <int dim, short unsigned int corder, typename InterpolationType>
    void DislocationSegment<dim,corder,InterpolationType>::addLoopLink(LoopLinkType* const pL)
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
                this->confinedObject().addGlidePlane(periodicPlanePatch->glidePlane.get());
            }
        }
        else
        {
            this->confinedObject().addGlidePlane(pL->loop->glidePlane.get());
        }


        updateSlipSystem();
    }
    
    
    
    /**********************************************************************/
    template <int dim, short unsigned int corder, typename InterpolationType>
    void DislocationSegment<dim,corder,InterpolationType>::removeLoopLink(LoopLinkType* const pL)
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
        
        this->confinedObject().clear();
        for(const auto& loopLink : this->loopLinks())
        {
            // std::cout<<"loopLink"<<loopLink->tag()<<" Adding glidePlane "<<std::endl;
            if(this->network().simulationParameters.isPeriodicSimulation())
            {
                const auto periodicPlanePatch(loopLink->periodicPlanePatch());
                if(periodicPlanePatch)
                {
                    this->confinedObject().addGlidePlane(periodicPlanePatch->glidePlane.get());
                }
            }
            else
            {
                this->confinedObject().addGlidePlane(loopLink->loop->glidePlane.get());
            }
        }
        
        updateSlipSystem();
        
    }
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    void DislocationSegment<dim,corder,InterpolationType>::updateGeometry()
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
        this->updateConfinement(this->source->get_P(),this->sink->get_P());
        straight.updateGeometry();
        VerboseDislocationSegment(2,"DislocationSegment "<<this->tag()<<", updateGeometry DONE"<<std::endl;);
    }
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    bool DislocationSegment<dim,corder,InterpolationType>::isGlissile() const
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
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    bool DislocationSegment<dim,corder,InterpolationType>::isSessile() const
    {
        return !isVirtualBoundarySegment() && !isGlissile();
    }
    
    /******************************************************************/
    template <int dim, short unsigned int corder, typename InterpolationType>
    void DislocationSegment<dim,corder,InterpolationType>::initFromFile(const std::string& fileName)
    {
        //        LinkType::alpha=TextFileParser(fileName).readScalar<double>("parametrizationExponent",true);
        //        assert((LinkType::alpha)>=0.0 && "parametrizationExponent MUST BE >= 0.0");
        //        assert((LinkType::alpha)<=1.0 && "parametrizationExponent MUST BE <= 1.0");
        quadPerLength=TextFileParser(fileName).readScalar<double>("quadPerLength",true);
        //            assembleWithTangentProjection=TextFileParser(fileName).readScalar<int>("assembleWithTangentProjection",true);
        assert((NetworkLinkType::quadPerLength)>=0.0 && "quadPerLength MUST BE >= 0.0");
        verboseDislocationSegment=TextFileParser("inputFiles/DD.txt").readScalar<int>("verboseDislocationSegment",true);
    }
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    typename DislocationSegment<dim,corder,InterpolationType>::MeshLocation DislocationSegment<dim,corder,InterpolationType>::meshLocation() const
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
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    bool DislocationSegment<dim,corder,InterpolationType>::isVirtualBoundarySegment() const
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
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    bool DislocationSegment<dim,corder,InterpolationType>::hasZeroBurgers() const
    {
        return BurgersNorm<FLT_EPSILON;
    }
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    const std::shared_ptr<SlipSystem>&  DislocationSegment<dim,corder,InterpolationType>::slipSystem() const
    {
        return _slipSystem;
    }
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    bool DislocationSegment<dim,corder,InterpolationType>::isBoundarySegment() const
    {/*!\returns true if both nodes are boundary nodes, and the midpoint is
      * on the boundary.
      */
        return this->isOnExternalBoundary();
    }
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    bool DislocationSegment<dim,corder,InterpolationType>::isGrainBoundarySegment() const
    {
        return this->isOnInternalBoundary();
    }
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    const typename DislocationSegment<dim,corder,InterpolationType>::VectorDim& DislocationSegment<dim,corder,InterpolationType>::burgers() const
    {
        return Burgers;
    }
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    void DislocationSegment<dim,corder,InterpolationType>::addToGlobalAssembly(std::deque<Eigen::Triplet<double> >& kqqT,
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
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    void DislocationSegment<dim,corder,InterpolationType>::updateQuadraturePointsSeg()
    {
        this->updateQuadraturePoints(*this,quadPerLength,false);
    }

    template <int dim, short unsigned int corder, typename InterpolationType>
    void DislocationSegment<dim,corder,InterpolationType>::assembleGlide()
    {
        VerboseDislocationSegment(2,"DislocationSegment "<<this->tag()<<", assembleGlide"<<std::endl;);
        this->updateForcesAndVelocities(*this);
        Fq= this->quadraturePoints().size()? this->nodalVelocityVector() : VectorNdof::Zero();
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
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    const typename DislocationSegment<dim,corder,InterpolationType>::VectorDim& DislocationSegment<dim,corder,InterpolationType>::glidePlaneNormal() const
    {
        return this->glidePlanes().size()==1? (*this->glidePlanes().begin())->unitNormal : zeroVector;
    }
    
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    const typename DislocationSegment<dim,corder,InterpolationType>::MatrixDim DislocationSegment<dim,corder,InterpolationType>::I=DislocationSegment<dim,corder,InterpolationType>::MatrixDim::Identity();
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    const typename DislocationSegment<dim,corder,InterpolationType>::VectorDim DislocationSegment<dim,corder,InterpolationType>::zeroVector=DislocationSegment<dim,corder,InterpolationType>::VectorDim::Zero();
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    double DislocationSegment<dim,corder,InterpolationType>::quadPerLength=0.2;
    
    //    template <int dim, short unsigned int corder, typename InterpolationType>
    //    bool DislocationSegment<Derived>::assembleWithTangentProjection=false;
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    int DislocationSegment<dim,corder,InterpolationType>::verboseDislocationSegment=0;
    
}
#endif
