/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DISLOCATIONNETWORKREMESH_H_
#define model_DISLOCATIONNETWORKREMESH_H_

#include <chrono>
#include <utility>  // for std::pair and "<" operator between them
#include <set>      // for std::set
#include <vector>      // for std::set
#include <assert.h>

#include <Eigen/Dense>
//#include <EdgeFinder.h>
//#include <GramSchmidt.h>
#include <TerminalColors.h>
#include <MPIcout.h>
#include <Simplex.h>
#include <LatticeMath.h>

namespace model
{
    
    /*! \brief Class template that handles the nodal remesh of the DislocationNetwork.
     */
    template <typename DislocationNetworkType>
    class DislocationNetworkRemesh
    {
        
        enum{dim=DislocationNetworkType::dim};
        typedef Eigen::Matrix<double,dim,1> VectorDimD;
        
        typedef typename DislocationNetworkType::LinkType LinkType;
        //        typedef typename DislocationNetworkType::NetworkLinkContainerType NetworkLinkContainerType;
        typedef typename DislocationNetworkType::NodeType NodeType;
        typedef LatticeVector<dim> LatticeVectorType;
        typedef typename DislocationNetworkType::IsNetworkLinkType IsNetworkLinkType;
        typedef typename DislocationNetworkType::IsConstNetworkLinkType IsConstNetworkLinkType;
        
        //! A reference to the DislocationNetwork
        DislocationNetworkType& DN;
        
        /**********************************************************************/
        std::pair<bool,double> mustBeContracted(const LinkType& segment) const
        {
            const double vTolcont=0.0;
            
            std::pair<bool,double> temp=std::make_pair(false,0.0);
            const VectorDimD chord(segment.chord()); // this is sink->get_P() - source->get_P()
            const double chordLength(chord.norm());
            const VectorDimD dv(segment.sink->get_V()-segment.source->get_V());
            bool endsAreApproaching( chord.dot(dv) < 0.0 );
            double currentLmin=Lmin;
            const int sourceNonZeroSize(segment.source->nonZeroNeighbors().size());
            const int   sinkNonZeroSize(segment.sink->nonZeroNeighbors().size());
            if((sourceNonZeroSize>2 && sinkNonZeroSize==2) || (sourceNonZeroSize==2 && sinkNonZeroSize>2))
            {
                currentLmin=1.0;
            }
            
            
            
            if (((endsAreApproaching || segment.isBoundarySegment())// ends are approaching
                 && dv.norm()*DN.simulationParameters.dt>vTolcont*chordLength // contraction is large enough compared to segment length
                 && chordLength<currentLmin // segment is small
                 )
                //                || segment.isSimpleBndSegment()
                //|| segment.isSimpleSessile()
                )
            {
                temp=std::make_pair(true,chordLength);
            }
            return temp;
        }
        
        /**********************************************************************/
        void remeshByRemoval()
        {
            const auto t0= std::chrono::system_clock::now();
            model::cout<<"		Remeshing network: removing... "<<std::flush;
            
            std::deque<size_t> toBeRemoved;
            for(const auto& node : DN.nodes())
            {
                //                std::cout<<"node "<<node.second->sID<<" "<<node.second->isSimpleBoundaryNode()<<" "<<node.second->isSimpleGrainBoundaryNode()<<std::endl;
                if(!node.second->masterNode)
                {
                    if(node.second->isRemovable(Lmin,cosRemove))
                    {
                        for(const auto& imagePair : node.second->imageSharedNodeContainer)
                        {
                            toBeRemoved.push_back(imagePair.second->sID);
                        }
                        toBeRemoved.push_back(node.second->sID); // insert image IDs BEFORE master ID, so that images are removed FIRST
                    }
                }
            }
            
            size_t Nremoved=0;
            for(const auto& nodeID : toBeRemoved)
            {
                const auto isNode=DN.node(nodeID);
                if(isNode.first)
                {// Removing may have deleted the node, check that it exists
                    if(isNode.second->isRemovable(Lmin,cosRemove))
                    {// Removing may have altered the isSimpleBoundaryNode, check again
                        Nremoved+=DN.remove(nodeID);
                    }
//                    if(!isNode.second->masterNode)
//                    {
//                        if(isNode.second->isRemovable(Lmin,cosRemove))
//                        {// Removing may have altered the isSimpleBoundaryNode, check again
//                            Nremoved+=DN.remove(nodeID);
//                        }
//                    }
                }
            }
            
            model::cout<<" ("<<Nremoved<<" removed)"<<std::flush;
            model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<defaultColor<<std::endl;
        }
        
        
//        /**********************************************************************/
//        void remeshByRemoval()
//        {
//            const auto t0= std::chrono::system_clock::now();
//            model::cout<<"		remeshing network: removing... "<<std::flush;
//            
//            std::deque<size_t> toBeRemoved;
//            for(const auto& node : DN.nodes())
//            {
////                std::cout<<"node "<<node.second->sID<<" "<<node.second->isSimpleBoundaryNode()<<" "<<node.second->isSimpleGrainBoundaryNode()<<std::endl;
//                if(   !node.second->isVirtualBoundaryNode
//                   && (node.second->isSimpleBoundaryNode() || node.second->isSimpleGrainBoundaryNode() || node.second->isSimpleSessileNode())
//                   )
//                {
//                    toBeRemoved.push_back(node.second->sID);
//                }
//            }
//            
//            size_t Nremoved=0;
//            for(const auto& nodeID : toBeRemoved)
//            {
//                const auto isNode=DN.node(nodeID);
//                if(isNode.first)
//                {// Removing may have deleted the node, check that it exists
//                    if(   isNode.second->isSimpleBoundaryNode()
//                       || isNode.second->isSimpleGrainBoundaryNode()
//                       || isNode.second->isSimpleSessileNode()
//                       || isNode.second->isRemovable(Lmin,cosRemove)
//                       )
//                    {// Removing may have altered the isSimpleBoundaryNode, check again
//                        Nremoved+=DN.remove(nodeID);
//                    }
//                }
//            }
//            
//            model::cout<<" ("<<Nremoved<<" removed)"<<std::flush;
//            model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<defaultColor<<std::endl;
//        }
        
        
//        /**********************************************************************/
//        void remeshByContraction()
//        {/*! Contract edges according to two criteria.
//          */
//            const auto t0= std::chrono::system_clock::now();
//            model::cout<<"		remeshing network: contracting... "<<std::flush;
//            
//            //            model::cout<<"contracting..."<<std::flush;
//            
//            
//            // Store segments to be contracted
//            std::set<std::pair<double,std::pair<size_t,size_t> > > toBeContracted; // order by increasing segment length
//            
//            for (const auto& linkIter : DN.links())
//            {
//                const LinkType& segment(*linkIter.second);
//                std::pair<bool,double> temp=mustBeContracted(segment);
//                if(temp.first)
//                {
//                    const bool inserted=toBeContracted.insert(std::make_pair(temp.second,segment.nodeIDPair)).second;
//                    assert(inserted && "COULD NOT INSERT IN SET.");
//                }
//            }
//            
//            // Call Network::contract
//            unsigned int Ncontracted(0);
//            for (const auto& smallIter : toBeContracted)
//            {
//                const size_t i(smallIter.second.first);
//                const size_t j(smallIter.second.second);
//                const IsConstNetworkLinkType Lij(DN.link(i,j));
//                
//                if (Lij.first )
//                {// previous contractions could have destroyed Lij, so check that Lij exists
//                    if(mustBeContracted(*Lij.second).first)
//                    {// previous contractions could have changed Lij, so check the contract conditions again
//                        Ncontracted+=DN.contract(Lij.second->source,Lij.second->sink);
//                    }
//                }
//            }
//            model::cout<<" ("<<Ncontracted<<" contracted)"<<std::flush;
//            model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<defaultColor<<std::endl;
//            
//        }
        
        /**********************************************************************/
        void remeshByExpansion()
        {
            const auto t0= std::chrono::system_clock::now();
            model::cout<<"		Remeshing network: expanding... "<<std::flush;
            
            //            model::cout<<"expanding..."<<std::flush;
            
//            const double cos_theta_max_crit = std::cos(M_PI-thetaDeg*M_PI/180.0);  /*critical angle */
            std::set<std::pair<size_t,size_t> > toBeExpanded;
            
            // Store the segments to be expanded
            //            for (typename NetworkLinkContainerType::const_iterator linkIter=DN.linkBegin();linkIter!=DN.linkEnd();++linkIter)
            for (const auto& linkIter : DN.links())
            {
                
                
                if( !linkIter.second->hasZeroBurgers()
                   //&& !linkIter.second->isSimpleSessile())
                   && !linkIter.second->isSessile()
                   && !linkIter.second->isBoundarySegment()
                   && !linkIter.second->isGrainBoundarySegment()
                   && !linkIter.second->isVirtualBoundarySegment()
                   //                   && !linkIter.second->isSimpleBndSegment()
                   )
                    
                {
//                    std::cout<<"Expanding "<<linkIter.second->source->sID<<"->"<<linkIter.second->sink->sID<<std::endl;
//                    std::cout<<linkIter.second->hasZeroBurgers()<<std::endl;
//                    std::cout<<linkIter.second->isSessile()<<std::endl;
//                    std::cout<<linkIter.second->isBoundarySegment()<<std::endl;
//                    std::cout<<linkIter.second->source->isBoundaryNode()<<std::endl;
//                    std::cout<<linkIter.second->sink->isBoundaryNode()<<std::endl;
//                    std::cout<<linkIter.second->boundingBoxSegments().contains(0.5*(linkIter.second->source->get_P()+linkIter.second->sink->get_P())).first<<std::endl;
                    
                    const VectorDimD chord(linkIter.second->chord()); // this is sink->get_P() - source->get_P()
                    const double chordLength(chord.norm());
                    //				const VectorDimD dv(linkIter->second.sink->get_V()-linkIter->second.source->get_V());
                    
                    
                    // Always expand single FR source segment
                    //                    if (linkIter.second->source->openOrder()==1 && linkIter.second->sink->openOrder()==1)
                    //                    {
                    //                        toBeExpanded.insert(linkIter->second.nodeIDPair);
                    //                    }
                    
                    // Expand pin points
                    if (   linkIter.second->source->constraintNormals().size()>2
                        && linkIter.second->  sink->constraintNormals().size()>2
                        && chordLength>3.0*Lmin)
                    {
                        toBeExpanded.insert(linkIter.second->nodeIDPair);
                    }
                    
                    if (!linkIter.second->source->isSimple() && !linkIter.second->sink->isSimple()
                        /*&& chord.dot(dv)>vTolexp*chordLength*dv.norm()*/ && chordLength>3.0*Lmin)
                    { // also expands a straight line to generate glissile segment
                        toBeExpanded.insert(linkIter.second->nodeIDPair);
                    }
                    
                    // Expand segments shorter than Lmax
                    if (chordLength>Lmax)
                    {
                        toBeExpanded.insert(linkIter.second->nodeIDPair);
                    }
                    
                    //                    // Angle criterion
                    //                    if (linkIter.second->source->is_simple())
                    //                    { //check angle criterion at source
                    //                        const VectorDimD c0(linkIter.second->source->openNeighborNode(0)->get_P()-linkIter.second->source->get_P());
                    //                        const VectorDimD c1(linkIter.second->source->openNeighborNode(1)->get_P()-linkIter.second->source->get_P());
                    //                        const double c0norm(c0.norm());
                    //                        const double c1norm(c1.norm());
                    //                        if(c0.dot(c1)>cos_theta_max_crit*c0norm*c1norm)
                    //                        {
                    //                            if (c0norm>3.0*Lmin /*&& c0.dot(v0)>vTolexp*c0norm*v0.norm()*/)
                    //                            {
                    //                                toBeExpanded.insert(linkIter.second->source->openNeighborLink(0)->nodeIDPair);
                    //                            }
                    //
                    //                            if (c1norm>3.0*Lmin /*&& c1.dot(v1)>vTolexp*c1norm*v1.norm()*/)
                    //                            {
                    //                                toBeExpanded.insert(linkIter.second->source->openNeighborLink(1)->nodeIDPair);
                    //                            }
                    //                        }
                    //                    }
                    //                    if (linkIter.second->sink->is_simple())
                    //                    { //check angle criterion at sink
                    //                        const VectorDimD c0(linkIter.second->sink->openNeighborNode(0)->get_P()-linkIter.second->sink->get_P());
                    //                        const VectorDimD c1(linkIter.second->sink->openNeighborNode(1)->get_P()-linkIter.second->sink->get_P());
                    //                        const double c0norm(c0.norm());
                    //                        const double c1norm(c1.norm());
                    //                        if(c0.dot(c1)>cos_theta_max_crit*c0norm*c1norm)
                    //                        {
                    //                            if (c0norm>3.0*Lmin /*&& c0.dot(v0)>vTolexp*c0norm*v0.norm()*/)
                    //                            {
                    //                                toBeExpanded.insert(linkIter.second->sink->openNeighborLink(0)->nodeIDPair);
                    //                            }
                    //                            if (c1norm>3.0*Lmin/* && c1.dot(v1)>vTolexp*c1norm*v1.norm()*/)
                    //                            {
                    //                                //														model::cout<<"Expanding 4"<<std::endl;
                    //                                toBeExpanded.insert(linkIter.second->sink->openNeighborLink(1)->nodeIDPair);
                    //                            }
                    //                        }
                    //                    }
                }
                
            }
            
            
            // Call Network::expand
            unsigned int Nexpanded(0);
            const double expand_at(0.5);
            for (std::set<std::pair<size_t,size_t> >::const_iterator expIter=toBeExpanded.begin(); expIter!=toBeExpanded.end(); ++expIter)
            {
                const size_t i(expIter->first);
                const size_t j(expIter->second);
                const IsConstNetworkLinkType Lij(DN.link(i,j));
                if(Lij.first)
                {
                    //                    std::cout<<"Expanding "<<i<<"->"<<j<<std::endl;
                    //VectorDimD expandPoint(Lij.second->get_r(expand_at));
                    VectorDimD expandPoint(Lij.second->get_r(expand_at));
                    
                    //                    LatticeVectorType expandPoint(Lij.second->glidePlane->snapToLattice(Lij.second->get_r(expand_at)));
                    //                    auto simplexCheckPair=DN.pointIsInsideMesh(expandPoint.cartesian(),Lij.second->source->includingSimplex());
                    //                    if(Lij.second->isSessile())
                    //                    {
                    //                        PlanePlaneIntersection ppi(Lij.second->glidePlane,Lij.second->sessilePlane);
                    //                        LatticeLine line(ppi.P,ppi.d);
                    //                        expandPoint=line.snapToLattice(expandPoint.cartesian());
                    //
                    //                    }
                    //                    else
                    //                    {
                    //                        if(simplexCheckPair.first && simplexCheckPair.second->region->regionID!=DN.node(i).second->grain.grainID)
                    //                        {
                    //                            const LatticePlane& GBplane(DN.shared.poly.grainBoundary(simplexCheckPair.second->region->regionID,DN.node(i).second->grain.grainID).latticePlane(DN.node(i).second->grain.grainID));
                    //                            const LatticePlane& glidePlane(Lij.second->glidePlane);
                    //                            const PlanePlaneIntersection ppi(GBplane,glidePlane);
                    //                            const LatticeLine line(ppi.P,ppi.d);
                    //                            expandPoint=line.snapToLattice(expandPoint.cartesian());
                    //
                    //                        }
                    //                    }
                    //
                    //                    simplexCheckPair=DN.pointIsInsideMesh(expandPoint.cartesian(),Lij.second->source->includingSimplex());
                    //                    assert(simplexCheckPair.second->region->regionID==DN.node(i).second->grain.grainID && simplexCheckPair.second->region->regionID==DN.node(j).second->grain.grainID && "EXPAND POINT IN INCORRECT REGION.");
                    
                    //                    if(!Lij.second->isSessile)
                    //                    {
                    //                        expandPoint=Lij.second->glidePlane.snapToLattice(expandPoint).cartesian();
                    //                    }
                    //                    else
                    //                    {
                    //                        PlanePlaneIntersection ppi(Lij.second->glidePlane,Lij.second->sessilePlane);
                    //                        LatticeLine line(ppi.P,ppi.d);
                    //                        expandPoint=line.snapToLattice(expandPoint).cartesian();
                    //                    }
                    
                    
                    
                    if(  (expandPoint-DN.node(i).second->get_P()).squaredNorm()
                       &&(expandPoint-DN.node(j).second->get_P()).squaredNorm() )
                    {
                        //                        if(simplexCheckPair.first)
                        //                        {
                        //                            //                        DN.expand(i,j,expandPoint);
                        //std::cout<<"Expanding "<<i<<"->"<<j<<std::endl;
                        DN.expand(i,j,expand_at);
                        Nexpanded++;
                        //                        }
                    }
                }
            }
            model::cout<<" ("<<Nexpanded<<" expanded)"<<std::flush;
            model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<defaultColor<<std::endl;
            
        }
        
        /**********************************************************************/
        void contract0chordSegments()
        {
            model::cout<<"		Contracting zero-chord segments... "<<std::flush;
            const auto t0= std::chrono::system_clock::now();
            
            std::set<std::pair<double,std::pair<size_t,size_t> > > toBeContracted; // order by increasing segment length
            
            //            for (typename NetworkLinkContainerType::const_iterator linkIter=DN.linkBegin();linkIter!=DN.linkEnd();++linkIter)
            for (const auto& linkIter : DN.links())
            {
                VectorDimD chord(linkIter.second->chord()); // this is sink->get_P() - source->get_P()
                double chordLength(chord.norm());
                if (chordLength<=FLT_EPSILON
                    && (linkIter.second->source->get_V()-linkIter.second->  sink->get_V()).squaredNorm()<FLT_EPSILON)
                {// toBeContracted part
                    toBeContracted.insert(std::make_pair(chordLength,linkIter.second->nodeIDPair));
                }
            }
            
            // Call Network::contract
            unsigned int Ncontracted(0);
            for (const auto& smallIter : toBeContracted)
            {
                const size_t i(smallIter.second.first);
                const size_t j(smallIter.second.second);
                const IsConstNetworkLinkType Lij(DN.link(i,j));
                
                if (Lij.first )
                {
                    Ncontracted+=DN.contract(Lij.second->source,Lij.second->sink);
                }
            }
            model::cout<<"("<<Ncontracted<<" contracted)"<<std::flush;
            model::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
            
        }
        
        /**********************************************************************/
        static double minMeshSize(const SimplicialMesh<dim>& mesh)
        {
            return std::min(mesh.xMax(0)-mesh.xMin(0),std::min(mesh.xMax(1)-mesh.xMin(1),mesh.xMax(2)-mesh.xMin(2)));
        }
        
    public:
        
        const double Lmax;
        const double Lmin;
        const double cosRemove;
//        const double thetaDeg;
//        const double neighborRadius;
        const short unsigned int remeshFrequency;
        
        /**********************************************************************/
        DislocationNetworkRemesh(DislocationNetworkType& DN_in) :
        /* init */ DN(DN_in)
        /* init */,Lmax(TextFileParser("inputFiles/DD.txt").readScalar<double>("Lmax",true)*minMeshSize(DN.mesh))
        /* init */,Lmin(TextFileParser("inputFiles/DD.txt").readScalar<double>("Lmin",true)*minMeshSize(DN.mesh))
        /* init */,cosRemove(cos(TextFileParser("inputFiles/DD.txt").readScalar<double>("nodeRemoveAngleDeg",true)*M_PI/180.0))
//        /* init */,thetaDeg;
//        /* init */,neighborRadius;
        /* init */,remeshFrequency(TextFileParser("inputFiles/DD.txt").readScalar<int>("remeshFrequency",true))
        {
            
            assert(Lmin<=Lmax);
            assert(Lmax>3.0*Lmin);
            assert(Lmin>=0.0);
//            assert(Lmin>=2.0*DDtimeIntegrator<0>::dxMax && "YOU MUST CHOOSE Lmin>2*dxMax.");

        }
        
        
        /**********************************************************************/
        void remesh(const long int& runID)
        {/*! Performs remeshByContraction and then remeshByExpansion.
          * This order guarantees that 2-vertex NetworkComponents are expanded.
          */
            if (remeshFrequency)
            {
                if(!(runID%remeshFrequency))
                {
                    remeshByRemoval();
//                    remeshByContraction();
                    remeshByExpansion();
                    contract0chordSegments();
                }
            }
        }
        
        
        
        
        
        
        
        /**********************************************************************/
        void loopInversion(const double& dt)
        {
            const auto t0= std::chrono::system_clock::now();
            model::cout<<"		Checking for loop inversions ... "<<std::flush;
            //! 3- Check and remove loop inversions
            std::vector<int> toBeErased;
            for (typename DislocationNetworkType::NetworkComponentContainerType::iterator snIter=DN.ABbegin(); snIter!=DN.ABend();++snIter)
            {
                typename DislocationNetworkType::DislocationNetworkComponentType dnC(*snIter->second);
                
                if (dnC.loopInversion(dt))
                {
                    model::cout<<"NetworkComponent "<<snIter->second->sID<<" containing "<<snIter->second->nodeOrder()<<" is an inverted loop"<<std::endl;
                    for (typename DislocationNetworkType::NetworkComponentType::NetworkComponentNodeContainerType::const_iterator nodeIter=snIter->second->nodeBegin();nodeIter!=snIter->second->nodeEnd();++nodeIter)
                    {
                        toBeErased.push_back(nodeIter->second->sID);
                    }
                }
            }
            model::cout<<" found "<<toBeErased.size()<<" inverted nodes ... ";
            for (unsigned int nn=0;nn<toBeErased.size();++nn)
            {
                DN.template removeVertex<true>(toBeErased[nn]);
            }
            model::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
        }
        
        
    };
    
    // static data
//    template <typename DislocationNetworkType>
//    double DislocationNetworkRemesh<DislocationNetworkType>::Lmin=10.0;
//
//    template <typename DislocationNetworkType>
//    double DislocationNetworkRemesh<DislocationNetworkType>::Lmax=100.0;
//
//    template <typename DislocationNetworkType>
//    double DislocationNetworkRemesh<DislocationNetworkType>::cosRemove=0.7071;
//
//    template <typename DislocationNetworkType>
//    double DislocationNetworkRemesh<DislocationNetworkType>::thetaDeg=45.0;
//
//    template <typename DislocationNetworkType>
//    double DislocationNetworkRemesh<DislocationNetworkType>::neighborRadius=0.001;
//
//    template <typename DislocationNetworkType>
//    short unsigned int DislocationNetworkRemesh<DislocationNetworkType>::remeshFrequency=1;
    
    
} // namespace model
#endif

