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
#include <model/Network/Operations/EdgeFinder.h>
//#include <model/Math/GramSchmidt.h>
#include <model/Utilities/TerminalColors.h>
#include <model/MPI/MPIcout.h>
#include <model/Mesh/Simplex.h>
#include <model/LatticeMath/LatticeMath.h>
#define VerboseRemesh(N,x) if(verboseRemesh>=N){model::cout<<x;}

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
        typedef typename DislocationNetworkType::NetworkLinkContainerType NetworkLinkContainerType;
        typedef typename DislocationNetworkType::NodeType NodeType;
        typedef LatticeVector<dim> LatticeVectorType;
        //! A reference to the DislocationNetwork
        DislocationNetworkType& DN;
        
    public:
        
        static double Lmax;
        static double Lmin;
        static double thetaDeg;
        static double neighborRadius;
        static int verboseRemesh;
        
        /**********************************************************************/
        DislocationNetworkRemesh(DislocationNetworkType& DN_in) :
        /* init list */ DN(DN_in)
        {
        }
        
        
        /**********************************************************************/
        void remesh()
        {/*! Performs remeshByContraction and then remeshByExpansion.
          * This order guarantees that 2-vertex NetworkComponents are expanded.
          */
            const auto t0= std::chrono::system_clock::now();
            model::cout<<"		Remeshing network... "<<std::flush;
            remeshByContraction();
            remeshByExpansion();
            model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<defaultColor<<std::endl;
        }
        
        /**********************************************************************/
        void remeshByContraction()
        {/*! Contract edges according to two criteria.
          */
            model::cout<<"contracting..."<<std::flush;
            
            const double vTolcont=0.0;
            
            // Store segments to be contracted
            std::set<std::pair<double,std::pair<size_t,size_t> > > toBeContracted; // order by increasing segment length
//            for (typename NetworkLinkContainerType::const_iterator linkIter=DN.linkBegin();linkIter!=DN.linkEnd();++linkIter)
//            {
//                const VectorDimD chord(linkIter->second.chord()); // this is sink->get_P() - source->get_P()
//                const double chordLength(chord.norm());
//                const VectorDimD dv(linkIter->second.sink->get_V()-linkIter->second.source->get_V());
//                //				bool endsAreApproaching( chord.dot(dv) < -vTolcont*chordLength*dv.norm() );
//                bool endsAreApproaching( chord.dot(dv) < 0.0 );
//                //				bool endsAreApproaching( chord.dot(dv) + dv.squaredNorm()*dt < 0.0 );
//                if ((endsAreApproaching || linkIter->second.is_boundarySegment())// ends are approaching
//                    && dv.norm()*DN.get_dt()>vTolcont*chordLength // contraction is large enough compared to segment length
//                    && chordLength<Lmin)
//                {
//                    assert(toBeContracted.insert(std::make_pair(chordLength,linkIter->second.nodeIDPair)).second && "COULD NOT INSERT IN SET.");
//                }
//            }

            for (const auto& linkIter : DN.links())
            {
                const LinkType& segment(linkIter.second);
                const VectorDimD chord(segment.chord()); // this is sink->get_P() - source->get_P()
                const double chordLength(chord.norm());
                const VectorDimD dv(segment.sink->get_V()-segment.source->get_V());
                bool endsAreApproaching( chord.dot(dv) < 0.0 );

                if ((endsAreApproaching || segment.is_boundarySegment())// ends are approaching
                    && dv.norm()*DN.get_dt()>vTolcont*chordLength // contraction is large enough compared to segment length
                    && chordLength<Lmin // segment is small
                    )
                {
                    VerboseRemesh(1,"contract candidate "<<segment.nodeIDPair.first<<" "<<segment.nodeIDPair.second<<std::endl;);
                    assert(toBeContracted.insert(std::make_pair(chordLength,segment.nodeIDPair)).second && "COULD NOT INSERT IN SET.");
                }
            }
            
            // Call Network::contract
            unsigned int Ncontracted(0);
//            for (std::set<std::pair<double,std::pair<size_t,size_t> > >::const_iterator smallIter =toBeContracted.begin();
//                 /*                                                                  */ smallIter!=toBeContracted.end();
//                 /*                                                                  */ smallIter++)
//            {
//                const size_t i(smallIter->second.first);
//                const size_t j(smallIter->second.second);
//                const typename EdgeFinder<LinkType>::isNetworkEdgeType Lij(DN.link(i,j));
//                
//                if (Lij.first )
//                {
//                    Ncontracted+=DN.contractWithConstraintCheck(DN.node(i),DN.node(j));
//                }
//            }
            for (const auto& smallIter : toBeContracted)
            {
                const size_t i(smallIter.second.first);
                const size_t j(smallIter.second.second);
                const typename EdgeFinder<LinkType>::isNetworkEdgeType Lij(DN.link(i,j));
                
                if (Lij.first )
                {
                    Ncontracted+=DN.contractWithConstraintCheck(DN.node(i),DN.node(j));
                }
            }
            model::cout<<" ("<<Ncontracted<<" contracted)"<<std::flush;
            
        }
        
        /**********************************************************************/
        void remeshByExpansion()
        {
            model::cout<<"expanding..."<<std::flush;
            
            double cos_theta_max_crit = std::cos(M_PI-thetaDeg*M_PI/180.0);  /*critical angle */
            std::set<std::pair<size_t,size_t> > toBeExpanded;
            
            // Store the segments to be expanded
            for (typename NetworkLinkContainerType::const_iterator linkIter=DN.linkBegin();linkIter!=DN.linkEnd();++linkIter)
            {
                
                const VectorDimD chord(linkIter->second.chord()); // this is sink->get_P() - source->get_P()
                const double chordLength(chord.norm());
                //				const VectorDimD dv(linkIter->second.sink->get_V()-linkIter->second.source->get_V());
                
                
                // Always expand single FR source segment
                if (linkIter->second.source->openOrder()==1 && linkIter->second.sink->openOrder()==1)
                {
                    VerboseRemesh(1,"expand candidate (FR) "<<linkIter->second.nodeIDPair.first<<"->"<<linkIter->second.nodeIDPair.second<<std::endl;);
                    toBeExpanded.insert(linkIter->second.nodeIDPair);
                }
                
                // Expand pin points
                if (   linkIter->second.source->constraintNormals().size()>2
                    && linkIter->second.  sink->constraintNormals().size()>2
                    && chordLength>3.0*Lmin)
                {
                    VerboseRemesh(1,"expand candidate (pinned points)"<<linkIter->second.nodeIDPair.first<<"->"<<linkIter->second.nodeIDPair.second<<std::endl;);
                    toBeExpanded.insert(linkIter->second.nodeIDPair);
                }
                
                if (!linkIter->second.source->is_simple() && !linkIter->second.sink->is_simple()
                    /*&& chord.dot(dv)>vTolexp*chordLength*dv.norm()*/ && chordLength>3.0*Lmin)
                { // also expands a straight line to generate glissile segment

                    VerboseRemesh(1,"expand candidate (not simple vertices)"<<linkIter->second.nodeIDPair.first<<"->"<<linkIter->second.nodeIDPair.second<<std::endl;);
                    toBeExpanded.insert(linkIter->second.nodeIDPair);
                }
                
                // Expand segments shorter than Lmax
                if (chordLength>Lmax)
                {
                    VerboseRemesh(1,"expand candidate (long)"<<linkIter->second.nodeIDPair.first<<"->"<<linkIter->second.nodeIDPair.second<<std::endl;);
                    toBeExpanded.insert(linkIter->second.nodeIDPair);
                }
                
                // Angle criterion
                if (linkIter->second.source->is_simple()
                    && linkIter->second.source->confiningPlanes().size()==1)
                { //check angle criterion at source
                    const VectorDimD c0(linkIter->second.source->openNeighborNode(0)->get_P()-linkIter->second.source->get_P());
                    const VectorDimD c1(linkIter->second.source->openNeighborNode(1)->get_P()-linkIter->second.source->get_P());
                    const double c0norm(c0.norm());
                    const double c1norm(c1.norm());
                    if(c0.dot(c1)>cos_theta_max_crit*c0norm*c1norm)
                    {
                        if (c0norm>3.0*Lmin /*&& c0.dot(v0)>vTolexp*c0norm*v0.norm()*/)
                        {
                            VerboseRemesh(1,"expand candidate (angle 1)"<<linkIter->second.nodeIDPair.first<<"->"<<linkIter->second.nodeIDPair.second<<std::endl;);
                            toBeExpanded.insert(linkIter->second.source->openNeighborLink(0)->nodeIDPair);
                        }
                        
                        if (c1norm>3.0*Lmin /*&& c1.dot(v1)>vTolexp*c1norm*v1.norm()*/)
                        {
                            VerboseRemesh(1,"expand candidate (angle 2)"<<linkIter->second.nodeIDPair.first<<"->"<<linkIter->second.nodeIDPair.second<<std::endl;);
                            toBeExpanded.insert(linkIter->second.source->openNeighborLink(1)->nodeIDPair);
                        }
                    }
                }
                if (linkIter->second.sink->is_simple()
                    && linkIter->second.sink->confiningPlanes().size()==1)
                { //check angle criterion at sink
                    const VectorDimD c0(linkIter->second.sink->openNeighborNode(0)->get_P()-linkIter->second.sink->get_P());
                    const VectorDimD c1(linkIter->second.sink->openNeighborNode(1)->get_P()-linkIter->second.sink->get_P());
                    const double c0norm(c0.norm());
                    const double c1norm(c1.norm());
                    if(c0.dot(c1)>cos_theta_max_crit*c0norm*c1norm)
                    {
                        if (c0norm>3.0*Lmin /*&& c0.dot(v0)>vTolexp*c0norm*v0.norm()*/)
                        {
                            VerboseRemesh(1,"expand candidate (angle 3)"<<linkIter->second.nodeIDPair.first<<"->"<<linkIter->second.nodeIDPair.second<<std::endl;);
                            toBeExpanded.insert(linkIter->second.sink->openNeighborLink(0)->nodeIDPair);
                        }
                        if (c1norm>3.0*Lmin/* && c1.dot(v1)>vTolexp*c1norm*v1.norm()*/)
                        {
                            //														model::cout<<"Expanding 4"<<std::endl;
                            VerboseRemesh(1,"expand candidate (angle 4)"<<linkIter->second.nodeIDPair.first<<"->"<<linkIter->second.nodeIDPair.second<<std::endl;);
                            toBeExpanded.insert(linkIter->second.sink->openNeighborLink(1)->nodeIDPair);
                        }
                    }
                }
                
            }
            
            
            // Call Network::expand
            unsigned int Nexpanded(0);
            const double expand_at(0.5);
            for (std::set<std::pair<size_t,size_t> >::const_iterator expIter=toBeExpanded.begin(); expIter!=toBeExpanded.end(); ++expIter)
            {
                const size_t i(expIter->first);
                const size_t j(expIter->second);
                const typename EdgeFinder<LinkType>::isNetworkEdgeType Lij(DN.link(i,j));
                if(Lij.first)
                {
                   //std::cout<<"Expanding "<<i<<"->"<<j<<std::endl;
                    VectorDimD expandPoint(Lij.second->get_r(expand_at));
                    
                    if(!Lij.second->isSessile)
                    {
                        expandPoint=Lij.second->glidePlane.snapToLattice(expandPoint);
                    }
                    else
                    {
                        PlanePlaneIntersection ppi(Lij.second->glidePlane,Lij.second->sessilePlane);
                        LatticeLine line(ppi.P,ppi.d);
                        expandPoint=line.snapToLattice(expandPoint);
                    }

                    if(  (LatticeVectorType(expandPoint)-DN.node(i).second->get_L()).squaredNorm()
                       &&(LatticeVectorType(expandPoint)-DN.node(j).second->get_L()).squaredNorm() )
                    {
                        if(DN.pointIsInsideMesh(expandPoint,Lij.second->source->includingSimplex()).first)
                        {
                            //                        DN.expand(i,j,expandPoint);
                            DN.expand(i,j,LatticeVectorType(expandPoint));
                            Nexpanded++;
                        }
                    }
                }
            }
            model::cout<<" ("<<Nexpanded<<" expanded)"<<std::flush;
        }
        
        /**********************************************************************/
        void contract0chordSegments()
        {
            model::cout<<"		contracting zero-chord segments... "<<std::flush;
            const auto t0= std::chrono::system_clock::now();
            
            std::set<std::pair<double,std::pair<size_t,size_t> > > toBeContracted; // order by increasing segment length
            
            for (typename NetworkLinkContainerType::const_iterator linkIter=DN.linkBegin();linkIter!=DN.linkEnd();++linkIter)
            {
                VectorDimD chord(linkIter->second.chord()); // this is sink->get_P() - source->get_P()
                double chordLength(chord.norm());
                if (chordLength<=FLT_EPSILON)
                {// toBeContracted part
                    VerboseRemesh(1,"0-chord contract candidate "<<linkIter->second.nodeIDPair.first<<" "<<linkIter->second.nodeIDPair.second<<std::endl;);
                    toBeContracted.insert(std::make_pair(chordLength,linkIter->second.nodeIDPair));
                }
            }
            
            // Call Network::contract
            int nContracted=0;
            for (std::set<std::pair<double,std::pair<size_t,size_t> > >::const_iterator smallIter=toBeContracted.begin(); smallIter!=toBeContracted.end(); ++smallIter)
            {
                const size_t i(smallIter->second.first);
                const size_t j(smallIter->second.second);
                typename EdgeFinder<LinkType>::isNetworkEdgeType Lij=DN.link(i,j);
                if (Lij.first)
                {
                    //DN.contractWithConstraintCheck(DN.node(i),DN.node(j));
                    nContracted+=DN.contractWithConstraintCheck(DN.node(i),DN.node(j));
                }
            }
            model::cout<<"("<<nContracted<<" contracted"<<std::flush;
            model::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
            
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
    template <typename DislocationNetworkType>
    double DislocationNetworkRemesh<DislocationNetworkType>::Lmin=10.0;
    
    template <typename DislocationNetworkType>
    double DislocationNetworkRemesh<DislocationNetworkType>::Lmax=100.0;
    
    template <typename DislocationNetworkType>
    double DislocationNetworkRemesh<DislocationNetworkType>::thetaDeg=45.0;
    
    template <typename DislocationNetworkType>
    double DislocationNetworkRemesh<DislocationNetworkType>::neighborRadius=0.001;
    
    // Declare Static Data
    template <typename DislocationNetworkType>
    int DislocationNetworkRemesh<DislocationNetworkType>::verboseRemesh=0;

    
} // namespace model
#endif



//        /**********************************************************************/
//        std::pair<bool,const Simplex<dim,dim>*> pointIsInsideMesh(const VectorDimD& P0, const Simplex<dim,dim>* const guess)
//        {/*!\param[in] P0 position vector
//          * \param[in] guess pointer of the Simplex where the search starts
//          * \returns true if P0 is inside the mesh
//          */
//            std::pair<bool,const Simplex<dim,dim>*> temp(true,NULL);
//            if (DN.shared.use_boundary)
//            {
//                temp=DN.shared.mesh.searchWithGuess(P0,guess);
//            }
//            return temp;
//        }


//        /**********************************************************************/
//        unsigned int contractWithCommonNeighborCheck(const NodeType& Ni,
//                                                     const NodeType& Nj,
//                                                     const VectorDimD& P0)
//        {/*! @param[in] Lij the edge i->j
//          * @param[in] P0 the position of the vertex resulting from contracting Lij
//          *
//          * Contracts the edge i->j into a new node located at P0, making sure
//          * that if P0 is occupied by a neighbor of either i or j, then no
//          * overlapping nodes are created.
//          */
//            unsigned int temp(0);
//            const size_t i(Ni.sID); // StaticID of the source node
//            const size_t j(Nj.sID); // StaticID of the sink   node
//
//            // collect all neighbors at P0 (but i and j)
//            std::set<size_t> neighbors;
//            Ni.neighborsAt(P0,neighbors,neighborRadius);
//            Nj.neighborsAt(P0,neighbors,neighborRadius);
//            //            neighbors.erase(i); // make sure
//            //            neighbors.erase(j); // make sure
//
//            if(neighbors.size())
//            {
//
//                // Remove all existing links among neighbors
//                for (std::set<size_t>::const_iterator nIter1=neighbors.begin();nIter1!=neighbors.end();++nIter1)
//                {
//                    for (std::set<size_t>::const_iterator nIter2=nIter1;nIter2!=neighbors.end();++nIter2)
//                    {
//                        if(nIter2!=nIter1)
//                        {
//                            if(DN.link(*nIter1,*nIter2).first || DN.link(*nIter2,*nIter1).first)
//                            {
//                                //std::cout<<"disconnecting "<<*nIter1<<" "<<*nIter2<<std::endl;
//                                DN.template disconnect<true>(*nIter1,*nIter2);
//                            }
//                        }
//                    }
//                }
//
//                // Contract
//                for (std::set<size_t>::const_iterator nIter=neighbors.begin();nIter!=neighbors.end();++nIter)
//                {
//                    if(nIter!=neighbors.begin())
//                    {
//                        if(DN.node(*nIter).first)
//                        {
//                            //std::cout<<"contractingSecond "<<*neighbors.begin()<<" "<<*nIter<<std::endl;
//                            DN.contractSecond(*neighbors.begin(),*nIter);
//                            temp++;
//                        }
//                    }
//                }
//
//                if(DN.node(i).first && *neighbors.begin()!=i)
//                {
//                    //std::cout<<"contractingSecond "<<*neighbors.begin()<<" "<<i<<std::endl;
//
//                    DN.contractSecond(*neighbors.begin(),i);
//                    temp++;
//                }
//
//                if(DN.node(j).first && *neighbors.begin()!=j)
//                {
//                    //std::cout<<"contractingSecond "<<*neighbors.begin()<<" "<<j<<std::endl;
//
//                    DN.contractSecond(*neighbors.begin(),j);
//                    temp++;
//                }
//            }
//            else // neither i nor j has a neighbor at P0
//            {
//                std::pair<bool,const Simplex<dim,dim>*> guess(DN.pointIsInsideMesh(P0,Ni.includingSimplex()));
//                if(guess.first) // check that P0 is inside mesh
//                {
//                    //std::cout<<"contracting "<<i<<" "<<j<<std::endl;
//
//                    DN.contract(i,j,P0,guess.second);
//                    temp++;
//                }
//            }
//
//            return temp;
//        }
//
//        /**********************************************************************/
//        unsigned int contractSecondWithCommonNeighborCheck(const NodeType& Ni,
//                                                           const NodeType& Nj)
//        {/*!@param[in] i StaticID of the first node (vertex i remains)
//          * @param[in] j StaticID of the second node (vertex j is contracted)
//          *
//          * Contracts  vertex j onto vertex i, making sure no other neighbors of j (but i)
//          * occupies the position of i.
//          */
//            unsigned int temp(0);
//            const int i(Ni.sID);
//            const int j(Nj.sID);
//
//            // collect all neighbors at Ni.get_P() (but j)
//            std::set<size_t> neighbors;
//            Nj.neighborsAt(Ni.get_P(),neighbors,100.0*FLT_EPSILON);
//            //            neighbors.erase(i);
//            //            neighbors.erase(j);
//
//            // Remove all existing links among neighbors
//            for (std::set<size_t>::const_iterator nIter1=neighbors.begin();nIter1!=neighbors.end();++nIter1)
//            {
//                for (std::set<size_t>::const_iterator nIter2=nIter1;nIter2!=neighbors.end();++nIter2)
//                {
//                    if(nIter2!=nIter1)
//                    {
//                        if(DN.link(*nIter1,*nIter2).first || DN.link(*nIter2,*nIter1).first)
//                        {
//                            //std::cout<<"disconnecting "<<*nIter1<<" "<<*nIter2<<std::endl;
//                            DN.template disconnect<true>(*nIter1,*nIter2);
//                        }
//                    }
//                }
//            }
//
//            // Contract
//            neighbors.erase(i);
//
//            //            if(neighbors.size())
//            //            {
//            for (std::set<size_t>::const_iterator nIter=neighbors.begin();nIter!=neighbors.end();++nIter)
//            {
//                if(DN.node(*nIter).first)
//                {
//                    //std::cout<<"contractingSecond "<<i<<" "<<*nIter<<std::endl;
//
//                    DN.contractSecond(i,*nIter);
//                    temp++;
//                }
//            }
//            //            }
//
//            if(DN.node(j).first)
//            {
//                //std::cout<<"contractingSecond (j) "<<i<<" "<<j<<std::endl;
//
//                DN.contractSecond(i,j);
//                temp++;
//            }
//
//            return temp;
//        }
//
//        /**********************************************************************/
//        unsigned int contractSecondWithCommonNeighborCheck(const int& i, const int& j)
//        {/*! @param[in] i StaticID of the first node (vertex i remains)
//          * @param[in] j StaticID of the second node (vertex j is contracted)
//          *
//          * Contracts  vertex j onto vertex i, making sure no other neighbors of j (but i)
//          * occupies the position of i.
//          */
//
//            const typename DislocationNetworkType::isNetworkNodeType Ni(DN.node(i));
//            assert(Ni.first && "NODE i DOES NOT EXIST");
//
//            const typename DislocationNetworkType::isNetworkNodeType Nj(DN.node(j));
//            assert(Nj.first && "NODE j DOES NOT EXIST");
//
//            return contractSecondWithCommonNeighborCheck(*Ni.second,*Nj.second);
//        }



//        /**********************************************************************/
//        unsigned int singleEdgeContract(const typename EdgeFinder<LinkType>::isNetworkEdgeType& Lij)
//        {
//            unsigned int Ncontracted(0);
//            //            if (Lij.first )
//            //            {
//            const size_t i(Lij.second->source->sID);
//            const size_t j(Lij.second->  sink->sID);
//
//
//            const typename DislocationNetworkType::NodeType::VectorOfNormalsType sourcePN(Lij.second->source->constraintNormals());
//            const typename DislocationNetworkType::NodeType::VectorOfNormalsType   sinkPN(Lij.second->  sink->constraintNormals());
//
//            const size_t sourcePNsize(sourcePN.size());
//            const size_t   sinkPNsize(  sinkPN.size());
//
//
//
//            //                if (Lij.second->source->meshLocation()==insideMesh && Lij.second->sink->meshLocation()==insideMesh)
//            //                { // source and sink are inside mesh
//            assert(sourcePNsize>0 && "source->planeNormals() CANNOT HAVE SIZE 0.");
//            assert(  sinkPNsize>0 && "  sink->planeNormals() CANNOT HAVE SIZE 0.");
//            if(sourcePNsize==1 && sinkPNsize==1)
//            {
//                //                        Ncontracted+=contractWithCommonNeighborCheck(Lij,Lij.second->get_r(0.5)); // PATCH FOR COMMON NEIGHBOR AND OUTSIDE-MESH
//                Ncontracted+=contractWithCommonNeighborCheck(*Lij.second->source,*Lij.second->sink,Lij.second->get_r(0.5)); // PATCH FOR COMMON NEIGHBOR AND OUTSIDE-MESH
//            }
//            else if(sourcePNsize==1 && sinkPNsize>1)
//            { // contract source
//                Ncontracted+=contractSecondWithCommonNeighborCheck(j,i);
//            }
//            else if(sourcePNsize>1 && sinkPNsize==1)
//            { // contract sink
//                Ncontracted+=contractSecondWithCommonNeighborCheck(i,j);
//            }
//            else
//            {
//                const VectorDimD P1(Lij.second->source->get_P());
//                const VectorDimD P2(Lij.second-> sink->get_P());
//                const VectorDimD C(P2-P1);
//                const double cNorm(C.norm());
//                if (cNorm<FLT_EPSILON)
//                { // nodes are on top of each other
//                    //                            Ncontracted+=contractWithCommonNeighborCheck(Lij,0.5*(P1+P2)); // PATCH FOR COMMON NEIGHBOR AND OUTSIDE-MESH
//                    Ncontracted+=contractWithCommonNeighborCheck(*Lij.second->source,*Lij.second->sink,0.5*(P1+P2)); // PATCH FOR COMMON NEIGHBOR AND OUTSIDE-MESH
//
//                }
//                else
//                {  // cNorm>=FLT_EPSILON
//                    if(sourcePNsize==2 && sinkPNsize==2)
//                    { // both sink and source move on a line, check if lines intersect
//                        // Compute first direction
//                        VectorDimD d1(sourcePN[0].cross(sourcePN[1]));
//                        const double d1norm(d1.norm());
//                        assert(d1norm>FLT_EPSILON && "DIRECTION d1 HAS ZERO NORM");
//                        d1/=d1norm;
//                        // Compute second direction
//                        VectorDimD d2(sinkPN[0].cross(sinkPN[1]));
//                        const double d2norm(d2.norm());
//                        assert(d2norm>FLT_EPSILON && "DIRECTION d2 HAS ZERO NORM");
//                        d2/=d2norm;
//                        // Necessary condition is plannarity: C*(d1xd2)=0
//                        const VectorDimD d3(d1.cross(d2));
//                        const double d3Norm2(d3.squaredNorm());
//                        if (d3Norm2<FLT_EPSILON)
//                        { // colinear or parallel
//                            if (d1.cross(C/cNorm).norm()<FLT_EPSILON)
//                            { // colinear
//                                //                                        Ncontracted+=contractWithCommonNeighborCheck(Lij,0.5*(P1+P2)); // PATCH FOR COMMON NEIGHBOR AND OUTSIDE-MESH
//
//                                Ncontracted+=contractWithCommonNeighborCheck(*Lij.second->source,*Lij.second->sink,0.5*(P1+P2)); // PATCH FOR COMMON NEIGHBOR AND OUTSIDE-MESH
//                            }
//                        }
//                        else
//                        { // coplanar or no-intersection
//                            if (std::fabs((C/cNorm).dot(d3))<FLT_EPSILON)
//                            { // coplanar
//                                const double u1=C.cross(d2).dot(d3)/d3Norm2;
//                                if(std::fabs(u1<Lmin))
//                                {
//                                    //                                            Ncontracted+=contractWithCommonNeighborCheck(Lij,P1+d1*u1); // PATCH FOR COMMON NEIGHBOR AND OUTSIDE-MESH
//                                    //const VectorDimD P0(P1+d1*u1);
//                                    if(std::fabs(u1)>FLT_EPSILON && std::fabs(std::fabs(u1)-cNorm)>FLT_EPSILON) // intersection point is neither P1 nor P2
//                                    {
//                                        Ncontracted+=contractWithCommonNeighborCheck(*Lij.second->source,*Lij.second->sink,P1+d1*u1);
//                                    }
//                                    else if(std::fabs(u1)>FLT_EPSILON && std::fabs(std::fabs(u1)-cNorm)<=FLT_EPSILON) // intersection point is P2
//                                    {
//                                        Ncontracted+=contractSecondWithCommonNeighborCheck(j,i);
//                                    }
//                                    else if(std::fabs(u1)<=FLT_EPSILON && std::fabs(std::fabs(u1)-cNorm)>FLT_EPSILON) // intersection point is P1
//                                    {
//                                        Ncontracted+=contractSecondWithCommonNeighborCheck(i,j);
//                                    }
//                                    else
//                                    {
//                                        assert(0 && "OOPS! THIS SHOULD NOT HAPPEN");
//                                    }
//                                }
//                            }
//                        }
//                    }
//                    else if(sourcePNsize==2 && sinkPNsize>2)
//                    { // source moves on a line and sink is fixed
//                        VectorDimD d1(sourcePN[0].cross(sourcePN[1]));
//                        double d1norm(d1.norm());
//                        assert(d1norm>FLT_EPSILON && "DIRECTION d1 HAS ZERO NORM");
//                        d1/=d1norm;
//                        if(d1.cross(C/cNorm).norm()<FLT_EPSILON)
//                        {
//                            Ncontracted+=contractSecondWithCommonNeighborCheck(j,i);
//                        }
//                    }
//                    else if(sourcePNsize>2 && sinkPNsize==2)
//                    { // source is fixed and sink moves on a line
//                        VectorDimD d2(sinkPN[0].cross(sinkPN[1]));
//                        double d2norm(d2.norm());
//                        assert(d2norm>FLT_EPSILON && "DIRECTION d2 HAS ZERO NORM");
//                        d2/=d2norm;
//                        if(d2.cross(C/cNorm).norm()<FLT_EPSILON)
//                        {
//                            Ncontracted+=contractSecondWithCommonNeighborCheck(i,j);
//                        }
//                    }
//                    else
//                    {
//                        // both are fixed, cannot contract
//                    }
//                } // end P1==P2
//            } // end case sourcePNsize>1 and sourcePNsize>1
//            //                } // end source and sink are inside mesh
//            //                else if (Lij.second->source->meshLocation()==insideMesh && Lij.second->sink->meshLocation()!=insideMesh)
//            //                { // source is inside mesh, sink in not
//            //                    switch (sourcePNsize)
//            //                    { // decide depending on size of source->planeNormals
//            //                        case 0:
//            //                            assert(0 && "source->planeNormals() CANNOT HAVE SIZE 0.");
//            //                            break;
//            //                        case 1:
//            //                            Ncontracted+=contractSecondWithCommonNeighborCheck(j,i);
//            //                            break;
//            //                        case 2: // source moves on a line
//            //                            // compute intersection of the line and the mesh and place the new node there
//            //                            break;
//            //                        default:
//            //                            // don't do anythig
//            //                            break;
//            //                    }
//            //                }
//            //                else if (Lij.second->source->meshLocation()!=insideMesh && Lij.second->sink->meshLocation()==insideMesh)
//            //                { // sink is inside mesh, source in not
//            //                    switch (sinkPNsize) { // decide depending on size of sink->planeNormals
//            //                        case 0:
//            //                            assert(0 && "sink->planeNormals() CANNOT HAVE SIZE 0.");
//            //                            break;
//            //                        case 1:
//            //                            Ncontracted+=contractSecondWithCommonNeighborCheck(i,j);
//            //                            break;
//            //                        case 2: // sink moves on a line
//            //                            // compute intersection of the line and the mesh and place the new node there
//            //                            break;
//            //
//            //                        default:
//            //                            // don't do anythig
//            //                            break;
//            //                    }
//            //                }
//            //                else{
//            //                    // both are on the mesh, don't do anything
//            //                }
//
//
//            //}
//            return Ncontracted;
//        }
