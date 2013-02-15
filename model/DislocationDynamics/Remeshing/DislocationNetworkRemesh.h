/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DISLOCATIONNETWORKREMESH_H_
#define model_DISLOCATIONNETWORKREMESH_H_

#include <utility>  // for std::pair and "<" operator between them
#include <set>      // for std::set
#include <assert.h>

#include <Eigen/Dense>
#include <model/Network/Operations/EdgeFinder.h>
#include <model/Math/GramSchmidt.h>
#include <model/BVP/SearchData.h>


namespace model {
	
    /*! \brief Class template that handles the nodal remesh of the DislocationNetwork.
     */
	template <typename DislocationNetworkType>
	class DislocationNetworkRemesh
    {
		
		enum{dim=DislocationNetworkType::dim};
		typedef Eigen::Matrix<double,dim,1> VectorDimD;
		
        typedef typename DislocationNetworkType::LinkType LinkType;
		typedef typename DislocationNetworkType::NetworkLinkContainerType NetworkLinkContainerType;
		
		DislocationNetworkType& DN;
        
        
        /* contractWithCommonNeighborCheck ************************************/
        unsigned int contractWithCommonNeighborCheck(const typename EdgeFinder<LinkType>::isNetworkEdgeType& Lij, const VectorDimD& P0)
        {/*! @param[in] Lij the edge i->j
          * @param[in] P0 the position of the vertex resulting from contracting Lij
          *
          * Contracts the edge i->j into a new node located at P0, making sure
          * that if P0 is occupied by a neighbor of either i or j, then no
          * overlapping nodes are created.
          */
            unsigned int temp(0);
            const size_t i(Lij.second->source->sID); // StaticID of the source node
            const size_t j(Lij.second->  sink->sID); // StaticID of the sink   node
            const std::pair<bool,size_t> isCNi(Lij.second->source->isNeighborAt(P0));
            const std::pair<bool,size_t> isCNj(Lij.second->  sink->isNeighborAt(P0));
            if(isCNi.first && isCNj.first) // both have a neighbor at P0
            {  
                assert(isCNi.second==isCNj.second && "THERE ARE TWO DISTINCT NEIGHBORS AT THE SAME POSITION.");
                DN.contractSecond(isCNi.second,i); 
                temp++;
                DN.contractSecond(isCNj.second,j);
                temp++;
            }
            else if(isCNi.first && !isCNj.first) // only i has a neighbor at P0
            {
                DN.contractSecond(isCNi.second,i);
                temp++;
                if(isCNi.second!=j)
                {
                    DN.contractSecond(isCNi.second,j); 
                    temp++;
                }
            }            
            else if(!isCNi.first && isCNj.first) // only j has a neighbor at P0
            {
                if(isCNj.second!=i)
                {
                    DN.contractSecond(isCNj.second,i); 
                    temp++;
                }
                DN.contractSecond(isCNj.second,j); 
                temp++;
            }
            else // neither i nor j has a neighbor at P0
            {
                if(pointIsInsideMesh(P0,Lij.second->source->meshID())) // check that P0 is inside mesh
                {
                    DN.contract(i,j,P0); 
                    temp++;
                }
            }
            return temp;
        }
        
        /* contractWithCommonNeighborCheck ************************************/
        unsigned int contractSecondWithCommonNeighborCheck(const int& i, const int& j)
        {/*! @param[in] i StaticID of the first node (vertex i remains)
          * @param[in] j StaticID of the second node (vertex j is contracted)
          *
          * Contracts  vertex j onto vertex i, making sure no other neighbors of j (but i)
          * occupies the position of i.
          */
            unsigned int temp(0);

            const typename DislocationNetworkType::isNetworkNodeType Ni(DN.node(i));
            assert(Ni.first && "NODE i DOES NOT EXIST");

            const typename DislocationNetworkType::isNetworkNodeType Nj(DN.node(j));
            assert(Nj.first && "NODE j DOES NOT EXIST");
            
            std::set<size_t> isCNj(Nj.second->areNeighborsAt(Ni.second->get_P()));
            assert(isCNj.erase(i)==1 && "node i must be found at Pi"); // remove i from the set
                        
            for (std::set<size_t>::const_iterator njIter=isCNj.begin(); njIter!=isCNj.end();++njIter)
            {
                const size_t k(*njIter);
                if (DN.node(k).first){
                    DN.contractSecond(i,k); // this could destroy j
                    temp++;
                }
            }
            if (DN.node(j).first) // j still exists
            { 
                DN.contractSecond(i,j); 
                temp++;
            }

            return temp;
        }
        
        
        /************************************************************/
        bool pointIsInsideMesh(const VectorDimD& P0, const size_t& startingMeshID){
            bool temp(true);
            if (DN.shared.boundary_type){
                SearchData<dim> SD(P0);
                DN.shared.domain.findIncludingTet(SD,startingMeshID);
                temp*=(SD.nodeMeshLocation==insideMesh);
            }
            return temp;
        }
        
        
		
	public:
		
        static double Lmax;
        static double Lmin;
        static double thetaDeg;
        
		/* Constructor ********************************************************/
		DislocationNetworkRemesh(DislocationNetworkType& DN_in) :
        /* init list */ DN(DN_in)
        {/*! Initializes the reference to the DislocationNetwork 
          */
        }
		
		/**********************************************************************/
		void remesh()
        {/*! Performs remeshByContraction and then remeshByExpansion. 
          * This order guarantees that 2-vertex subnetworks are expanded.
          */
            remeshByContraction();
			remeshByExpansion();
        }
		
		/* remeshByContraction ************************************************/
		void remeshByContraction()
        {/*! Contract edges according to two criteria.
          */
			
			const double vTolcont=0.0;
			
			std::set<std::pair<double,std::pair<size_t,size_t> > > toBeContracted; // order by increasing segment length
			for (typename NetworkLinkContainerType::const_iterator linkIter=DN.linkBegin();linkIter!=DN.linkEnd();++linkIter){
				const VectorDimD chord(linkIter->second->chord()); // this is sink->get_P() - source->get_P()
				const double chordLength(chord.norm());
				const VectorDimD dv(linkIter->second->sink->get_V()-linkIter->second->source->get_V());				
				bool endsAreApproaching( chord.dot(dv) < -vTolcont*chordLength*dv.norm() );				
                if (endsAreApproaching && chordLength<Lmin){// toBeContracted part                    
					assert(toBeContracted.insert(std::make_pair(chordLength,linkIter->second->nodeIDPair)).second && "COULD NOT INSERT IN SET.");
				}
			}
			
            
			// Call Network::contract 
			unsigned int Ncontracted(0); 
			for (std::set<std::pair<double,std::pair<size_t,size_t> > >::const_iterator smallIter=toBeContracted.begin(); smallIter!=toBeContracted.end(); ++smallIter) {
				const size_t i(smallIter->second.first);
				const size_t j(smallIter->second.second);
				const typename EdgeFinder<LinkType>::isNetworkEdgeType Lij(DN.link(i,j));
                Ncontracted+=singleEdgeContract(Lij);
			}
			std::cout<<" ("<<Ncontracted<<" contracted)"<<std::flush;
			
		}
		
		
		
		/* remeshByExpansion ****************************************************************/
//		void remeshByExpansion(const double& Lmax, const double& Lmin,const double& thetaDeg){
        void remeshByExpansion()
        {
			
			
			double cos_theta_max_crit = std::cos(M_PI-thetaDeg*M_PI/180.0);  /*critical angle */
			std::set<std::pair<size_t,size_t> > toBeExpanded;
			
			// Store the segments to be expanded
			for (typename NetworkLinkContainerType::const_iterator linkIter=DN.linkBegin();linkIter!=DN.linkEnd();++linkIter){
				
				const VectorDimD chord(linkIter->second->chord()); // this is sink->get_P() - source->get_P()
				const double chordLength(chord.norm());
                //				const VectorDimD dv(linkIter->second->sink->get_V()-linkIter->second->source->get_V());				
				
				
				// Always expand single FR source segment
                if (linkIter->second->source->openOrder()==1 && linkIter->second->sink->openOrder()==1)
                {
                    toBeExpanded.insert(linkIter->second->nodeIDPair);
                }
				
                // Expand single FR source segment
                if (   linkIter->second->source->constraintNormals().size()>2
                    && linkIter->second->  sink->constraintNormals().size()>2
                    && chordLength>3.0*Lmin){
                    toBeExpanded.insert(linkIter->second->nodeIDPair);
				}
				
				
				
				if (chordLength>Lmax){
					toBeExpanded.insert(linkIter->second->nodeIDPair);
				}
				
				// Angle criterion
				//double vTolexp=0.5;
				
				if (linkIter->second->source->is_simple()){ //check angle criterion at source
					const VectorDimD c0(linkIter->second->source->openNeighborNode(0)->get_P()-linkIter->second->source->get_P());
					const VectorDimD c1(linkIter->second->source->openNeighborNode(1)->get_P()-linkIter->second->source->get_P());
					const double c0norm(c0.norm());
					const double c1norm(c1.norm());
					if(c0.dot(c1)>cos_theta_max_crit*c0norm*c1norm){
						if (c0norm>3.0*Lmin /*&& c0.dot(v0)>vTolexp*c0norm*v0.norm()*/){
							toBeExpanded.insert(linkIter->second->source->openNeighborLink(0)->nodeIDPair);
						} 
						
						if (c1norm>3.0*Lmin /*&& c1.dot(v1)>vTolexp*c1norm*v1.norm()*/){
							toBeExpanded.insert(linkIter->second->source->openNeighborLink(1)->nodeIDPair);
						}
					}
				}
				if (linkIter->second->sink->is_simple()){ //check angle criterion at sink
					const VectorDimD c0(linkIter->second->sink->openNeighborNode(0)->get_P()-linkIter->second->sink->get_P());
					const VectorDimD c1(linkIter->second->sink->openNeighborNode(1)->get_P()-linkIter->second->sink->get_P());
					const double c0norm(c0.norm());
					const double c1norm(c1.norm());
					if(c0.dot(c1)>cos_theta_max_crit*c0norm*c1norm){
						if (c0norm>3.0*Lmin /*&& c0.dot(v0)>vTolexp*c0norm*v0.norm()*/){
							toBeExpanded.insert(linkIter->second->sink->openNeighborLink(0)->nodeIDPair);
						} 
						if (c1norm>3.0*Lmin/* && c1.dot(v1)>vTolexp*c1norm*v1.norm()*/){
							//														std::cout<<"Expanding 4"<<std::endl;
							toBeExpanded.insert(linkIter->second->sink->openNeighborLink(1)->nodeIDPair);
						}							
					}
				}
				
				if (!linkIter->second->source->is_simple() && !linkIter->second->sink->is_simple() 
					/*&& chord.dot(dv)>vTolexp*chordLength*dv.norm()*/ && chordLength>3.0*Lmin){ // also expands a straight line to generate glissile segment 
					toBeExpanded.insert(linkIter->second->nodeIDPair);
				}
			}
			
			
			// Call Network::expand
			unsigned int Nexpanded(0);
			const double expand_at(0.5);
			for (std::set<std::pair<size_t,size_t> >::const_iterator expIter=toBeExpanded.begin(); expIter!=toBeExpanded.end(); ++expIter) {
				const size_t i(expIter->first);
				const size_t j(expIter->second);	
                const typename EdgeFinder<LinkType>::isNetworkEdgeType Lij(DN.link(i,j));
				if(Lij.first){
                    const VectorDimD expandPoint(Lij.second->get_r(expand_at));
                    //                    bool expandPointInsideMesh(true);
                    //                pointIsInsideMesh(expandPoint,Lij.second->source->meshID());
                    //                    
                    //                    
                    //                    if (DN.shared.boundary_type){
                    //                        SearchData<dim> SD(expandPoint);
                    //                        DN.shared.domain.findIncludingTet(SD,Lij.second->source->meshID());
                    //                        expandPointInsideMesh*=(SD.nodeMeshLocation==insideMesh);
                    //                    }
                    if(pointIsInsideMesh(expandPoint,Lij.second->source->meshID())){
                        DN.expand(i,j,expandPoint);	
                        Nexpanded++;	
                    }
				}
			}
			std::cout<<" ("<<Nexpanded<<" expanded)"<<std::flush;
			
		}
		
		
		
		
		
		/* remeshByContraction **************************************/
		void contract0chordSegments(){
			std::set<std::pair<double,std::pair<size_t,size_t> > > toBeContracted; // order by increasing segment length
			
			for (typename NetworkLinkContainerType::const_iterator linkIter=DN.linkBegin();linkIter!=DN.linkEnd();++linkIter){								
				VectorDimD chord(linkIter->second->chord()); // this is sink->get_P() - source->get_P()
				double chordLength(chord.norm());
				if (chordLength<=FLT_EPSILON){// toBeContracted part
					toBeContracted.insert(std::make_pair(chordLength,linkIter->second->nodeIDPair));
				}
			}
			
			// Call Network::contract 
			for (std::set<std::pair<double,std::pair<size_t,size_t> > >::const_iterator smallIter=toBeContracted.begin(); smallIter!=toBeContracted.end(); ++smallIter) {
				const size_t i(smallIter->second.first);
				const size_t j(smallIter->second.second);
				typename EdgeFinder<LinkType>::isNetworkEdgeType Lij=DN.link(i,j);				
				if (Lij.first ){
//					DN.contractSecond(i,j);	
                    contractSecondWithCommonNeighborCheck(i,j);
				}
			}
		}
        
        
        
        /*************************************************************/
        unsigned int singleEdgeContract(const typename EdgeFinder<LinkType>::isNetworkEdgeType& Lij)
        {
            unsigned int Ncontracted(0);
            if (Lij.first ){
                const size_t i(Lij.second->source->sID);
                const size_t j(Lij.second->  sink->sID);
                //                    const typename DislocationNetworkType::NodeType::VectorOfNormalsType sourcePN(GramSchmidt<dim>(Lij.second->source->planeNormals())); // THIS CAN BE A REFERENCE
                //                    const typename DislocationNetworkType::NodeType::VectorOfNormalsType   sinkPN(GramSchmidt<dim>(Lij.second->  sink->planeNormals())); // THIS CAN BE A REFERENCE
                const typename DislocationNetworkType::NodeType::VectorOfNormalsType sourcePN(GramSchmidt<dim>(Lij.second->source->constraintNormals())); // THIS CAN BE A REFERENCE
                const typename DislocationNetworkType::NodeType::VectorOfNormalsType   sinkPN(GramSchmidt<dim>(Lij.second->  sink->constraintNormals())); // THIS CAN BE A REFERENCE
                
                
                
                
                const size_t sourcePNsize(sourcePN.size());
                const size_t   sinkPNsize(  sinkPN.size());
                if (Lij.second->source->nodeMeshLocation==insideMesh && Lij.second->sink->nodeMeshLocation==insideMesh){ // source and sink are inside mesh
                    //const bool sourceIsBalanced(Lij.second->source->is_balanced());
                    //const bool   sinkIsBalanced(Lij.second->  sink->is_balanced());
                    assert(sourcePNsize>0 && "source->planeNormals() CANNOT HAVE SIZE 0.");
                    assert(  sinkPNsize>0 && "  sink->planeNormals() CANNOT HAVE SIZE 0.");
                    if(sourcePNsize==1 && sinkPNsize==1){
                        //                            std::cout<<"Contract case 1: contracting "<<Lij.second->source->sID<<"->"<<Lij.second->sink->sID<<std::endl;
                        //                            DN.contract(i,j,Lij.second->get_r(0.5));
                        //                            Ncontracted++;
                        Ncontracted+=contractWithCommonNeighborCheck(Lij,Lij.second->get_r(0.5)); // PATCH FOR COMMON NEIGHBOR AND OUTSIDE-MESH
                    }
                    else if(sourcePNsize==1 && sinkPNsize>1){ // contract source
                        //                            std::cout<<"Contract case 2: contracting "<<Lij.second->source->sID<<"->"<<Lij.second->sink->sID<<std::endl;
                        //DN.contractSecond(j,i);
                        //Ncontracted++;
                        Ncontracted+=contractSecondWithCommonNeighborCheck(j,i);
                    }
                    else if(sourcePNsize>1 && sinkPNsize==1){ // contract sink
                        //                            std::cout<<"Contract case 3: contracting "<<Lij.second->source->sID<<"->"<<Lij.second->sink->sID<<std::endl;
                        //                            DN.contractSecond(i,j);
                        //                            Ncontracted++;
                        Ncontracted+=contractSecondWithCommonNeighborCheck(i,j);
                    }
                    else{
                        const VectorDimD P1(Lij.second->source->get_P());
                        const VectorDimD P2(Lij.second-> sink->get_P());
                        const VectorDimD C(P2-P1);
                        const double cNorm(C.norm());
                        if (cNorm<FLT_EPSILON){ // nodes are on top of each other
                            //                                std::cout<<"Contract case 4: contracting "<<Lij.second->source->sID<<"->"<<Lij.second->sink->sID<<std::endl;
                            //                                DN.contract(i,j,0.5*(P1+P2));
                            //                                Ncontracted++;
                            Ncontracted+=contractWithCommonNeighborCheck(Lij,0.5*(P1+P2)); // PATCH FOR COMMON NEIGHBOR AND OUTSIDE-MESH
                        }
                        else{  // cNorm>=FLT_EPSILON
                            if(sourcePNsize==2 && sinkPNsize==2){ // both sink and source move on a line, check if lines intersect
                                // Compute first direction
                                VectorDimD d1(sourcePN[0].cross(sourcePN[1]));
                                const double d1norm(d1.norm());
                                assert(d1norm>FLT_EPSILON && "DIRECTION d1 HAS ZERO NORM");
                                d1/=d1norm;
                                // Compute second direction
                                VectorDimD d2(sinkPN[0].cross(sinkPN[1]));
                                const double d2norm(d2.norm());
                                assert(d2norm>FLT_EPSILON && "DIRECTION d2 HAS ZERO NORM");
                                d2/=d2norm;
                                // Necessary condition is plannarity: C*(d1xd2)=0
                                const VectorDimD d3(d1.cross(d2));
                                const double d3Norm2(d3.squaredNorm());
                                if (d3Norm2<FLT_EPSILON){ // colinear or parallel
                                    if (d1.cross(C/cNorm).norm()<FLT_EPSILON){ // colinear
                                        //                                            std::cout<<"Contract case 5: contracting "<<Lij.second->source->sID<<"->"<<Lij.second->sink->sID<<std::endl;
                                        //DN.contract(i,j,0.5*(P1+P2));
                                        //Ncontracted++;
                                        Ncontracted+=contractWithCommonNeighborCheck(Lij,0.5*(P1+P2)); // PATCH FOR COMMON NEIGHBOR AND OUTSIDE-MESH
                                    }
                                }
                                else{ // coplanar or no-intersection
                                    if (std::fabs((C/cNorm).dot(d3))<FLT_EPSILON){ // coplanar
                                        const double u1=C.cross(d2).dot(d3)/d3Norm2;
                                        if(std::fabs(u1<Lmin)){
                                            //                                                std::cout<<"Contract case 6: contracting "<<Lij.second->source->sID<<"->"<<Lij.second->sink->sID<<std::endl;
                                            //contract(i,j,P1+d1*u1);
                                            //Ncontracted++;
                                            Ncontracted+=contractWithCommonNeighborCheck(Lij,P1+d1*u1); // PATCH FOR COMMON NEIGHBOR AND OUTSIDE-MESH
                                        }
                                    }
                                }
                            }
                            else if(sourcePNsize==2 && sinkPNsize>2){ // source moves on a line and sink is fixed
                                VectorDimD d1(sourcePN[0].cross(sourcePN[1]));
                                double d1norm(d1.norm());
                                assert(d1norm>FLT_EPSILON && "DIRECTION d1 HAS ZERO NORM");
                                d1/=d1norm;
                                if(d1.cross(C/cNorm).norm()<FLT_EPSILON){
                                    //                                        std::cout<<"Contract case 7: contracting "<<Lij.second->source->sID<<"->"<<Lij.second->sink->sID<<std::endl;
                                    //DN.contractSecond(j,i);
                                    //Ncontracted++;
                                    Ncontracted+=contractSecondWithCommonNeighborCheck(j,i);
                                    
                                }
                            }
                            else if(sourcePNsize>2 && sinkPNsize==2){ // source is fixed and sink moves on a line
                                VectorDimD d2(sinkPN[0].cross(sinkPN[1]));
                                double d2norm(d2.norm());
                                assert(d2norm>FLT_EPSILON && "DIRECTION d2 HAS ZERO NORM");
                                d2/=d2norm;
                                if(d2.cross(C/cNorm).norm()<FLT_EPSILON){
                                    //                                        std::cout<<"Contract case 8: contracting "<<Lij.second->source->sID<<"->"<<Lij.second->sink->sID<<std::endl;
                                    //                                        DN.contractSecond(i,j);
                                    //                                        Ncontracted++;
                                    Ncontracted+=contractSecondWithCommonNeighborCheck(i,j);
                                }
                            }
                            else{
                                // both are fixed, cannot contract
                            }
                        } // end P1==P2
                    } // end case sourcePNsize>1 and sourcePNsize>1
                } // end source and sink are inside mesh
                else if (Lij.second->source->nodeMeshLocation==insideMesh && Lij.second->sink->nodeMeshLocation!=insideMesh){ // source is inside mesh, sink in not
                    switch (sourcePNsize) { // decide depending on size of source->planeNormals
                        case 0:
                            assert(0 && "source->planeNormals() CANNOT HAVE SIZE 0.");
                            break;
                        case 1:
                            //                               std::cout<<"Contract case 9: contracting "<<Lij.second->source->sID<<"->"<<Lij.second->sink->sID<<std::endl;
                            //                               DN.contractSecond(j,i);
                            //                               Ncontracted++;
                            Ncontracted+=contractSecondWithCommonNeighborCheck(j,i);
                            break;
                            
                        case 2: // source moves on a line
                            // compute intersection of the line and the mesh and place the new node there
                            break;
                            
                        default:
                            // don't do anythig
                            break;
                    }
                }
                else if (Lij.second->source->nodeMeshLocation!=insideMesh && Lij.second->sink->nodeMeshLocation==insideMesh){ // sink is inside mesh, source in not
                    switch (sinkPNsize) { // decide depending on size of sink->planeNormals
                        case 0:
                            assert(0 && "sink->planeNormals() CANNOT HAVE SIZE 0.");
                            break;
                        case 1:
                            //                               std::cout<<"Contract case 10: contracting "<<Lij.second->source->sID<<"->"<<Lij.second->sink->sID<<std::endl;
                            //                                DN.contractSecond(i,j);
                            //                                Ncontracted++;
                            Ncontracted+=contractSecondWithCommonNeighborCheck(i,j);
                            break;
                            
                        case 2: // sink moves on a line
                            // compute intersection of the line and the mesh and place the new node there
                            break;
                            
                        default:
                            // don't do anythig
                            break;
                    }
                }
                else{
                    // both are on the mesh, don't do anything
                }
                
                
            }
            return Ncontracted;
        }
        
        
        
		
	};
    
    
    // static data
    template <typename DislocationNetworkType>
	double DislocationNetworkRemesh<DislocationNetworkType>::Lmin=10.0;
    
    template <typename DislocationNetworkType>
    double DislocationNetworkRemesh<DislocationNetworkType>::Lmax=100.0;
    
    template <typename DislocationNetworkType>
    double DislocationNetworkRemesh<DislocationNetworkType>::thetaDeg=45.0;
	
	//////////////////////////////////////////////////////////////
} // namespace model
#endif







//				if (Lij.first ){
////                    const typename DislocationNetworkType::NodeType::VectorOfNormalsType sourcePN(GramSchmidt<dim>(Lij.second->source->planeNormals())); // THIS CAN BE A REFERENCE
////                    const typename DislocationNetworkType::NodeType::VectorOfNormalsType   sinkPN(GramSchmidt<dim>(Lij.second->  sink->planeNormals())); // THIS CAN BE A REFERENCE
//                    const typename DislocationNetworkType::NodeType::VectorOfNormalsType sourcePN(GramSchmidt<dim>(Lij.second->source->constraintNormals())); // THIS CAN BE A REFERENCE
//                    const typename DislocationNetworkType::NodeType::VectorOfNormalsType   sinkPN(GramSchmidt<dim>(Lij.second->  sink->constraintNormals())); // THIS CAN BE A REFERENCE
//
//
//
//
//                    const size_t sourcePNsize(sourcePN.size());
//                    const size_t   sinkPNsize(  sinkPN.size());
//                    if (Lij.second->source->nodeMeshLocation==insideMesh && Lij.second->sink->nodeMeshLocation==insideMesh){ // source and sink are inside mesh
//                        //const bool sourceIsBalanced(Lij.second->source->is_balanced());
//                        //const bool   sinkIsBalanced(Lij.second->  sink->is_balanced());
//                        assert(sourcePNsize>0 && "source->planeNormals() CANNOT HAVE SIZE 0.");
//                        assert(  sinkPNsize>0 && "  sink->planeNormals() CANNOT HAVE SIZE 0.");
//                        if(sourcePNsize==1 && sinkPNsize==1){
////                            std::cout<<"Contract case 1: contracting "<<Lij.second->source->sID<<"->"<<Lij.second->sink->sID<<std::endl;
//                            //                            DN.contract(i,j,Lij.second->get_r(0.5));
//                            //                            Ncontracted++;
//                            Ncontracted+=contractWithCommonNeighborCheck(Lij,Lij.second->get_r(0.5)); // PATCH FOR COMMON NEIGHBOR AND OUTSIDE-MESH
//                        }
//                        else if(sourcePNsize==1 && sinkPNsize>1){ // contract source
////                            std::cout<<"Contract case 2: contracting "<<Lij.second->source->sID<<"->"<<Lij.second->sink->sID<<std::endl;
//                            //DN.contractSecond(j,i);
//                            //Ncontracted++;
//                            Ncontracted+=contractSecondWithCommonNeighborCheck(j,i);
//                        }
//                        else if(sourcePNsize>1 && sinkPNsize==1){ // contract sink
////                            std::cout<<"Contract case 3: contracting "<<Lij.second->source->sID<<"->"<<Lij.second->sink->sID<<std::endl;
////                            DN.contractSecond(i,j);
////                            Ncontracted++;
//                            Ncontracted+=contractSecondWithCommonNeighborCheck(i,j);
//                        }
//                        else{
//                            const VectorDimD P1(Lij.second->source->get_P());
//                            const VectorDimD P2(Lij.second-> sink->get_P());
//                            const VectorDimD C(P2-P1);
//                            const double cNorm(C.norm());
//                            if (cNorm<FLT_EPSILON){ // nodes are on top of each other
////                                std::cout<<"Contract case 4: contracting "<<Lij.second->source->sID<<"->"<<Lij.second->sink->sID<<std::endl;
//                                //                                DN.contract(i,j,0.5*(P1+P2));
//                                //                                Ncontracted++;
//                                Ncontracted+=contractWithCommonNeighborCheck(Lij,0.5*(P1+P2)); // PATCH FOR COMMON NEIGHBOR AND OUTSIDE-MESH
//                            }
//                            else{  // cNorm>=FLT_EPSILON
//                                if(sourcePNsize==2 && sinkPNsize==2){ // both sink and source move on a line, check if lines intersect
//                                    // Compute first direction
//                                    VectorDimD d1(sourcePN[0].cross(sourcePN[1]));
//                                    const double d1norm(d1.norm());
//                                    assert(d1norm>FLT_EPSILON && "DIRECTION d1 HAS ZERO NORM");
//                                    d1/=d1norm;
//                                    // Compute second direction
//                                    VectorDimD d2(sinkPN[0].cross(sinkPN[1]));
//                                    const double d2norm(d2.norm());
//                                    assert(d2norm>FLT_EPSILON && "DIRECTION d2 HAS ZERO NORM");
//                                    d2/=d2norm;
//                                    // Necessary condition is plannarity: C*(d1xd2)=0
//                                    const VectorDimD d3(d1.cross(d2));
//                                    const double d3Norm2(d3.squaredNorm());
//                                    if (d3Norm2<FLT_EPSILON){ // colinear or parallel
//                                        if (d1.cross(C/cNorm).norm()<FLT_EPSILON){ // colinear
////                                            std::cout<<"Contract case 5: contracting "<<Lij.second->source->sID<<"->"<<Lij.second->sink->sID<<std::endl;
//                                            //DN.contract(i,j,0.5*(P1+P2));
//                                            //Ncontracted++;
//                                            Ncontracted+=contractWithCommonNeighborCheck(Lij,0.5*(P1+P2)); // PATCH FOR COMMON NEIGHBOR AND OUTSIDE-MESH
//                                        }
//                                    }
//                                    else{ // coplanar or no-intersection
//                                        if (std::fabs((C/cNorm).dot(d3))<FLT_EPSILON){ // coplanar
//                                            const double u1=C.cross(d2).dot(d3)/d3Norm2;
//                                            if(std::fabs(u1<Lmin)){
////                                                std::cout<<"Contract case 6: contracting "<<Lij.second->source->sID<<"->"<<Lij.second->sink->sID<<std::endl;
//                                                //contract(i,j,P1+d1*u1);
//                                                //Ncontracted++;
//                                                Ncontracted+=contractWithCommonNeighborCheck(Lij,P1+d1*u1); // PATCH FOR COMMON NEIGHBOR AND OUTSIDE-MESH
//                                            }
//                                        }
//                                    }
//                                }
//                                else if(sourcePNsize==2 && sinkPNsize>2){ // source moves on a line and sink is fixed
//                                    VectorDimD d1(sourcePN[0].cross(sourcePN[1]));
//                                    double d1norm(d1.norm());
//                                    assert(d1norm>FLT_EPSILON && "DIRECTION d1 HAS ZERO NORM");
//                                    d1/=d1norm;
//                                    if(d1.cross(C/cNorm).norm()<FLT_EPSILON){
////                                        std::cout<<"Contract case 7: contracting "<<Lij.second->source->sID<<"->"<<Lij.second->sink->sID<<std::endl;
//                                        //DN.contractSecond(j,i);
//                                        //Ncontracted++;
//                                        Ncontracted+=contractSecondWithCommonNeighborCheck(j,i);
//
//                                    }
//                                }
//                                else if(sourcePNsize>2 && sinkPNsize==2){ // source is fixed and sink moves on a line
//                                    VectorDimD d2(sinkPN[0].cross(sinkPN[1]));
//                                    double d2norm(d2.norm());
//                                    assert(d2norm>FLT_EPSILON && "DIRECTION d2 HAS ZERO NORM");
//                                    d2/=d2norm;
//                                    if(d2.cross(C/cNorm).norm()<FLT_EPSILON){
////                                        std::cout<<"Contract case 8: contracting "<<Lij.second->source->sID<<"->"<<Lij.second->sink->sID<<std::endl;
////                                        DN.contractSecond(i,j);
////                                        Ncontracted++;
//                                        Ncontracted+=contractSecondWithCommonNeighborCheck(i,j);
//                                    }
//                                }
//                                else{
//                                    // both are fixed, cannot contract
//                                }
//                            } // end P1==P2
//                        } // end case sourcePNsize>1 and sourcePNsize>1
//                    } // end source and sink are inside mesh
//                    else if (Lij.second->source->nodeMeshLocation==insideMesh && Lij.second->sink->nodeMeshLocation!=insideMesh){ // source is inside mesh, sink in not
//                        switch (sourcePNsize) { // decide depending on size of source->planeNormals
//                            case 0:
//                                assert(0 && "source->planeNormals() CANNOT HAVE SIZE 0.");
//                                break;
//                            case 1:
// //                               std::cout<<"Contract case 9: contracting "<<Lij.second->source->sID<<"->"<<Lij.second->sink->sID<<std::endl;
// //                               DN.contractSecond(j,i);
// //                               Ncontracted++;
//                                Ncontracted+=contractSecondWithCommonNeighborCheck(j,i);
//                                break;
//
//                            case 2: // source moves on a line
//                                // compute intersection of the line and the mesh and place the new node there
//                                break;
//
//                            default:
//                                // don't do anythig
//                                break;
//                        }
//                    }
//                    else if (Lij.second->source->nodeMeshLocation!=insideMesh && Lij.second->sink->nodeMeshLocation==insideMesh){ // sink is inside mesh, source in not
//                        switch (sinkPNsize) { // decide depending on size of sink->planeNormals
//                            case 0:
//                                assert(0 && "sink->planeNormals() CANNOT HAVE SIZE 0.");
//                                break;
//                            case 1:
// //                               std::cout<<"Contract case 10: contracting "<<Lij.second->source->sID<<"->"<<Lij.second->sink->sID<<std::endl;
////                                DN.contractSecond(i,j);
////                                Ncontracted++;
//                                Ncontracted+=contractSecondWithCommonNeighborCheck(i,j);
//                                break;
//
//                            case 2: // sink moves on a line
//                                // compute intersection of the line and the mesh and place the new node there
//                                break;
//
//                            default:
//                                // don't do anythig
//                                break;
//                        }
//                    }
//                    else{
//                        // both are on the mesh, don't do anything
//                    }
//
//
//				}

