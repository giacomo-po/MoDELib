/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2011 by Benjamin Ramirez<ramirezbrf@gmail.com>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SPLINENODE_CATMULLROM_H_
#define model_SPLINENODE_CATMULLROM_H_

//#include <vector>
#include <assert.h>
#include <iterator>
#include <Eigen/Dense>
#include <model/Network/NetworkNode.h>
#include <model/Network/Operations/EdgeExpansion.h>
//#include <model/Math/CompileTimeMath/Pow.h>
#include <model/Geometry/Splines/SplineNodeBase.h>


namespace model
{
    

    
    template <typename Derived, short unsigned int dim>
    class SplineNode<Derived, dim,1,CatmullRom> :
    /* inherits */  public NetworkNode<Derived>,
    /* inherits */  public SplineNodeBase<Derived,dim,1>
    {
        
        
        
    public:
        
        //	Eigen::VectorXi edgeConfiguration;
        
        
        enum  {corder=1};
        enum  {NdofXnode=dim};
        
        typedef Eigen::Matrix<double, NdofXnode, 1> VectorDofType;
        
//#include "model/Geometry/Splines/SplineNodeBase_corder1.h"
        typedef SplineNode<Derived, dim,corder,CatmullRom> NodeType;
        
        
        typedef NetworkNode<Derived> NetworkNodeType;
        typedef typename NetworkNodeType::LinkType LinkType;
        typedef typename NetworkNodeType::NeighborType NeighborType;
        typedef typename NetworkNodeType::NeighborContainerType NeighborContainerType;
        
        
        typedef typename NeighborContainerType::iterator NeighborIteratorType;
        typedef typename NeighborContainerType::const_iterator constNeighborIteratorType;
        
        typedef SplineNodeBase<Derived,dim,1> SplineNodeBaseType;
        typedef typename SplineNodeBaseType::VectorDim VectorDim;
        typedef typename SplineNodeBaseType::MatrixDim MatrixDim;
        
        
    public:
        
        SplineNode(const VectorDim& P_in) :
        /* init list */ SplineNodeBaseType(P_in,VectorDim::Zero())
        {
//            std::cout<<"Verify SplineNode_CatmullRom get_nodeDof() 1"<<std::endl;
            set(this->get_P()); // trigger calculation of tangents
        }
        
        SplineNode(const EdgeRef<LinkType>& pL, const double& u) :
        /* init list */ NetworkNode<Derived>::NetworkNode(pL),
        /* init list */ SplineNodeBaseType(pL.E.get_r(u),VectorDim::Zero())
        {
//            std::cout<<"Verify SplineNode_CatmullRom get_nodeDof() 2"<<std::endl;
            set(this->get_P()); // trigger calculation of tangents
        }
        
        SplineNode(const EdgeRef<LinkType>& pL, const VectorDim& P_in) :
        /* init list */ NetworkNode<Derived>::NetworkNode(pL),
        /* init list */ SplineNodeBaseType(P_in,VectorDim::Zero())
        {
//            std::cout<<"SplineNode_CatmullRom from ExpadingLink B "<<this->sID<<std::endl;

//            std::cout<<"Verify SplineNode_CatmullRom get_nodeDof() 3"<<std::endl;
            set(this->get_P()); // trigger calculation of tangents
        }
        
        /*************************************************/
        double neighborLinksLength() const
        {
            double temp=0.0;
            for (constNeighborIteratorType neighborIter=this->neighborhood().begin();neighborIter!=this->neighborhood().end();++neighborIter)
            {
                if (std::get<2>(neighborIter->second)!=0)
                {
                    temp+=std::get<1>(neighborIter->second)->chordLength();
                }
            }
            return temp;
        }
        
        /*************************************************/
        void set(const VectorDim& P_in)
        {
            /*! Because of the CatmullRom rule for parametric tangents, changing the position of this
             *  CatmullRom node must also change its parametric tangent and the
             *  parametric tangents of its first-neighboring nodes. This in turn affects the shape
             *  of all SplineSegments attached to the first-neighboring nodes. The update strategy is as follows:
             */
            
            //! 1- Sets P=P_in
            //this->P=P_in;
            SplineNodeBaseType::set_P(P_in);
            
            //! 2- Transmits 'make_T' on the neighbor nodes of level (corder+1=2), that is this node and its first-neighbors
            typedef void (Derived::*node_member_function_pointer_type)(void);
            node_member_function_pointer_type Nmfp;
            Nmfp=&Derived::make_T;
            this->depthFirstNodeExecute(Nmfp,corder+1); // transmit to 1st neighbors (corder+1=2)
        }
        
        /*************************************************/
        void make_T()
        {
            //! Calls make_CR2H()
            make_CR2H();
            
            //! computes T=(CR2H*VectorDof).segment(dim,dim)
            SplineNodeBaseType::set_T((CR2H*VectorDof).template segment<dim>(dim));
        }
        
        /*************************************************/
        VectorDim get_nodeDof() const
        {
            return this->get_P();
        }
        
        /*************************************************/
        Eigen::Matrix<int,dim,1> node_dofID() const {
            /*! The IDs of DOFs of this node in the subnetwork
             */
            return ( (Eigen::Array<int,dim,1>() << 0, 1, 2).finished()+this->snID()*dim).matrix();
        }
        
        /*************************************************/
        Eigen::VectorXd get_dof() const
        {
            return VectorDof;
        }
        
        /*************************************************/
        Eigen::VectorXi dofID() const
        {
            /*! The IDs of DOFs of this and the neighbor nodes in the subnetwork
             */
            Eigen::VectorXi VectorDofID;
            VectorDofID.resize(Ndof);
            
            size_t k=0;
            for (constNeighborIteratorType neighborIter=this->neighborhood().begin(); neighborIter != this->neighborhood().end(); ++neighborIter) {
                VectorDofID.segment<dim>(k*dim)	=std::get<0>(neighborIter->second)->node_dofID();
                ++k;
            }
            
            return VectorDofID;
        }
        
        /*************************************************/
        const Eigen::Matrix<double, dim*(corder+1), Eigen::Dynamic> & W2H() const
        {
            return CR2H;
        }
        
        
        /*************************************************/
        Eigen::Matrix<double, dim*(corder+1), Eigen::Dynamic> W2Ht() const
        {
            Eigen::Matrix<double, dim*(corder+1), Eigen::Dynamic> temp(CR2H);
            
//            LinkType
//            const double alpha(0.5); //! CHANGE THIS IN CATMULLROM
            const double alpha(LinkType::alpha); //! CHANGE THIS IN CATMULLROM

            
            if ( !this->is_isolated() && this->is_balanced() ){
                //	double CPL;
                //	double CPLDPinv=0.0;
                //	double CPLARinv=0.0;
                double CPLT=0.0;
                double sjT(0.0);
                double sjOverGjT(0.0);
                
                MatrixDim cP02(MatrixDim::Zero());
                VectorDim cP03a(VectorDim::Zero());
                VectorDim cP03b(VectorDim::Zero());
                
                
                //         int sgnID(0);
                for (constNeighborIteratorType neighborIter=this->neighborhood().begin();neighborIter!=this->neighborhood().end();++neighborIter){
                    
                    const int dir(std::get<2>(neighborIter->second));
                    if (dir!=0)
                    {
                        int edgeConfig(0);
                        switch ( dir )
                        {
                            case  1:	// out
                                edgeConfig=std::get<1>(neighborIter->second)->sourceTfactor;
                                break;
                            case -1:    // in
                                edgeConfig=std::get<1>(neighborIter->second)->  sinkTfactor;
                                break;
                            default:
                                assert(0);
                                break;
                        }
                        const double CPL=std::get<1>(neighborIter->second)->chordParametricLength();	// chord parametric length
                        CPLT+=CPL;
                        sjT+=edgeConfig;
                        sjOverGjT+=edgeConfig/CPL;
                        
                        const VectorDim ci( (std::get<0>(neighborIter->second)->get_P()-this->get_P())/CPL );
                        cP02+= alpha * edgeConfig * ci*ci.transpose() / CPL;
                        cP03a+= edgeConfig*ci;
                        cP03b+= ci*alpha/CPL;
                        
                        
                        //                 sgnID++;
                    }
                    
                    
                }
                
                //          assert(sgnID==edgeConfiguration.size());
                
                MatrixDim cP03( - cP03a/CPLT/CPLT * cP03b.transpose());
                
                
                double CPLTinv=1.0/CPLT;
                //Ndof=dim*this->neighborhood().size();
                //		short int dir;
                int k=0;
                //sgnID=0;
                int sgnID=0;
                
                for (constNeighborIteratorType neighborIter=this->neighborhood().begin();neighborIter!=this->neighborhood().end();++neighborIter){
                    int dir=std::get<2>(neighborIter->second);
                    switch ( dir ) {
                        case  0:	// self
                            temp.template block<dim,dim>(0,k*dim).setIdentity();
                            temp.template block<dim,dim>(dim,k*dim)=this->prjM*((sjT/CPLT -  sjOverGjT)*MatrixDim::Identity()
                                                                                +cP02
                                                                                +cP03)/(this->neighborhood().size()-2);
                            break;
                            //                    default:	// departing or arriving
                            //                        const double CPL=std::get<1>(neighborIter->second)->chordParametricLength();	// chord parametric length
                            //                        const VectorDim ci( (std::get<0>(neighborIter->second)->get_P()-this->get_P())/CPL );
                            //
                            //                        temp.template block<dim,dim>(dim,k*dim)= this->prjM*(edgeConfiguration(sgnID)*( 1.0/CPL-CPLTinv)*MatrixDim::Identity()
                            //                                                                             -ci*ci.transpose()*edgeConfiguration(sgnID)*alpha/CPL
                            //                                                                             +cP03a/CPLT/CPLT*ci.transpose()*alpha)/(this->neighborhood().size()-2);
                            //                        sgnID++;
                            //                        break;
                        case  1:	// departing
                        {
                            const double CPL=std::get<1>(neighborIter->second)->chordParametricLength();	// chord parametric length
                            const VectorDim ci( (std::get<0>(neighborIter->second)->get_P()-this->get_P())/CPL );
                            int edgeConfig=std::get<1>(neighborIter->second)->sourceTfactor;
                            
                            temp.template block<dim,dim>(dim,k*dim)= this->prjM*(edgeConfig*( 1.0/CPL-CPLTinv)*MatrixDim::Identity()
                                                                                 -ci*ci.transpose()*edgeConfig*alpha/CPL
                                                                                 +cP03a/CPLT/CPLT*ci.transpose()*alpha)/(this->neighborhood().size()-2);
                            //    sgnID++;
                        }
                            break;
                        case -1:	// arriving
                        {
                            const double CPL=std::get<1>(neighborIter->second)->chordParametricLength();	// chord parametric length
                            const VectorDim ci( (std::get<0>(neighborIter->second)->get_P()-this->get_P())/CPL );
                            int edgeConfig=std::get<1>(neighborIter->second)->sinkTfactor;
                            
                            temp.template block<dim,dim>(dim,k*dim)= this->prjM*(edgeConfig*( 1.0/CPL-CPLTinv)*MatrixDim::Identity()
                                                                                 -ci*ci.transpose()*edgeConfig*alpha/CPL
                                                                                 +cP03a/CPLT/CPLT*ci.transpose()*alpha)/(this->neighborhood().size()-2);
                            //    sgnID++;
                        }
                            break;
                    }
                    ++k;
                }
                
                //            assert(sgnID==edgeConfiguration.size());
                
            }
            return temp;
        }
        
        
        
        
        
    private:
        
        
        //	size_t Nneighbors;
        size_t Ndof;
        
        Eigen::VectorXd VectorDof;
        //	Eigen::VectorXi VectorDofID;
        
        Eigen::Matrix<double, dim*(corder+1), Eigen::Dynamic> CR2H;
        
        
        
        /* findEdgeConfiguration (possibly overwritten by Derived) ************************/
        void findEdgeConfiguration(){
            //
            assert(0 && "NEED TO RE-IMPLEMENT ORIGINAL CR-RULE");
        }
        
        
        
        
        /* make_CRneighbors ************************************/
        void make_CRneighbors()
        {
            //size_t Nneighbors=this->neighborhood().size();
            Ndof=dim*this->neighborhood().size();
            VectorDof.resize(Ndof);
            size_t k=0;
            for (constNeighborIteratorType neighborIter=this->neighborhood().begin(); neighborIter != this->neighborhood().end(); ++neighborIter) {
                VectorDof.segment<dim>(k*dim)	=std::get<0>(neighborIter->second)->get_P();
                ++k;
            }
        }
        
        /*************************************************/
        void make_CR2H(){
            
            //! 1- Call make_CRneighbors()
            make_CRneighbors();
            
            //! 2- Resize CR2H to dim*(corder+1) x Ndof
            CR2H.setZero(dim*(corder+1),Ndof);
            
            
            
            if (this->is_isolated())
            {
                //            std::cout<<"Node "<<this->sID<<" is isolated"<<std::endl;
                CR2H.template block<dim,dim>(0,0).setIdentity();
                CR2H.template block<dim,dim>(dim,0).setZero();
            }
            else
            {
                if (this->is_balanced())
                {
                    //                        std::cout<<"Node "<<this->sID<<" is not isolated and balanced "<<std::endl;
                    make_CR2H_central();
                }
                else
                {
                    //             std::cout<<"Node "<<this->sID<<" is not isolated and not balanced "<<std::endl;
                    
                    int k=0;
                    for (constNeighborIteratorType neighborIter=this->neighborhood().begin();neighborIter!=this->neighborhood().end();++neighborIter)
                    {
                        const short int dir(std::get<2>(neighborIter->second));
                        switch ( dir )
                        {
                            case  0:	// self
                                CR2H.template block<dim,dim>(0,k*dim).setIdentity();
                                //							CR2H.template block<dim,dim>(dim,k*dim)=this->prjM*( (int(this->outOrder())-int(this->inOrder() ))/CPLT - CPLDPinv + CPLARinv)/(CRneighbors.size()-2);
                                break;
                        }
                        ++k;
                    }
                }
            }
            
        }
        
        /*************************************************/
        void make_CR2H_central()
        {
            
            //	double CPL;
            //	double CPLDPinv=0.0;
            //	double CPLARinv=0.0;
            double CPLT=0.0;
            double sjT(0.0);
            double sjOverGjT(0.0);
            
            //		int sgnID(0);
            for (constNeighborIteratorType neighborIter=this->neighborhood().begin();neighborIter!=this->neighborhood().end();++neighborIter)
            {
                const int dir(std::get<2>(neighborIter->second));
                if (dir!=0) // not self
                {
                    int edgeConfig(0);
                    switch ( dir )
                    {
                        case  1:	// out
                            edgeConfig=std::get<1>(neighborIter->second)->sourceTfactor;
                            break;
                        case -1:    // in
                            edgeConfig=std::get<1>(neighborIter->second)->  sinkTfactor;
                            break;
                        default:
                            assert(0);
                            break;
                    }
                    const double CPL(std::get<1>(neighborIter->second)->chordParametricLength());	// chord parametric length
                    CPLT+=CPL;
                    //                assert(edgeConfiguration(sgnID)==edgeConfig);
                    //                sjT+=edgeConfiguration(sgnID);
                    //                sjOverGjT+=edgeConfiguration(sgnID)/CPL;
                    //                sgnID++;
                    
                    //                if (this->sID==88)
                    //                {
                    //                    std::cout<<std::get<1>(neighborIter->second)->source->sID<<std::endl;
                    //                    std::cout<<std::get<1>(neighborIter->second)->sink->sID<<std::endl;
                    //                    std::cout<<"edgeConfig="<<edgeConfig<<std::endl;
                    //                }
                    
                    sjT+=edgeConfig;
                    sjOverGjT+=edgeConfig/CPL;
                }
                
                //			switch ( dir )
                //            {
                //				case  0:	// self
                //					break;
                //				default:	// neighbor
                //					double CPL=std::get<1>(neighborIter->second)->chordParametricLength();	// chord parametric length
                //					CPLT+=CPL;
                //					sjT+=edgeConfiguration(sgnID);
                //					sjOverGjT+=edgeConfiguration(sgnID)/CPL;
                //					sgnID++;
                //					break;
                //			}
            }
            
            double CPLTinv=1.0/CPLT;
            //Ndof=dim*this->neighborhood().size();
            //		short int dir;
            int k=0;
            
            //        int sgnID(0);
            
            //		sgnID=0;
            
            for (constNeighborIteratorType neighborIter=this->neighborhood().begin();neighborIter!=this->neighborhood().end();++neighborIter){
                const int dir(std::get<2>(neighborIter->second));
                switch ( dir ) {
                    case  0:	// self
                        CR2H.template block<dim,dim>(0,k*dim).setIdentity();
                        CR2H.template block<dim,dim>(dim,k*dim)=this->prjM*( sjT/CPLT -  sjOverGjT)/(this->neighborhood().size()-2);						
                        break;
                        //				default:	// departing or arriving
                        //					CR2H.template block<dim,dim>(dim,k*dim)= edgeConfiguration(sgnID)*this->prjM*( 1.0/std::get<1>(neighborIter->second)->chordParametricLength()-CPLTinv)/(this->neighborhood().size()-2);
                        //					sgnID++;
                        //					break;
                    case  1:	// departing or arriving
                        //                    const int edgeConfig=std::get<1>(neighborIter->second)->sourceTfactor;
                        //					CR2H.template block<dim,dim>(dim,k*dim)= edgeConfiguration(sgnID)*this->prjM*( 1.0/std::get<1>(neighborIter->second)->chordParametricLength()-CPLTinv)/(this->neighborhood().size()-2);
                        CR2H.template block<dim,dim>(dim,k*dim)= std::get<1>(neighborIter->second)->sourceTfactor*this->prjM*( 1.0/std::get<1>(neighborIter->second)->chordParametricLength()-CPLTinv)/(this->neighborhood().size()-2);
                        
                        //					sgnID++;
                        break;
                    case -1:	// departing or arriving
                        //                    const int edgeConfig=std::get<1>(neighborIter->second)->sinkTfactor;
                        //					CR2H.template block<dim,dim>(dim,k*dim)= edgeConfiguration(sgnID)*this->prjM*( 1.0/std::get<1>(neighborIter->second)->chordParametricLength()-CPLTinv)/(this->neighborhood().size()-2);
                        CR2H.template block<dim,dim>(dim,k*dim)= std::get<1>(neighborIter->second)->sinkTfactor*this->prjM*( 1.0/std::get<1>(neighborIter->second)->chordParametricLength()-CPLTinv)/(this->neighborhood().size()-2);
                        //					sgnID++;
                        break;
                    default:
                        assert(0);
                        break;
                }
                ++k;
            }
        }
        
    };
    
}

#endif
