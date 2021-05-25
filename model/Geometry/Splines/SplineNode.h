/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SplineNode_H_
#define model_SplineNode_H_

#include <Eigen/Dense>
//#include <LoopNode.h>
#include <CatmullRom.h>
#include <Hermite.h>

namespace model
{
    
//    double chordal    =1.0;
//    double centripetal=0.5;
//    double uniform    =0.0;
    
    
    /**************************************************************************/
    /**************************************************************************/
    template <typename Derived,	short unsigned int dim, short unsigned int corder,typename TangentRule>
    class SplineNode
    {
        
    };
    
    /**************************************************************************/
    /* Template specialization for corder=0                                   */
    /**************************************************************************/
    template <typename Derived, short unsigned int dim>
    class SplineNode<Derived, dim,0,Hermite> //: public LoopNode<Derived>
    {
        
    public:
        constexpr static int corder=0;
        constexpr static int NdofXnode=dim;
        
        typedef LoopNode<Derived> LoopNodeType;
        typedef typename TypeTraits<Derived>::LoopNetworkType LoopNetworkType;
        typedef Eigen::Matrix<double, dim, 1>   VectorDim;
        typedef Eigen::Matrix<double, dim, dim> MatrixDim;
        
    private:
        
        VectorDim P;
        
        
    public:
        
        /*************************************************/
        SplineNode(const VectorDim& P_in) :
        /* init list */ P(P_in)
        {
            
        }
        
        /*************************************************/
        void set_P(const VectorDim& P_in)
        {
            P=P_in;
            static_cast<Derived*>(this)->updateGeometry();
//            for(const auto& neighboor : this->neighbors())
//            {
//                std::get<1>(neighboor.second)->updateGeometry();
//            }
            
//            for(const auto& loop : this->loops())
//            {
//                loop->updateGeometry();
//            }
        }
        
        /*************************************************/
        const VectorDim& get_P() const
        {
            return P;
        }
        
    };
    
    
    
//    /**************************************************************************/
//    /* Template specialization for corder=1                                   */
//    /**************************************************************************/
//    template <typename Derived, short unsigned int dim,typename TangentRule>
//    class SplineNode<Derived, dim,1,TangentRule> :
//    /*                              */ public SplineNode<Derived, dim,0,Hermite>
//    {
//
//    public:
//
//        typedef LoopNode<Derived> LoopNodeType;
//        typedef typename LoopNode<Derived>::LinkByLoopContainerType LinkByLoopContainerType;
//        typedef typename LoopNode<Derived>::LoopLinkType LoopLinkType;
//        typedef typename TypeTraits<Derived>::LoopNetworkType LoopNetworkType;
//
//
//        constexpr static int corder=1;
//        constexpr static int NdofXnode=dim;
//
//        typedef SplineNode<Derived, dim,0,Hermite> Base;
//        typedef typename Base::VectorDim VectorDim;
//        typedef typename Base::MatrixDim MatrixDim;
//
//    private:
//
//
//        /**********************************************************************/
//        void computeLoopTangents()
//        {
//
//            loopTangents.clear();
//            const LinkByLoopContainerType temp=this->linksByLoopID();
//            for(const auto& pair : temp)
//            {
//                loopTangents.emplace(pair.first,this->prjM*TangentRule::loopTangent(pair.second));
//            }
//        }
//
//        std::map<size_t,VectorDim> loopTangents;
//
//    public:
//
//        SplineNode(LoopNetworkType* const ln,
//                   const VectorDim& P_in) :
//        /* init list */ Base(ln,P_in)
//        {
//
//        }
//
//        /**********************************************************************/
//        void addLoopLink(LoopLinkType* const pL)
//        {/*@param[in] pL LoopLink pointer
//          *
//          * This functin overrides LoopNode::addLoopLink
//          */
//            LoopNodeType::addLoopLink(pL); // forward to base class
//            computeLoopTangents(); // update tangents due to new connectivity
//        }
//
//        /**********************************************************************/
//        void removeLoopLink(LoopLinkType* const pL)
//        {/*@param[in] pL LoopLink pointer
//          *
//          * This functin overrides LoopNode::removeLoopLink
//          */
//            LoopNodeType::removeLoopLink(pL); // forward to base class
//            computeLoopTangents(); // update tangents due to new connectivity
//        }
//
//        /*************************************************/
//        void set_P(const VectorDim& P_in)
//        {/*@param[in] P_in new position
//          *
//          * This functin overrides SplineNode<Derived, dim,0>::set_P
//          */
//
//            Base::set_P(P_in); // forward to base class
//
//            // Tangents of neighbors may change due to change in P, so updated them
//            for(const auto& neighboor : this->neighbors())
//            {
//                std::get<0>(neighboor.second)->computeLoopTangents();
//            }
//        }
//
//        /**********************************************************************/
//        const std::map<size_t,VectorDim>& tangents() const
//        {
//            return loopTangents;
//        }
//
//        /**********************************************************************/
//        std::map<size_t,std::map<size_t,std::pair<double,VectorDim>>> loopTangentCoeffs() const
//        {
//
//            std::map<size_t,std::map<size_t,std::pair<double,VectorDim>>> temp;
//
//            std::map<size_t,std::set<LoopLinkType*>> linksbyID=this->linksByLoopID();
//            for(const auto& pair : linksbyID)
//            {
//                temp.emplace(pair.first,TangentRule::loopTangentCoeffs(pair.second));
//            }
//
//            return temp;
//        }
//
//
//        /*************************************************/
//        Eigen::Matrix<int,dim,1> node_dofID() const
//        {
//            /*! The IDs of DOFs of this node in the subnetwork
//             */
//            return ( (Eigen::Array<int,dim,1>() << 0, 1, 2).finished()+this->snID()*dim).matrix();
//        }
//    };
    
    
}
#endif

