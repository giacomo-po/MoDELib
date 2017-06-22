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
#include <model/LoopNetwork/LoopNode.h>
#include <model/Geometry/Splines/CatmullRom.h>

//#include <model/Math/RoundEigen.h>
//#include <model/DislocationDynamics/Materials/CrystalOrientation.h>

namespace model
{
    
    double chordal    =1.0;
    double centripetal=0.5;
    double uniform    =0.0;

    
    /**************************************************************************/
    /**************************************************************************/
    template <typename Derived,	short unsigned int dim, short unsigned int corder>
    class SplineNode
    {
        
    };
    
    /**************************************************************************/
    /* Template specialization for corder=0                                   */
    /**************************************************************************/
    template <typename Derived, short unsigned int dim>
    class SplineNode<Derived, dim,0> : public LoopNode<Derived>
    {
        
    public:
        constexpr static int corder=0;
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
//            P= roundP(P_in);
        }
        
        /*************************************************/
        const VectorDim& get_P() const
        {
            return P;
        }
        
    };
    
    /**************************************************************************/
    /* Template specialization for corder=1                                   */
    /**************************************************************************/
    template <typename Derived, short unsigned int dim>
    class SplineNode<Derived, dim,1> :
    /*                              */ public SplineNode<Derived, dim,0>
    {
        
    public:
        
        typedef LoopNode<Derived> LoopNodeType;
        typedef typename LoopNode<Derived>::LinkByLoopContainerType LinkByLoopContainerType;
        typedef typename LoopNode<Derived>::LoopLinkType LoopLinkType;

        
        constexpr static int corder=1;
        typedef SplineNode<Derived, dim,0> Base;
        typedef typename Base::VectorDim VectorDim;
        typedef typename Base::MatrixDim MatrixDim;
        
    private:

        
        /**********************************************************************/
        void computeLoopTangents()
        {
            std::cout<<"Node "<<this->sID<<" computeLoopTangents"<<std::endl;
            
            loopTangents.clear();
            const LinkByLoopContainerType temp=this->linksByLoopID();
            for(const auto& pair : temp)
            {

                loopTangents.emplace(pair.first,CatmullRom<dim>::loopTangent(pair));
                
//                loopTangents.emplace(pair.first,VectorDim::Zero());
                
//                double cT=0.0;
//                for(const auto& link : pair.second)
//                {
//                    const double cL=link->pLink->parametricChordLength();
//                    cT+=cL;
//                }
//                
//                for(const auto& link : pair.second)
//                {
//                    const double cL=link->pLink->parametricChordLength();
//                    loopTangents[pair.first]+=(link->sink()->get_P()-link->source()->get_P())/cL*(cT-cL)/cT;
//                }
                
            }
        }
        
    public:
        std::map<size_t,VectorDim> loopTangents;
        MatrixDim prjM;	//! the projection matrix. THIS SHOULD BE PRIVATE
        
        SplineNode(const VectorDim& P_in) :
        /* init list */ Base(P_in),
        /* init list */ prjM(MatrixDim::Identity())
        {
            
        }
        
        /**********************************************************************/
        void addLoopLink(LoopLinkType* const pL)
        {/*@param[in] pL LoopLink pointer
          *
          * This functin overrides LoopNode::addLoopLink
          */
            LoopNodeType::addLoopLink(pL); // forward to base class
            computeLoopTangents(); // update tangents due to new connectivity
        }
        
        /**********************************************************************/
        void removeLoopLink(LoopLinkType* const pL)
        {/*@param[in] pL LoopLink pointer
          *
          * This functin overrides LoopNode::removeLoopLink
          */
            LoopNodeType::removeLoopLink(pL); // forward to base class
            computeLoopTangents(); // update tangents due to new connectivity
        }
        
        /*************************************************/
        void set_P(const VectorDim& P_in)
        {/*@param[in] P_in new position
          *
          * This functin overrides SplineNode<Derived, dim,0>::set_P
          */
            Base::set_P(P_in); // forward to base class

            // Tangents of neighbors may change due to change in P, so updated them
            
        }
        
        /**********************************************************************/
        const std::map<size_t,VectorDim>& tangents() const
        {
            return loopTangents;
        }
        
    };
    
//    /**************************************************************************/
//    /**************************************************************************/
//    template <typename Derived, short unsigned int dim>
//    class SplineNode<Derived, dim,2> : public SplineNode<Derived, dim,1>
//    {
//        
//        typedef Eigen::Matrix<double, dim, 1>	VectorDim;
//        typedef Eigen::Matrix<double, dim, dim> MatrixDim;
//        
//        VectorDim K;
//        
//    public:
//        
//        const VectorDim& get_K() const
//        {
//            return K;
//        }
//        
//    };
    
}
#endif

