/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2011 by Benjamin Ramirez<ramirezbrf@gmail.com>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_CatmullRomNode_H_
#define model_CatmullRomNode_H_

#include <model/LoopNetwork/LoopNode.h>
#include <model/Geometry/Splines/SplineNodeBase.h>


namespace model
{
    
    /**************************************************************************/
    /* class template SplineNodeBase: C1-continuous nodes                     */
    /**************************************************************************/
    template <typename Derived,	short unsigned int dim>
    class CatmullRomNode :
    /* inherits */  public LoopNode<Derived>,
    /* inherits */  public SplineNodeBase<Derived,dim,1>
    
    {
        
        typedef LoopNode<Derived> LoopNodeType;
        typedef typename LoopNode<Derived>::LinkByLoopContainerType LinkByLoopContainerType;
        typedef typename LoopNode<Derived>::LoopLinkType LoopLinkType;
        
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        
        
        std::map<size_t,VectorDim> loopTangents;
        
        constexpr static int corder=1;
        typedef SplineNodeBase<Derived,dim,1> SplineNodeBaseType;
        
        /**********************************************************************/
        void computeLoopTangents()
        {
            std::cout<<"Node "<<this->sID<<" computeLoopTangents"<<std::endl;
            
            loopTangents.clear();
            const LinkByLoopContainerType temp=this->linksByLoopID();
            for(const auto& pair : temp)
            {
                loopTangents.emplace(pair.first,VectorDim::Zero());
                
                
                double cT=0.0;
                for(const auto& link : pair.second)
                {
                    const double cL=link->pLink->parametricChordLength();
                    cT+=cL;
                }
                
                for(const auto& link : pair.second)
                {
                    const double cL=link->pLink->parametricChordLength();
                    loopTangents[pair.first]+=(link->sink()->get_P()-link->source()->get_P())/cL*(cT-cL)/cT;
                }
                
            }
        }
        
    public:
        
//        static double alpha; // parametrization exponent

        
        /**********************************************************************/
        CatmullRomNode(const VectorDim& P_in) :
        /* init list */ SplineNodeBaseType(P_in,VectorDim::Zero())
        //        /* init list */ P(roundP(P_in))
        {
            
        }
        
        /**********************************************************************/
        void addLoopLink(LoopLinkType* const pL)
        {
            LoopNodeType::addLoopLink(pL); // forward to base class
            computeLoopTangents();
            this->set_T(loopTangents[0]);

        }
        
        /**********************************************************************/
        void removeLoopLink(LoopLinkType* const pL)
        {
            LoopNodeType::removeLoopLink(pL); // forward to base class
            computeLoopTangents();
            this->set_T(loopTangents[0]);
        }
        
        /**********************************************************************/
        const std::map<size_t,VectorDim>& tangents() const
        {
            return loopTangents;
        }
        
        
    };
    
//    template <typename Derived,	short unsigned int dim>
//    double CatmullRomNode<Derived,dim>::alpha=0.5;

    
}

//#include "model/Geometry/Splines/SplineNode_Hermite.h"
//#include "model/Geometry/Splines/SplineNodeCatmullRom.h"

#endif
