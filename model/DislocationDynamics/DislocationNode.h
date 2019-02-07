/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po         <gpo@ucla.edu>.
 * Copyright (C) 2011 by Benjamin Ramirez   <ramirezbrf@gmail.com>.
 * Copyright (C) 2011 by Tamer Crsoby       <tamercrosby@gmail.com>,
 * Copyright (C) 2011 by Can Erel           <canerel55@gmail.com>,
 * Copyright (C) 2011 by Mamdouh Mohamed    <msm07d@fsu.edu>
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DISLOCATIONNODE_H_
#define model_DISLOCATIONNODE_H_

#include <PlanarDislocationNode.h>


//#ifndef NDEBUG
//#define VerboseDislocationNode(N,x) if(verboseDislocationNode>=N){std::cout<<x;}
//#else
//#define VerboseDislocationNode(N,x)
//#endif

namespace model
{
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    class DislocationNode : public PlanarDislocationNode<DislocationNode<dim,corder,InterpolationType>,InterpolationType>
    {

    public:

        
        typedef DislocationNode<dim,corder,InterpolationType> NodeType;
        typedef DislocationSegment<dim,corder,InterpolationType> LinkType;
        typedef PlanarDislocationNode<NodeType,InterpolationType> NodeBaseType;
        typedef typename TypeTraits<NodeType>::LoopType LoopType;
        typedef typename TypeTraits<NodeType>::LoopNetworkType LoopNetworkType;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        constexpr static int NdofXnode=NodeBaseType::NdofXnode;
        typedef Eigen::Matrix<double,NdofXnode,1> VectorDofType;

        
        
        
        /**********************************************************************/
        DislocationNode(LoopNetworkType* const ln,
                        const VectorDim& Pin,
                        const VectorDofType& Vin,
                        const double& vrc) :
        /* base constructor */ NodeBaseType(ln,Pin,Vin,vrc)
        {/*! Constructor from DOF
          */
        }
        
        /**********************************************************************/
        DislocationNode(const LinkType& pL,
                        const double& u) :
        /* base constructor */ NodeBaseType(pL,u)
        {/*! Constructor from ExpandingEdge and DOF
          */
        }
        
        /**********************************************************************/
        DislocationNode(LoopNetworkType* const ln,
                        const VectorDim& Pin,
                        const NodeType* const master) :
        /* base constructor */ NodeBaseType(ln,Pin,master)
        {/*! Constructor from DOF
          */
        }
        
        
        /**********************************************************************/
        template <class T>
        friend T& operator << (T& os, const NodeType& ds)
        {
            os<< DislocationNodeIO<dim>(ds);
            return os;
        }
        
    };
    

    
}
#endif
