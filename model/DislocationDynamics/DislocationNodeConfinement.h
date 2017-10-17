/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationNodeConfinement_H_
#define model_DislocationNodeConfinement_H_

#include <memory>
#include <set>

#include <model/Utilities/TypeTraits.h>
#include <model/DislocationDynamics/GlidePlanes/GlidePlane.h>
#include <model/DislocationDynamics/BoundingLineSegments.h>


namespace model
{
    template <typename NodeType>
    struct DislocationNodeConfinement :
    /*         */ private std::set<const GlidePlane<typename TypeTraits<NodeType>::LoopType>*>
    {

        static constexpr int dim=TypeTraits<NodeType>::dim;
        typedef typename TypeTraits<NodeType>::LoopType LoopType;
        typedef GlidePlane<LoopType> GlidePlaneType;
        typedef std::set<const GlidePlaneType*> GlidePlaneContainerType;
        
        BoundingLineSegments<dim> _boundingBoxSegments;

        
        /**********************************************************************/
        const NodeType* const node;
        
        /**********************************************************************/
        DislocationNodeConfinement(const NodeType* const pN) :
        /* init */node(pN)
        {
        
        }
        
        /**********************************************************************/
        void clear()
        {
            glidePlanes().clear();
        }
        
        /**********************************************************************/
        bool addGlidePlane(const GlidePlaneType& gp)
        {
            const bool success=glidePlanes().insert(&gp).second;
            if(success)
            {
                assert(gp.contains(node->get_P()) && "Glide Plane does not contain DislocationNode");
//                _isGlissile*=pL->loop()->isGlissile;
                _boundingBoxSegments.updateWithGlidePlane(gp); // Update _boundingBoxSegments. This must be called before updateGlidePlaneIntersections
//                updateGlidePlaneIntersections(pL->loop()->glidePlane);
//                grainSet.insert(&(pL->loop()->grain)); // Insert new grain in grainSet
//                if(grainSet.size()>1)
//                {
//                    std::cout<<"WARNING: CHECK THAT NODE IS ON REGION BND"<<std::endl;
//                }
//                
//                const VectorDim bbP(_boundingBoxSegments.snap(this->get_P()));
//                if((this->get_P()-bbP).squaredNorm()<FLT_EPSILON)
//                {
//                    set_P(bbP);
//                    _isOnBoundingBox=true;
//                }
            }
            return success;
        }
        
        /**********************************************************************/
        const GlidePlaneContainerType& glidePlanes() const
        {
            return *this;
        }
        
        /**********************************************************************/
        GlidePlaneContainerType& glidePlanes()
        {
            return *this;
        }
        
        /**********************************************************************/
        const BoundingLineSegments<dim>& boundingBoxSegments() const
        {
            return _boundingBoxSegments;
        }
        
    };
}

#endif
