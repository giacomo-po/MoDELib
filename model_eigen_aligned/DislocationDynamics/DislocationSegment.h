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

#ifndef model_DISLOCATIONSEGMENT_H
#define model_DISLOCATIONSEGMENT_H

#include <PlanarDislocationSegment.h>

#ifndef NDEBUG
#define VerboseDislocationSegment(N,x) if(verboseDislocationSegment>=N){model::cout<<x;}
#else
#define VerboseDislocationSegment(N,x)
#endif

namespace model
{
    template <int dim, short unsigned int corder, typename InterpolationType>
    struct DislocationSegment : public PlanarDislocationSegment<DislocationSegment<dim,corder,InterpolationType>>
    {
        typedef typename TypeTraits<DislocationSegment<dim,corder,InterpolationType>>::NodeType NodeType;
        typedef PlanarDislocationSegment<DislocationSegment<dim,corder,InterpolationType>> PlanarDislocationSegmentType;
        static constexpr int Ncoeff=PlanarDislocationSegmentType::Ncoeff;
        static constexpr int pOrder=PlanarDislocationSegmentType::pOrder;
        typedef typename PlanarDislocationSegmentType::VectorDim VectorDim;
        typedef typename PlanarDislocationSegmentType::MatrixDim MatrixDim;
        typedef DislocationParticle<dim> DislocationParticleType;
        
#ifdef _MODEL_GREATWHITE_
#include <DislocationSegmentGreatWhite.h>
#endif

        /**********************************************************************/
        DislocationSegment(const std::shared_ptr<NodeType>& nI,
                           const std::shared_ptr<NodeType>& nJ) :
        /* init */ PlanarDislocationSegmentType(nI,nJ)
        {

        }

        /**********************************************************************/
        const MatrixDim& midPointStress() const __attribute__ ((deprecated))
        {/*!\returns The stress matrix for the centre point over this segment.*/
            return this->quadraturePoints().size()? quadraturePoint(this->quadraturePoints().size()/2).stress : MatrixDim::Zero();            
        }
    };
}
#endif
