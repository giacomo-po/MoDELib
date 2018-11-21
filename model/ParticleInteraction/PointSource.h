/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2012 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _model_PointSource_h
#define _model_PointSource_h

#include <SpatialCellParticle.h>
#include <SingleSourcePoint.h>


namespace model
{

    
    /*!\brief A class template that...
     */
    template<typename Derived,
    /*    */ short unsigned int _dim,
    /*    */ typename ...FieldTypes>
    struct PointSource;
    
    // template specialization: one or more fields
    template<typename Derived,
    /*    */ short unsigned int _dim,
    /*    */ typename FieldType,
    /*    */ typename ...MoreFieldTypes>
    struct PointSource<Derived,_dim,FieldType,MoreFieldTypes...> :
    /* inheritance  */ public PointSource<Derived,_dim,MoreFieldTypes...>,
    /* inheritance  */ public SingleSourcePoint<Derived,FieldType>
    {
        
        typedef typename PointSource<Derived,_dim,MoreFieldTypes...>::VectorDimD VectorDimD;
        
        /**********************************************************************/
        template<typename...T>
        PointSource(const VectorDimD& Pin, const bool& enb, const T&...moreEnb) :
        /* base init */ PointSource<Derived,_dim,MoreFieldTypes...>(Pin,moreEnb...),
        /* base init */ SingleSourcePoint<Derived,FieldType>(enb)
        {/*! @param[in] Pin the position of this PointSource
          * Constructor initializes SpatialCellParticle
          */
        }
        
    };
    
    // template specialization: zero fields
    template<typename Derived,
    /*    */ short unsigned int _dim>
    struct PointSource<Derived,_dim> :
    /* inheritance  */ public SpatialCellParticle<Derived,_dim>
    {
        
        typedef SpatialCellParticle<Derived,_dim> SpatialCellParticleType;
        typedef typename SpatialCellParticleType::VectorDimD VectorDimD;
        
        /**********************************************************************/
        PointSource(const VectorDimD& Pin) :
        /* base init */ SpatialCellParticleType(Pin)
        {/*! @param[in] Pin the position of this PointSource
          * Constructor initializes SpatialCellParticle
          */
        }
        
    };
    
    
} // end namespace
#endif
