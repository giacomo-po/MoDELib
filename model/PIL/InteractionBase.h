/* This file is part of PIL, the Particle Interaction Library.
 *
 * Copyright (C) 2012 by Giacomo Po <gpo@ucla.edu>
 * Copyright (C) 2012 by Shao-Ching Huang <sch@ucla.edu>
 * Copyright (C) 2012 by Tajendra Singh <tvsingh@ucla.edu>
 * Copyright (C) 2012 by Tamer Crosby <tcrosby@ucla.edu>
 *
 * PIL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _pil_InteractionBase_h
#define _pil_InteractionBase_h

//#include <pil/MpiDataType.h>
#include <pil/SpatialParticle.h>


namespace pil {
    
    
//    template <typename Derived, int _dim>
//    struct SpatialParticle

    
    //template <typename _ParticleType, typename UserSystemProperties = SystemProperties<> >
    template<typename _DataType, int _DataPerParticle>
    struct InteractionBase
    {
        //enum{DataType=MpiDataType<_DataType>::DataType};
        enum{DataPerParticle=_DataPerParticle};
        
        static std::vector<_DataType> resultVector;

//        /*****************************************/
//        template <typename DerivedParticle,int dim>
//        static ResultType getResult(const SpatialParticle<DerivedParticle,dim>& p){
//            const int rid(p.rID);
//            
//            return (ResultType()<<resultVector[rid*3+0],resultVector[rid*3+1],resultVector[rid*3+2]).finished();
//            
//            //resultVector[rid*3+1]
//            //resultVector[rid*3+2]
//        }
        
    };
    
    
    // declare statica data
    template<typename _DataType, int _DataPerParticle>
    std::vector<_DataType> InteractionBase<_DataType,_DataPerParticle>::resultVector;


    
} // end namespace
#endif