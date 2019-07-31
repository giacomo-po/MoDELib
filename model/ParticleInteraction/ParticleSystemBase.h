/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2012 by Giacomo Po <gpo@ucla.edu>
 * Copyright (C) 2012 by Shao-Ching Huang <sch@ucla.edu>
 * Copyright (C) 2012 by Tajendra Singh <tvsingh@ucla.edu>
 * Copyright (C) 2012 by Tamer Crosby <tcrosby@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _MODEL_ParticleSystemBase_h_
#define _MODEL_ParticleSystemBase_h_

#include <deque>
//#include <modelMacros.h> // model_checkInput
#include <FieldPoint.h>
#include <PointSource.h>
#include <SpatialCellObserver.h>


namespace model
{
    
    
    
    template <typename _ParticleType>
    class ParticleSystemBase :
    /* inheritance          */  public std::deque<_ParticleType> // iteration and fast insertion are priority
    {
        
    public:
        
        typedef _ParticleType ParticleType; // make ParticleType available outside the class
        typedef typename _ParticleType::PositionType PositionType;
        typedef SpatialCellObserver<_ParticleType,_ParticleType::dim> SpatialCellObserverType;
        typedef typename SpatialCellObserverType::SpatialCellType SpatialCellType;
        typedef std::deque<_ParticleType> ParticleContainerType;
        
        /**********************************************************************/
        template <typename ...AdditionalConstructorTypes>
        ParticleType* addParticle(const PositionType& p, const AdditionalConstructorTypes&... args)
        {/*! Adds a ParticleType particle with position p and an arbitrary (variadic)
          * number of additional constructor parameters.
          *\returns a pointer to the newly created particle
          */
            this->emplace_back(p,args...);
            return &*this->rbegin();
        }
        
        /**********************************************************************/
        void clearParticles()
        {/*! Removes all particles from the ParticleSystem
          */
            this->clear();
        }
        
        /**********************************************************************/
        const ParticleContainerType& particles() const
        {/*\returns the const base container of particles
          */
            return *static_cast<const ParticleContainerType* const>(this);
        }
        
        /**********************************************************************/
        const typename SpatialCellObserverType::CellMapType& cells() const
        {
            return SpatialCellObserverType::cells();
        }
        
        /**********************************************************************/
        static void setCellSize(const double& cellSize)
        {
            assert(cellSize>0.0);
            //            SpatialCellObserverType::cellSize=cellSize;
            SpatialCellObserverType::setCellSize(cellSize);
        }
        
        /**********************************************************************/
        template <typename OtherParticleType, typename FieldType,typename... OtherSourceTypes>
        void computeField(OtherParticleType& part, const OtherSourceTypes&... otherSources) const
        {/*!@param[in] part
          * @param[in] otherSources
          * \brief computes the field FieldType for the only particle part,
          * without parallelization
          */
            if(static_cast<FieldPointBase<OtherParticleType,FieldType>* const>(&part)->enabled)
            {
                // Nearest-neighbor interaction
                typename SpatialCellType::CellMapType neighborCells(part.template neighborCells<_ParticleType>());
                
                //! -2 loop over neighbor cells of current particle
                for (typename SpatialCellType::CellMapType::const_iterator cIter =neighborCells.begin();
                     /*                                                 */ cIter!=neighborCells.end();
                     /*                                               */ ++cIter)
                {
                    //! -3 loop over particles in the current neighbor cell
                    for(typename SpatialCellType::ParticleContainerType::const_iterator qIter =cIter->second->particleBegin();
                        /*                                                           */ qIter!=cIter->second->particleEnd();
                        /*                                                         */ ++qIter)
                    {
                        //part.template field<FieldType>() += FieldType::compute(**qIter,part);
                        if(static_cast<const SingleSourcePoint<_ParticleType,FieldType>* const>(*qIter)->enabled)
                        {// source is enabled
                            *static_cast<FieldPointBase<OtherParticleType,FieldType>* const>(&part) += FieldType::compute(**qIter,part);
                        }
                    }
                }
                
                // Non-nearest-neighbor interaction
                if(FieldType::use_multipole)
                {
                    typename SpatialCellType::CellMapType farCells(part.template farCells<_ParticleType>());
                    //
                    //                    part.template field<FieldType>() += FieldType::multipole(part,farCells);
                    *static_cast<FieldPointBase<OtherParticleType,FieldType>* const>(&part) += FieldType::multipole(part,farCells);
                }
                
                // Add contribution of other sources
                if(sizeof...(OtherSourceTypes))
                {
                    *static_cast<FieldPointBase<OtherParticleType,FieldType>* const>(&part) += FieldType::addSourceContribution(part,otherSources...);
                }
            }
            
        }
        
        /**********************************************************************/
        template <class T>
        friend T& operator<< (T& os, const ParticleSystemBase<ParticleType>& system)
        {/*! Operator << uses ParticleType-specific operator <<
          */
            for (const auto& particle : system)
            {
                os << particle <<"\n";
            }
            return os;
        }
        
    };
    
    
    
} // end namespace
#endif

