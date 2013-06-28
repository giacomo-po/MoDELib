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
#include <model/Utilities/modelMacros.h> // model_checkInput
#include <model/ParticleInteraction/FieldPoint.h>
#include <model/ParticleInteraction/PointSource.h>
#include <model/SpaceDecomposition/SpatialCellObserver.h>


namespace model {
    
    
    
    template <typename _ParticleType>
    class ParticleSystemBase :
    /* inheritance          */  protected std::deque<_ParticleType> // iteration and fast insertion are priority
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
            SpatialCellObserverType::cellSize=cellSize;
        }
        
    };
    
    
    
} // end namespace
#endif

