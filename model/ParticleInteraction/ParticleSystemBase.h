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




#include <memory> // for auto_ptr
#include <boost/ptr_container/ptr_map.hpp>

#include <model/Utilities/modelMacros.h> // model_checkInput
#include <model/ParticleInteraction/InteractionBase.h>
#include <model/ParticleInteraction/SystemProperties.h>

#include <model/SpaceDecomposition/SpatialCellObserver.h>


namespace model {
    
    
    // TO DO: REMOVE cellSize from Constructor and simplify constructors
    
    template <typename _ParticleType, typename UserSystemProperties>
    class ParticleSystemBase :
    /* inheritance */  protected boost::ptr_map<size_t,_ParticleType>,
    /* inheritance */  private UserSystemProperties
    {
        
        //        size_t particleCounter;
        
    public:
        
        typedef _ParticleType ParticleType; // make ParticleType available outside the class
        typedef typename _ParticleType::PositionType PositionType;
        typedef typename _ParticleType::SpatialCellType SpatialCellType;
        typedef typename _ParticleType::SpatialCellObserverType SpatialCellObserverType;
        typedef boost::ptr_map<size_t,_ParticleType> ParticleContainerType;
        
        static bool useCellPartitioner;
        
        /**********************************************************************/
        ParticleSystemBase(const double& cellSize=1.0) //:
        //        /* init list */ particleCounter(0)
        {
            SpatialCellObserverType::cellSize=cellSize;
        }
        
        
        /**********************************************************************/
        template <typename ...AdditionalConstructorTypes>
        int addParticle(const PositionType& p, const AdditionalConstructorTypes&... args)
        {/*! Adds a ParticleType particle with position p and an arbitrary (variadic)
          * number of additional constructor parameters.
          */
            // 1- Generate a new ParticleType objects using std::auto_ptr
            std::auto_ptr<ParticleType> pN (new ParticleType(p,args...) );
            // 2- Obtain the sID (StaticID) of the new particle
            const size_t sID(pN->sID);
            
            //            const size_t sID(particleCounter);
            //            particleCounter++;
            
            // 3- Insert the particle in the ParticleSystem using sID as the key. Assert succesful insertion
//            assert(this->insert(sID , pN ).second && "CANNOT INSERT PARTICLE IN PARTICLE SYSTEM.");
            model_execAssert(this->insert(sID , pN ),.second,"CANNOT INSERT PARTICLE IN PARTICLE SYSTEM.");
            // 4- Returns the sID of the particle.
            return sID;
        }
        
        /**********************************************************************/
        void clearParticles()
        {
            this->clear();
        }
        
        /*****************************************/
        template <typename Property>
        const typename Property::PropertyType& getProperty() const
        {/*! Compute full particle interaction according to InteractionType
          */
            return *static_cast<const Property*>(this);
        }
        
        /**********************************************************************/
        const ParticleContainerType& particles() const
        {
            return *static_cast<const ParticleContainerType* const>(this);
        }
        
        /**********************************************************************/
        const typename SpatialCellObserverType::CellMapType& cells() const
        {
            return SpatialCellObserverType::cells();
        }
        
//        /*****************************************/
//        template <class T>
//        friend T& operator<< (T& os, const ParticleContainerType& pS)
//        {/*! Operator << use ParticleType specific << operator
//          */
//            for (typename ParticleContainerType::const_iterator pIter=pS.begin();pIter!=pS.end();++pIter)
//            {
//                os<<(*pIter->second)<<std::endl;
//            }
//            return os;
//        }
        
        //        /*****************************************/
        //        template <class T>
        //        friend T& operator<< (T& os, const typename SpatialCellObserverType::CellMapType& cM)
        //        {/*! Operator << use ParticleType specific << operator
        //          */
        //            for (typename SpatialCellObserverType::CellMapType::const_iterator cIter =cM.cellBegin();
        //                 /*                                                         */ cIter!=cM.cellEnd();
        //                 /*                                                         */ cIter++)
        //            {
        //                os<<(*cIter->second)<<std::endl;
        //            }
        //            return os;
        //        }
        
    };
    
    // declare static data
    template <typename _ParticleType, typename UserSystemProperties>
    bool ParticleSystemBase<_ParticleType,UserSystemProperties>::useCellPartitioner=true;
    
    
} // end namespace
#endif

