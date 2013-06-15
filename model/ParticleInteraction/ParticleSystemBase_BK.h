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
//#include <boost/ptr_container/ptr_vector.hpp>
//#include <boost/ptr_container/ptr_deque.hpp> // SHOULD USE THIS TO STORE PARTICLES, BUT REQUIRES CHANGING THE PARTITIONING ALGORITHM


#include <model/Utilities/modelMacros.h> // model_checkInput
#include <model/ParticleInteraction/InteractionBase.h>
#include <model/ParticleInteraction/SystemProperties.h>

#include <model/SpaceDecomposition/SpatialCellObserver.h>


namespace model {
    
    
    
    template <typename _ParticleType, typename UserSystemProperties>
    class ParticleSystemBase :
//    /* inheritance          */ public SpatialCellObserver<_ParticleType,_ParticleType::dim>,
    /* inheritance          */  protected boost::ptr_map<size_t,_ParticleType>,
    /* inheritance          */  private UserSystemProperties
    {
        
        //        size_t particleCounter;
        
    public:
        
        typedef _ParticleType ParticleType; // make ParticleType available outside the class
        typedef typename _ParticleType::PositionType PositionType;
        typedef SpatialCellObserver<_ParticleType,_ParticleType::dim> SpatialCellObserverType;
        
        typedef typename SpatialCellObserverType::SpatialCellType SpatialCellType;
        //        typedef typename _ParticleType::SpatialCellType SpatialCellType;
        typedef boost::ptr_map<size_t,_ParticleType> ParticleContainerType;
        
        static bool useCellPartitioner;
        
        /**********************************************************************/
        ~ParticleSystemBase() //:
        //        /* init list */ particleCounter(0)
        {
            
            std::cout<<"~ParticleSystemBase destructor: There are now "<< SpatialCellObserverType::totalCells()<< "cells. ";
//            clearParticles();
            std::cout<<" done"<<std::endl;

        }
        
        
        /**********************************************************************/
        template <typename ...AdditionalConstructorTypes>
        ParticleType* addParticle(const PositionType& p, const AdditionalConstructorTypes&... args)
        {/*! Adds a ParticleType particle with position p and an arbitrary (variadic)
          * number of additional constructor parameters.
          */
            // 1- Generate a new ParticleType objects using std::auto_ptr
            std::auto_ptr<ParticleType> pN (new ParticleType(p,args...) );
            
//            // 2- Obtain the sID (StaticID) of the new particle
            const size_t sID(pN->sID);
//            // 3- Insert the particle in the ParticleSystem using sID as the key. Assert succesful insertion
//            std::pair<typename ParticleContainerType::iterator,bool> pIb(this->insert(sID , pN ));
            
            
            // 2- Insert the particle in the ParticleSystem using sID as the key. Assert succesful insertion
            std::pair<typename ParticleContainerType::iterator,bool> pIb(this->insert(sID , pN ));
            model_removeAssert(pIb.second && "CANNOT INSERT PARTICLE IN PARTICLE SYSTEM.");

            
            // 4- Returns the sID of the particle.
            return pIb.first->second;
        }
        
        /**********************************************************************/
        void clearParticles()
        {/*! Removes all particles from the ParticleSystem
          */
            this->clear();
        }
        
        /*****************************************/
        template <typename Property>
        const typename Property::PropertyType& getProperty() const
        {/*!
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
        
        static void setCellSize(const double& cellSize)
        {
            assert(cellSize>0.0);
            SpatialCellObserverType::cellSize=cellSize;
        }
        
        /**********************************************************************/
        template <typename InteractionType>
        void computeMoment0() const
        {
            
            // TO DO: remove from here and parallelize in ParticlesSystemSerial and ParticleSystemParallel
            
            
            //#ifdef _OPENMP
            InteractionType::momentVector.resize(SpatialCellObserverType::totalCells()*InteractionType::DataPerCellMoment);
            
            //#else
            int kk=0;
            for (typename SpatialCellType::CellMapType::const_iterator cIter =SpatialCellObserverType::cellBegin();
                 /*                                                 */ cIter!=SpatialCellObserverType::cellEnd();
                 /*                                               */ ++cIter)
            {
//                InteractionType::computeMoment0(*cIter->second);
                InteractionType::momentVector[kk]=InteractionType::computeMoment0(*cIter->second);

            }
            //#endif
        }
        
        /**********************************************************************/
        template <typename InteractionType>
        void computeFarInteraction()
        {
            
            double t0(clock());
            std::cout<<"Computing far field"<<std::endl;
            
            // TO DO: remove from here and parallelize in ParticlesSystemSerial and ParticleSystemParallel
            
            
//            //#ifdef _OPENMP
//            
//            //#else
//            //int kk=0;
//            
//            
//            
//            
            for (typename ParticleContainerType::iterator pIter=this->begin();pIter!=this->end();++pIter)
            {
                // loop over neighbor cells
                for (typename SpatialCellType::CellMapType::const_iterator cIter =SpatialCellObserverType::cellBegin();
                     /*                                                 */ cIter!=SpatialCellObserverType::cellEnd();
                     /*                                               */ ++cIter)
                {

                        InteractionType(*pIter->second,*cIter->second);  // the constructor of InteractionType actually computes the binary interaction between *iterI and *iterJ
                }
            }
//            //#endif
            
            
            std::cout<<std::scientific<<" ["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]."<<std::endl;

        }



        
    };
    
    // declare static data
    template <typename _ParticleType, typename UserSystemProperties>
    bool ParticleSystemBase<_ParticleType,UserSystemProperties>::useCellPartitioner=true;
    
    
} // end namespace
#endif


//        /**********************************************************************/
//        template <typename CellProperty>
//        void computeCellProperty() const
//        {
//
//            // TO DO: remove from here and parallelize in ParticlesSystemSerial and ParticleSystemParallel
//
//
//            //#ifdef _OPENMP
//
//
//            //#else
//            for (typename SpatialCellType::CellMapType::const_iterator cIter =SpatialCellObserverType::cellBegin();
//                 /*                                                 */ cIter!=SpatialCellObserverType::cellEnd();
//                 /*                                               */ ++cIter)
//            {
//                CellProperty(*cIter->second);
//            }
//            //#endif
//        }


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
