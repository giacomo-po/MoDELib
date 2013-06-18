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

#ifndef _MODEL_ParticleSystemSerial_h_
#define _MODEL_ParticleSystemSerial_h_


#ifdef _OPENMP
#include <omp.h>
#define EIGEN_DONT_PARALLELIZE // disable Eigen Internal openmp Parallelization
#else
#define omp_get_thread_num() 0
#define omp_get_num_threads() 1
//#define omp_get_max_threads() 1
#endif


#include <memory> // for auto_ptr
#include <time.h>
#include <iomanip> // std::scientific
#include <utility> // std::pair
#include <vector> // std::vector

#include <model/Utilities/SequentialOutputFile.h>

//#include <model/SpaceDecomposition/SpatialCellProperty.h>

#include <model/ParticleInteraction/ParticleSystemBase.h>
#include <model/ParticleInteraction/SystemProperties.h>
#include <model/ParticleInteraction/PointSource.h>
#include <model/ParticleInteraction/FieldPoint.h>



namespace model {
    
    
    
    
    /**************************************************************************/
    /**************************************************************************/
    /*! \brief Serial Version
     */
    template <typename _ParticleType, typename UserSystemProperties = SystemProperties<> >
    class ParticleSystemSerial :
    /* inheritance */  public ParticleSystemBase<_ParticleType,UserSystemProperties>
    {
        typedef ParticleSystemBase<_ParticleType,UserSystemProperties> ParticleSystemBaseType;
        typedef _ParticleType ParticleType; // make ParticleType available outside the class
        typedef typename ParticleSystemBaseType::SpatialCellType SpatialCellType;
        typedef typename ParticleSystemBaseType::ParticleContainerType ParticleContainerType;
        
        
        
    public:

        
        /**********************************************************************/
        template <typename FieldType>
        void computeNeighborField()
        {
#ifdef _OPENMP
#pragma omp parallel for
#endif
            //! -1 loop over all particles in the ParticleSystem (parallelized in OpenMP)
            for (unsigned int k=0; k<this->size();++k)
            {
                //! -2 loop over neighbor cells of current particle
                for (typename SpatialCellType::CellMapType::const_iterator cIter =this->operator[](k).neighborCellsBegin();
                     /*                                                 */ cIter!=this->operator[](k).neighborCellsEnd();
                     /*                                               */ ++cIter)
                {
                    //! -3 loop over particles in the current neighbor cell
                    for(typename SpatialCellType::ParticleContainerType::const_iterator qIter =cIter->second->particleBegin();
                        /*                                                           */ qIter!=cIter->second->particleEnd();
                        /*                                                         */ ++qIter)
                    {
//                        FieldType(this->operator[](k),**qIter);  // the constructor of FieldType actually computes the binary interaction between *iterI and *iterJ
//                        FieldType::compute(*pIter->second,**qIter);  // the constructor of FieldType actually computes the binary interaction between *iterI and *iterJ
//                        *static_cast<*FieldPointParticle<FieldType> >(this->operator[](k)) += (*qIter)->computeAt(this->operator[](k));
                        
//                        *static_cast<*FieldPointParticle<FieldType> >(this->operator[](k)) += static_cast<SourcePoint<FieldType<ParticleType>>*>(*qIter)->computeAt(this->operator[](k));
//                        *static_cast<*FieldPointParticle<FieldType> >(this->operator[](k)) +=
                        *static_cast<FieldPointBase<ParticleType,FieldType>* const>(&this->operator[](k)) += FieldType::compute(**qIter,this->operator[](k));

//                        FieldType::operate
                        
//                        static_cast<FieldPoint<ParticleType,3,FieldType>* const>(&this->operator[](k));
                        
//                        static_cast<const PointSource<ParticleType,3,FieldType>* const>(*qIter);

                    }
                }
            }
            //                pIter->second->computeStress(); 
//            std::cout<<std::scientific<<" ["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]."<<std::endl;
        }
        
        
//        /**********************************************************************/
//        template <typename FieldType>
//        void computeAtFieldPoints(std::deque<FieldPoint<FieldType>*>& fieldDeq)
//        {
////            // loop over all particles
//#ifdef _OPENMP
//#pragma omp parallel for
//#endif
//            for (unsigned int k=0; k<fieldDeq.size();++k)
//            {// loop over neighbor cells
//            
//            
//                for (typename SpatialCellType::CellMapType::const_iterator cIter =this->operator[](k).neighborCellsBegin();
//                     /*                                                 */ cIter!=this->operator[](k).neighborCellsEnd();
//                     /*                                               */ ++cIter)
//                {// loop over particles in the neighbor cell
//                    for(typename SpatialCellType::ParticleContainerType::const_iterator qIter =cIter->second->particleBegin();
//                        /*                                                           */ qIter!=cIter->second->particleEnd();
//                        /*                                                         */ ++qIter)
//                    {
//                        //                        FieldType(this->operator[](k),**qIter);  // the constructor of FieldType actually computes the binary interaction between *iterI and *iterJ
//                        //                        FieldType::compute(*pIter->second,**qIter);  // the constructor of FieldType actually computes the binary interaction between *iterI and *iterJ
//                        *static_cast<*FieldPointParticle<FieldType> >(this->operator[](k)) += (*qIter)->computeAt(this->operator[](k));
//                    }
//                }
//            }
//            //                pIter->second->computeStress();
//            //            std::cout<<std::scientific<<" ["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]."<<std::endl;
//        }
        
 
        
        /**********************************************************************/
        template <class T>
        friend T& operator<< (T& os, const ParticleContainerType& pS)
        {/*! Operator<< use ParticleType-specific operator<<
          */
            for (typename ParticleContainerType::const_iterator pIter=pS.begin();pIter!=pS.end();++pIter)
            {
                os<<(*pIter)<<std::endl;
            }
            return os;
        }
        
    };
    
    
    
    
    
} // end namespace
#endif
                
                
                
                //        /**********************************************************************/
                //        template <typename CellProperty>
                //        void computeCellProperty() const
                //        {
                ////#ifdef _OPENMP
                //
                //
                ////#else
                //            for (typename SpatialCellType::CellMapType::const_iterator cIter=SpatialCellObserverType::cellBegin();cIter!=SpatialCellObserverType::cellEnd();++cIter)
                //            {
                //                CellProperty(*cIter->second);
                //            }
                ////#endif
                //        }

