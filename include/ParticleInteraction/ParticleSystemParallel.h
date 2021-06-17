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

#ifndef _MODEL_MPI_
#error YOU NEED TO DEFINE _MODEL_MPI_ TO USE THIS HEADER
#endif

#ifndef _MODEL_ParticleSystemParallel_h_
#define _MODEL_ParticleSystemParallel_h_

#ifdef _MODEL_MPI_
#include <ModelMPIbase.h>
//#include <metis.h> // partitioner
#endif

#ifdef _OPENMP
#include <omp.h>
#define EIGEN_DONT_PARALLELIZE // disable Eigen Internal openmp Parallelization
#endif


#include <vector>
#include <time.h>

#include <modelMacros.h> // model_checkInput
#include <SequentialOutputFile.h>
#include <CompareVectorsByComponent.h>
#include <PointSource.h>
#include <FieldPoint.h>
#include <ParticleSystemBase.h>
#include <LPTpartitioner.h>


namespace model {
    
    /**************************************************************************/
    /**************************************************************************/
    /*! \brief Parallel Version (openMPI + openMP)
     */
    template <typename _ParticleType>
    class ParticleSystemParallel :
    /* inheritance */  public ParticleSystemBase<_ParticleType>,
    /* inheritance */  public ModelMPIbase
    {
        
        typedef ParticleSystemBase<_ParticleType> ParticleSystemBaseType;
        typedef _ParticleType ParticleType; // make ParticleType available outside the class
        typedef typename ParticleSystemBaseType::SpatialCellObserverType SpatialCellObserverType;
        
    public:
        
        typedef typename ParticleSystemBaseType::SpatialCellType SpatialCellType;
        typedef typename ParticleSystemBaseType::ParticleContainerType ParticleContainerType;
        
        
        /**********************************************************************/
        ParticleSystemParallel(int argc, char* argv[]) :
        /* init list       */  ModelMPIbase(argc,argv)
        {
            
        }
        
        /**********************************************************************/
        ParticleSystemParallel()
        {
            
        }
        
        /**********************************************************************/
        template <typename ...AdditionalConstructorTypes>
        ParticleType* addParticle(const typename ParticleSystemBaseType::PositionType& p, const AdditionalConstructorTypes&... args)
        {/*! Adds a ParticleType particle with position p and an arbitrary (variadic)
          * number of additional constructor parameters.
          */
            assert(mpiInitialized() && "CANNOT INSERT PARTICLES BEFORE CALLING ModelMPIbase::initMPI(argc,argv).");
            //           partitionIsValid=false;
            return ParticleSystemBaseType::addParticle(p,args...);
        }
        
        
        //        /**********************************************************************/
        //        template <typename T>
        //        void partionSystem(LPTpartitioner<T>& partitioner)
        //        {/*! @param[in] useCellPartitioner Flag that determines whether the
        //          *  ParticleSystem is partitioned by particles or by cells.
        //          *
        //          *  Partitions particles among processors.
        //          */
        //
        //
        //            // Clear and call the LPT partitioner
        //            partitioner.clear();
        //
        //            for (typename ParticleContainerType::iterator pIter =this->begin();
        //                 /*                                    */ pIter!=this->end();
        //                 /*                                  */ ++pIter)
        //            {
        //                partitioner.insert(pIter->pCell->n2Weight(),&*pIter);
        //            }
        //
        //            partitioner.partition(this->mpiProcs());
        //        }
        
        /**********************************************************************/
        template <typename FieldType,typename... OtherSourceTypes>
        void computeNeighborField(const OtherSourceTypes&... otherSources)
        {/*! Compute nearest-neighbor particle interaction according to FieldType
          *\returns the number of interactions computed
          *
          * This member function performs the following operations:
          */
            
            MPI_Barrier(MPI_COMM_WORLD);
            
            
            //! -1 Partitioning:
            //! -1.1 creates a LPTpartitioner
            typedef LPTpartitioner<ParticleType> PartitionerType;
            PartitionerType lpt;
            
            //! -1.2 populate the LPTpartitioner using the particles stored in *this
            for (typename ParticleContainerType::iterator pIter =this->begin();
                 /*                                    */ pIter!=this->end();
                 /*                                  */ ++pIter)
            {
                //lpt.insert(pIter->pCell->n2Weight(),&*pIter);
                lpt.insert(pIter->pCell->neighborSize(),&*pIter);
            }
            //! -1.3 partition particles among MPI processes
            lpt.partition(this->mpiProcs());
            
            // Resize FieldPointBase<ParticleType,FieldType>::resultVector
            FieldPointBase<ParticleType,FieldType>::resize(this->size(),0.0);
            
            // Assign mpiID
            size_t binOffset(0);
            for (size_t b=0;b<lpt.bins();++b)
            {
                const size_t binSize(lpt.bin(b).size());
                for (size_t p=0;p<binSize;++p)
                {
                    lpt.bin(b)[p]->set_mpiID(binOffset+p);
                }
                binOffset+=binSize;
            }
            
            
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (unsigned int k=0; k<lpt.bin(this->mpiRank()).size();++k)
            {
                // Nearest-neighbor interaction
                static_cast<FieldPointBase<ParticleType,FieldType>* const>(lpt.bin(this->mpiRank())[k]) -> setZero();
                
                // loop over neighbor cells of the current particle
                for (typename SpatialCellType::CellMapType::const_iterator cIter=lpt.bin(this->mpiRank())[k]->neighborCellsBegin();cIter!=lpt.bin(this->mpiRank())[k]->neighborCellsEnd();++cIter)
                {
                    // loop over source particles in the neighbor cell
                    for(typename SpatialCellType::ParticleContainerType::const_iterator sIter=cIter->second->particleBegin();sIter!=cIter->second->particleEnd();++sIter) // loop over neighbor particles
                    {
                        //                        lpt.bin(this->mpiRank())[k]->template field<FieldType>() += FieldType::compute(**sIter,*lpt.bin(this->mpiRank())[k]);
                        if(static_cast<const SingleSourcePoint<ParticleType,FieldType>* const>(*sIter)->enabled)
                        {
                            *static_cast<FieldPointBase<ParticleType,FieldType>* const>(lpt.bin(this->mpiRank())[k]) += FieldType::compute(**sIter,*lpt.bin(this->mpiRank())[k]);
                        }
                    }
                }
                
                // Non nearest-neighbor interaction
                if(FieldType::use_multipole)
                {
                    typename SpatialCellType::CellMapType farCells(lpt.bin(this->mpiRank())[k]->template farCells<_ParticleType>());
                    //                    lpt.bin(this->mpiRank())[k]->template field<FieldType>() += FieldType::multipole(*lpt.bin(this->mpiRank())[k],farCells);
                    *static_cast<FieldPointBase<ParticleType,FieldType>* const>(lpt.bin(this->mpiRank())[k]) += FieldType::multipole(*lpt.bin(this->mpiRank())[k],farCells);
                }
                
                // Add contribution of other sources
                if(sizeof...(OtherSourceTypes))
                {
                    *static_cast<FieldPointBase<ParticleType,FieldType>* const>(lpt.bin(this->mpiRank())[k]) += FieldType::addSourceContribution(*lpt.bin(this->mpiRank())[k],otherSources...);
                }
            }
            
            //! 4- syncronize FieldType::resultVector among processors
            std::vector<int> interactionSizeVector(this->mpiProcs());
            std::vector<int> interactionRankOffsetVector(this->mpiProcs());
            
            size_t dofOffset(0);
            for (int k=0;k<this->mpiProcs();++k)
            {
                // Store the number
                interactionSizeVector[k]=lpt.bin(k).size()*FieldPointBase<ParticleType,FieldType>::DataPerParticle;
                interactionRankOffsetVector[k]=dofOffset;
                dofOffset+=lpt.bin(k).size()*FieldPointBase<ParticleType,FieldType>::DataPerParticle;
            }
            
            MPI_Allgatherv(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,
                           &FieldPointBase<ParticleType,FieldType>::resultVector[0],
                           &interactionSizeVector[0],&interactionRankOffsetVector[0],
                           MPI_DOUBLE,MPI_COMM_WORLD);
        }
        
        
        /**********************************************************************/
        template <typename OtherParticleType, typename FieldType,typename... OtherSourceTypes>
        void computeField(std::deque<OtherParticleType>& fpDeq,const OtherSourceTypes&... otherSources) const
        {/*! Compute nearest-neighbor particle interaction according to FieldType
          *\returns the number of interactions computed
          */
            
            MPI_Barrier(MPI_COMM_WORLD);
            
            //            LPTpartitioner<ParticleType> lpt;
            typedef LPTpartitioner<OtherParticleType> PartitionerType;
            PartitionerType lpt;
            //            partionSystem(lpt);
            
            for (typename std::deque<OtherParticleType>::iterator pIter =fpDeq.begin();
                 /*                                    */ pIter!=fpDeq.end();
                 /*                                  */ ++pIter)
            {
                typename SpatialCellType::CellMapType cellMap(pIter->template neighborCells<_ParticleType>());
                double weight(0.0);
                for (typename SpatialCellType::CellMapType::const_iterator cIter=cellMap.begin();cIter!=cellMap.end();++cIter)
                {
                    weight += cIter->second->size();
                }
                
                lpt.insert(weight,&*pIter);
            }
            
            lpt.partition(this->mpiProcs());
            
            
            
            // Resize FieldPointBase<OtherParticleType,FieldType>::resultVector
            FieldPointBase<OtherParticleType,FieldType>::resize(fpDeq.size(),0.0);
            
            // Assign mpiID
            size_t binOffset(0);
            for (size_t b=0;b<lpt.bins();++b)
            {
                const size_t binSize(lpt.bin(b).size());
                for (size_t p=0;p<binSize;++p)
                {
                    lpt.bin(b)[p]->set_mpiID(binOffset+p);
                }
                binOffset+=binSize;
            }
            
            
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (unsigned int k=0; k<lpt.bin(this->mpiRank()).size();++k)
            {
                // Nearest-neighbor interaction
                static_cast<FieldPointBase<OtherParticleType,FieldType>* const>(lpt.bin(this->mpiRank())[k]) -> setZero();
                
                typename SpatialCellType::CellMapType cellMap((lpt.bin(this->mpiRank())[k])->template neighborCells<_ParticleType>());
                
                // loop over neighbor cells of the current particle
                for (typename SpatialCellType::CellMapType::const_iterator cIter =cellMap.begin();
                     /*                                                 */ cIter!=cellMap.end();
                     /*                                               */ ++cIter)
                {
                    // loop over source particles in the neighbor cell
                    for(typename SpatialCellType::ParticleContainerType::const_iterator sIter=cIter->second->particleBegin();sIter!=cIter->second->particleEnd();++sIter) // loop over neighbor particles
                    {
                        //                        lpt.bin(this->mpiRank())[k]->template field<FieldType>() += FieldType::compute(**sIter,*lpt.bin(this->mpiRank())[k]);
                        if(static_cast<const SingleSourcePoint<ParticleType,FieldType>* const>(*sIter)->enabled)
                        {
                            *static_cast<FieldPointBase<OtherParticleType,FieldType>* const>(lpt.bin(this->mpiRank())[k]) += FieldType::compute(**sIter,*lpt.bin(this->mpiRank())[k]);
                        }
                    }
                }
                
                // Non nearest-neighbor interaction
                if(FieldType::use_multipole)
                {
                    typename SpatialCellType::CellMapType farCells(lpt.bin(this->mpiRank())[k]->template farCells<_ParticleType>());
                    //                    lpt.bin(this->mpiRank())[k]->template field<FieldType>() += FieldType::multipole(*lpt.bin(this->mpiRank())[k],farCells);
                    *static_cast<FieldPointBase<OtherParticleType,FieldType>* const>(lpt.bin(this->mpiRank())[k]) += FieldType::multipole(*lpt.bin(this->mpiRank())[k],farCells);
                    
                }
                
                // Add contribution of other sources
                if(sizeof...(OtherSourceTypes))
                {
                    *static_cast<FieldPointBase<OtherParticleType,FieldType>* const>(lpt.bin(this->mpiRank())[k]) += FieldType::addSourceContribution(*lpt.bin(this->mpiRank())[k],otherSources...);
                }
                
            }
            
            
            std::vector<int> interactionSizeVector(this->mpiProcs());
            std::vector<int> interactionRankOffsetVector(this->mpiProcs());
            
            size_t dofOffset(0);
            for (int k=0;k<this->mpiProcs();++k)
            {
                // Store the number
                interactionSizeVector[k]=lpt.bin(k).size()*FieldPointBase<OtherParticleType,FieldType>::DataPerParticle;
                interactionRankOffsetVector[k]=dofOffset;
                dofOffset+=lpt.bin(k).size()*FieldPointBase<OtherParticleType,FieldType>::DataPerParticle;
            }
            
            //! Syncronize FieldType::resultVector among processors
            MPI_Allgatherv(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,
                           &FieldPointBase<OtherParticleType,FieldType>::resultVector[0],
                           &interactionSizeVector[0],&interactionRankOffsetVector[0],
                           MPI_DOUBLE,MPI_COMM_WORLD);
        }
        
        /**********************************************************************/
        template <typename OtherParticleType, typename FieldType,typename... OtherSourceTypes>
        void computeField(OtherParticleType& part, const OtherSourceTypes&... otherSources) const
        {/*!@param[in] part
          * @param[in] otherSources
          * \brief computes the field FieldType for the only particle part,
          * without parallelization
          */
            ParticleSystemBaseType::template computeField<OtherParticleType,FieldType,OtherSourceTypes...>(part,otherSources...);
        }
        
        /**********************************************************************/
        template <class T>
        friend T& operator<< (T& os, const ParticleContainerType& pS)
        {/*! Operator << use ParticleType specific << operator
          */
            if (ModelMPIbase::mpiRank()==0)
            {
                for (typename ParticleContainerType::const_iterator pIter=pS.begin();pIter!=pS.end();++pIter)
                {
                    os<<(*pIter);
                }
            }
            return os;
        }
        
    };
    
    
} // end namespace
#endif


//        /**********************************************************************/
//        template <typename FieldType>
//        typename FieldType::ResultType getInteractionResult(const int& pID)
//        {
//            return FieldType::getResult(*(this->find(pID)->second));
//        }
//
//        /**********************************************************************/
//        template <typename FieldType>
//        void resetInteraction()
//        {/*! Compute full particle interaction according to FieldType
//          */
//            for (typename ParticleContainerType::iterator iterJ=this->begin();iterJ!=this->end();++iterJ) // FOR SYMMETRIC INTERACTION iterJ STARTS AT iterI
//            {
//                FieldType::reset(*iterJ->second);  // the constructor of FieldType actually computes the binary interaction between *iterI and *iterJ
//            }
//        }



/**********************************************************************/
//        void partionByCells()
//        {
//            typedef typename SpatialCellObserverType::CellIdType CellIdType;
//            typedef std::map<CellIdType,const int,CompareVectorsByComponent<int,_ParticleType::dim> > ijk2idMapType;
//            //            typedef std::map<Eigen::Matrix<int,_ParticleType::dim,1>,const int,CompareVectorsByComponent<int,_ParticleType::dim> > ijk2idMapType;
//            ijk2idMapType  ijk2idMap;;
//            typedef std::map<int,const CellIdType> id2ijkMapType;
//            //            typedef std::map<int,Eigen::Matrix<int,_ParticleType::dim,1> > id2ijkMapType;
//            id2ijkMapType id2ijkMap;;
//
//            // Populate ijk2idMap and id2ijkMap
//            int k(0);
//            for (typename SpatialCellObserverType::CellMapType::const_iterator cIter =SpatialCellObserverType::cellBegin();
//                 /*                                                         */ cIter!=SpatialCellObserverType::cellEnd();
//                 /*                                                         */ ++cIter)
//            {
//                ijk2idMap.insert(std::pair<Eigen::Matrix<int,_ParticleType::dim,1>,const int>(cIter->second->cellID,k));
//                id2ijkMap.insert(std::pair<int,Eigen::Matrix<int,_ParticleType::dim,1> >(k,cIter->second->cellID));
//                ++k;
//            }
//
//
//
//            // Compute weightScaleFactor to avoid int overflow in std::vector<int> vWeights
//            double totalWeightD(0.0);
//            for (typename SpatialCellObserverType::CellMapType::const_iterator cIter =SpatialCellObserverType::cellBegin();
//                 /*                                                         */ cIter!=SpatialCellObserverType::cellEnd();
//                 /*                                                       */ ++cIter)
//            {
//                totalWeightD+=cIter->second->n2WeightD();
//            }
//            double weightScaleFactor(1.0);
//            if (totalWeightD>=INT_MAX) // need to rescale actual "int" weights
//            {
//                weightScaleFactor=totalWeightD/INT_MAX;
//                if (this->mpiRank() == 0)
//                {
//
//                    std::cout<<"detected weight overflow. Rescaling weights: weightScaleFactor="<<weightScaleFactor<<std::endl;
//                }
//            }
//
//            // Assemble std::vector<int> vWeights
//            std::vector<int> vWeights(SpatialCellObserverType::totalCells(),0); // vector containing the (scaled) weight of each cell
//            int kk(0);
//            int totalWeight(0);
//            for (typename SpatialCellObserverType::CellMapType::const_iterator cIter =SpatialCellObserverType::cellBegin();
//                 /*                                                         */ cIter!=SpatialCellObserverType::cellEnd();
//                 /*                                                       */ ++cIter)
//            {
//                vWeights[kk]=cIter->second->n2Weight()/weightScaleFactor;
//                totalWeight+=vWeights[kk];
//                ++kk;
//            }
//
//            // Check overflow
//            if (totalWeightD>=INT_MAX) // actual totalWeight will overflow: need to rescale actual "int" weights
//            {
//                const double actualScaleFactor(totalWeightD/totalWeight);
//                if (this->mpiRank() == 0)
//                {
//                    std::cout<<"totalWeight="<<totalWeight<<", totalWeightD="<<totalWeightD<<", actualScaleFactor="<<actualScaleFactor<<std::endl;
//                }
//                model_removeAssert(std::fabs(actualScaleFactor-weightScaleFactor)/weightScaleFactor<0.01 && "WEIGHT SCALE FACTOR MISMATCH");
//            }
//
//
//            // - Populate vector of neighbors and offsets
//            std::vector<int> vNeighbors;
//            std::vector<int> vOffsets;
//
//            int offSet(0);
//            for (typename SpatialCellObserverType::CellMapType::const_iterator cIter =SpatialCellObserverType::cellBegin();
//                 /*                                                         */ cIter!=SpatialCellObserverType::cellEnd();
//                 /*                                                       */ ++cIter)
//            {
//                vOffsets.push_back(offSet);
//                for (typename SpatialCellObserverType::CellMapType::const_iterator nIter =cIter->second->neighborCellsBegin();
//                     /*                                                         */ nIter!=cIter->second->neighborCellsEnd();
//                     /*                                                       */ ++nIter)
//                {
//                    if (cIter->second->cellID!=nIter->second->cellID)
//                    {
//                        typename ijk2idMapType::const_iterator mIter(ijk2idMap.find(nIter->second->cellID.transpose()));
//                        vNeighbors.push_back(mIter->second);
//                        ++offSet;
//                    }
//                }
//            }
//            vOffsets.push_back(offSet); // push back an extra offset (required by Metis)
//
//
//            // Create a vector to store the output of the partitioner. Initialize at process=0
//            std::vector<int> partitionerRankVector(SpatialCellObserverType::totalCells(),0);
//
//            // partition into the same number of MPI processes
//            int numPart = this->mpiProcs();
//
//            if (numPart > 1)
//            {
//                int metisOptions[METIS_NOPTIONS];
//                int numVert = SpatialCellObserverType::totalCells(); // number of nodes
//                int ncon = 1;  // number of weights per node
//                int objval;
//                float ubvec = 1.001; // load imbalance tolerance
//                METIS_SetDefaultOptions(metisOptions);
//                METIS_PartGraphKway(&numVert,&ncon,&vOffsets[0],&vNeighbors[0],
//                                    &vWeights[0],NULL,NULL,&numPart,NULL,
//                                    &ubvec,metisOptions,&objval,&partitionerRankVector[0]);
//            }
//
//
//
//            // Assign each cell the corresponding processor, as found by metis
//            //            int assignedParticleSize(0);
//            for (unsigned int i=0; i<SpatialCellObserverType::totalCells(); ++i)
//            {
//                typename id2ijkMapType::const_iterator idIter(id2ijkMap.find(i));
//                const CellIdType cellID(idIter->second);
//
//                typename SpatialCellObserverType::isCellType isC(SpatialCellObserverType::isCell(cellID));
//                model_removeAssert(isC.first && "CELL MUST EXIST!");
//                isC.second->assignedRank=partitionerRankVector[i];
//
//                if (this->mpiRank()==partitionerRankVector[i])
//                {
//                    //                    assignedCells.push_back(isC.second);
//                    // loop over particles in the current cell
//                    for (typename SpatialCellType::ParticleContainerType::const_iterator pIter =isC.second->particleBegin();
//                         /*                                                           */ pIter!=isC.second->particleEnd();
//                         /*                                                         */ ++pIter)
//                    {
//
//                        const typename ParticleContainerType::iterator nIter(this->find((*pIter)->sID));
//                        //                        assert(nIter!=this->end());
//                        assignedParticles.insert(std::make_pair((*pIter)->sID,nIter->second));
//                        //                        assignedParticles.push_back(nIter->second);
//                    }
//
//
//                    //                    assignedParticleSize+=isC.second->size();
//                }
//            }
//
//            // evaluate Performance
//            evaluatePerformance(partitionerRankVector,vWeights,totalWeight);
//
//            //            return assignedParticles.size();
//            //            return assignedParticleSize;
//        }

/**********************************************************************/
//        void partionByParticles()
//        {
//            // Initialize the vector containing the weight of each particle
//            std::vector<int> vWeights(this->size(),1);
//
//            // Populate the vector containing the weight of each particle
//            int kk(0);
//            for (typename ParticleContainerType::const_iterator pIter =this->begin();
//                 /*                                          */ pIter!=this->end();
//                 /*                                        */ ++pIter)
//            {
//                // The weight of the kk-th particle is the n2Weight() of its cell
//                vWeights[kk]=pIter->second->pCell->n2Weight();
//                ++kk;
//            }
//
//            // Initialize METIS vectors
//            std::vector<int> vNeighbors;
//            std::vector<int> vOffsets;
//
//            // Populate METIS vectors
//            int offSet(0);
//            int pp(0);
//            for (typename ParticleContainerType::const_iterator pIter =this->begin();
//                 /*                                          */ pIter!=this->end();
//                 /*                                        */ ++pIter)
//            {
//                vOffsets.push_back(offSet);
//                vNeighbors.push_back(pp); // each particle is its only neighbor
//                ++offSet;
//                ++pp;
//            }
//            vOffsets.push_back(offSet); // push back an extra offset (required by Metis)
//
//
//            // partition into the same number of MPI processes
//            int numPart = this->mpiProcs();
//
//            // Create a vector to store the output of the partitioner. Initialize at process=0
//            std::vector<int> partitionerRankVector(this->size(),0); //
//
//            // Call METIS
//            if (numPart > 1)
//            {
//                int metisOptions[METIS_NOPTIONS];
//                int numVert = this->size(); // number of nodes
//                int ncon = 1;  // number of weights per node
//                int objval;
//                float ubvec = 1.001; // load imbalance tolerance
//                METIS_SetDefaultOptions(metisOptions);
//                METIS_PartGraphKway(&numVert,&ncon,&vOffsets[0],&vNeighbors[0],
//                                    &vWeights[0],NULL,NULL,&numPart,NULL,
//                                    &ubvec,metisOptions,&objval,&partitionerRankVector[0]);
//            }
//
//            //            int assignedParticleSize(0);
//
//            for (unsigned int i=0; i<this->size(); ++i)
//            {
//                if (this->mpiRank()==partitionerRankVector[i])
//                {
//                    typename ParticleContainerType::iterator pIter(this->begin());
//                    std::advance(pIter,i);
//                    //                    assignedParticles.push_back(pIter->second);
//                    assignedParticles.insert(std::make_pair(pIter->second->sID,pIter->second));
//
//                    //                    assignedParticleSize+=1;
//                }
//            }
//
//            //          return assignedParticleSize;
//        }

