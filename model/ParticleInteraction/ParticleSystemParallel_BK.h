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
#include <model/MPI/ModelMPIbase.h>
#include <metis.h> // partitioner
#endif

#include <utility> // for std::pair
#include <map>
#include <vector>
#include <deque>
#include <time.h>
#include <iomanip> // std::scientific

#include <Eigen/Core>

#include <model/Utilities/modelMacros.h> // model_checkInput
#include <model/Utilities/SequentialOutputFile.h>
#include <model/Utilities/CompareVectorsByComponent.h>
#include <model/ParticleInteraction/ParticleSystemBase.h>
#include <model/ParticleInteraction/SystemProperties.h>
#include <model/MPI/LongestProcessingTime.h>

//#include <model/SpaceDecomposition/SpatialCellObserver.h>


namespace model {
    
    /**************************************************************************/
    /**************************************************************************/
    /*! \brief Parallel Version (MPI)
     */
    template <typename _ParticleType, typename UserSystemProperties>
    class ParticleSystemParallel :
    /* inheritance */  public ParticleSystemBase<_ParticleType,UserSystemProperties>,
    /* inheritance */  public ModelMPIbase
    {
        //enum {dim=_ParticleType::dim};
        //        typedef SpatialCell<_ParticleType,_ParticleType::dim> SpatialCellType;
        
        typedef ParticleSystemBase<_ParticleType,UserSystemProperties> ParticleSystemBaseType;
        typedef _ParticleType ParticleType; // make ParticleType available outside the class
        typedef typename ParticleSystemBaseType::SpatialCellType SpatialCellType;
        typedef typename ParticleSystemBaseType::ParticleContainerType ParticleContainerType;
        typedef typename ParticleSystemBaseType::SpatialCellObserverType SpatialCellObserverType;
        //        typedef std::deque<SpatialCellType*> AssignedCellContainerType;
        typedef std::map<size_t,_ParticleType*>   AssignedParticleContainerType;
        
        //! The SpatialCell cells assigned to this MPI process
        //        AssignedCellContainerType assignedCells;
        AssignedParticleContainerType assignedParticles;
        
        
        //        bool useCellPartitioner;
        
        bool partitionIsValid;
        std::vector<int> mpiIDoffsetVector;
        std::vector<int> particleSizeVector;
        
        /**********************************************************************/
        void evaluatePerformance(const std::vector<int>& partitionerRankVector,
                                 const std::vector<int>& vWeights,
                                 const int& totalWeight) const
        {
            std::vector<int> procWeights(this->mpiProcs(),0); // vector containing the weight of each processor
            
            for (unsigned int w=0;w<partitionerRankVector.size();++w)
            {
                procWeights[partitionerRankVector[w]]+=vWeights[w];
            }
            
            if (mpiRank()==0)
            {
                std::cout<<"Processor weights: ";
                for (unsigned int w=0;w<procWeights.size();++w)
                {
                    std::cout<<static_cast<float>(procWeights[w])/totalWeight*100.0f<<"%, ";
//                    double speedup(this->mpiProcs());
//                    float temp(static_cast<float>(totalWeight)/procWeights[w]);
//                    if(speedup>temp)
//                    {
//                        speedup=temp;
//                    }
                }
                std::cout<<std::endl;

            }
            
            //            std::cout<<"Expected speed-up="<<speedup<<std::endl;
        }
        
        
        /**********************************************************************/
        void partionByCells()
        {
            typedef typename SpatialCellObserverType::CellIdType CellIdType;
            typedef std::map<CellIdType,const int,CompareVectorsByComponent<int,_ParticleType::dim> > ijk2idMapType;
            //            typedef std::map<Eigen::Matrix<int,_ParticleType::dim,1>,const int,CompareVectorsByComponent<int,_ParticleType::dim> > ijk2idMapType;
            ijk2idMapType  ijk2idMap;;
            typedef std::map<int,const CellIdType> id2ijkMapType;
            //            typedef std::map<int,Eigen::Matrix<int,_ParticleType::dim,1> > id2ijkMapType;
            id2ijkMapType id2ijkMap;;
            
            // Populate ijk2idMap and id2ijkMap
            int k(0);
            for (typename SpatialCellObserverType::CellMapType::const_iterator cIter =SpatialCellObserverType::cellBegin();
                 /*                                                         */ cIter!=SpatialCellObserverType::cellEnd();
                 /*                                                         */ ++cIter)
            {
                ijk2idMap.insert(std::pair<Eigen::Matrix<int,_ParticleType::dim,1>,const int>(cIter->second->cellID,k));
                id2ijkMap.insert(std::pair<int,Eigen::Matrix<int,_ParticleType::dim,1> >(k,cIter->second->cellID));
                ++k;
            }
            
            
            
            // Compute weightScaleFactor to avoid int overflow in std::vector<int> vWeights
            double totalWeightD(0.0);
            for (typename SpatialCellObserverType::CellMapType::const_iterator cIter =SpatialCellObserverType::cellBegin();
                 /*                                                         */ cIter!=SpatialCellObserverType::cellEnd();
                 /*                                                       */ ++cIter)
            {
                totalWeightD+=cIter->second->n2WeightD();
            }
            double weightScaleFactor(1.0);
            if (totalWeightD>=INT_MAX) // need to rescale actual "int" weights
            {
                weightScaleFactor=totalWeightD/INT_MAX;
                if (this->mpiRank() == 0)
                {
                    
                    std::cout<<"detected weight overflow. Rescaling weights: weightScaleFactor="<<weightScaleFactor<<std::endl;
                }
            }
            
            // Assemble std::vector<int> vWeights
            std::vector<int> vWeights(SpatialCellObserverType::totalCells(),0); // vector containing the (scaled) weight of each cell
            int kk(0);
            int totalWeight(0);
            for (typename SpatialCellObserverType::CellMapType::const_iterator cIter =SpatialCellObserverType::cellBegin();
                 /*                                                         */ cIter!=SpatialCellObserverType::cellEnd();
                 /*                                                       */ ++cIter)
            {
                vWeights[kk]=cIter->second->n2Weight()/weightScaleFactor;
                totalWeight+=vWeights[kk];
                ++kk;
            }
            
            // Check overflow
            if (totalWeightD>=INT_MAX) // actual totalWeight will overflow: need to rescale actual "int" weights
            {
                const double actualScaleFactor(totalWeightD/totalWeight);
                if (this->mpiRank() == 0)
                {
                    std::cout<<"totalWeight="<<totalWeight<<", totalWeightD="<<totalWeightD<<", actualScaleFactor="<<actualScaleFactor<<std::endl;
                }
                model_removeAssert(std::fabs(actualScaleFactor-weightScaleFactor)/weightScaleFactor<0.01 && "WEIGHT SCALE FACTOR MISMATCH");
            }
            
            
            // - Populate vector of neighbors and offsets
            std::vector<int> vNeighbors;
            std::vector<int> vOffsets;
            
            int offSet(0);
            for (typename SpatialCellObserverType::CellMapType::const_iterator cIter =SpatialCellObserverType::cellBegin();
                 /*                                                         */ cIter!=SpatialCellObserverType::cellEnd();
                 /*                                                       */ ++cIter)
            {
                vOffsets.push_back(offSet);
                for (typename SpatialCellObserverType::CellMapType::const_iterator nIter =cIter->second->neighborCellsBegin();
                     /*                                                         */ nIter!=cIter->second->neighborCellsEnd();
                     /*                                                       */ ++nIter)
                {
                    if (cIter->second->cellID!=nIter->second->cellID)
                    {
                        typename ijk2idMapType::const_iterator mIter(ijk2idMap.find(nIter->second->cellID.transpose()));
                        vNeighbors.push_back(mIter->second);
                        ++offSet;
                    }
                }
            }
            vOffsets.push_back(offSet); // push back an extra offset (required by Metis)
            
            
            // Create a vector to store the output of the partitioner. Initialize at process=0
            std::vector<int> partitionerRankVector(SpatialCellObserverType::totalCells(),0);
            
            // partition into the same number of MPI processes
            int numPart = this->mpiProcs();
            
            if (numPart > 1)
            {
                int metisOptions[METIS_NOPTIONS];
                int numVert = SpatialCellObserverType::totalCells(); // number of nodes
                int ncon = 1;  // number of weights per node
                int objval;
                float ubvec = 1.001; // load imbalance tolerance
                METIS_SetDefaultOptions(metisOptions);
                METIS_PartGraphKway(&numVert,&ncon,&vOffsets[0],&vNeighbors[0],
                                    &vWeights[0],NULL,NULL,&numPart,NULL,
                                    &ubvec,metisOptions,&objval,&partitionerRankVector[0]);
            }
            
            
            
            // Assign each cell the corresponding processor, as found by metis
            //            int assignedParticleSize(0);
            for (unsigned int i=0; i<SpatialCellObserverType::totalCells(); ++i)
            {
                typename id2ijkMapType::const_iterator idIter(id2ijkMap.find(i));
                const CellIdType cellID(idIter->second);
                
                typename SpatialCellObserverType::isCellType isC(SpatialCellObserverType::isCell(cellID));
                model_removeAssert(isC.first && "CELL MUST EXIST!");
                isC.second->assignedRank=partitionerRankVector[i];
                
                if (this->mpiRank()==partitionerRankVector[i])
                {
                    //                    assignedCells.push_back(isC.second);
                    // loop over particles in the current cell
                    for (typename SpatialCellType::ParticleContainerType::const_iterator pIter =isC.second->particleBegin();
                         /*                                                           */ pIter!=isC.second->particleEnd();
                         /*                                                         */ ++pIter)
                    {
                        
                        const typename ParticleContainerType::iterator nIter(this->find((*pIter)->sID));
                        //                        assert(nIter!=this->end());
                        assignedParticles.insert(std::make_pair((*pIter)->sID,nIter->second));
                        //                        assignedParticles.push_back(nIter->second);
                    }
                    
                    
                    //                    assignedParticleSize+=isC.second->size();
                }
            }
            
            // evaluate Performance
            evaluatePerformance(partitionerRankVector,vWeights,totalWeight);
            
            //            return assignedParticles.size();
            //            return assignedParticleSize;
        }
        
        /**********************************************************************/
        void partionByParticles()
        {
            // Initialize the vector containing the weight of each particle
            std::vector<int> vWeights(this->size(),1);
            
            // Populate the vector containing the weight of each particle
            int kk(0);
            for (typename ParticleContainerType::const_iterator pIter =this->begin();
                 /*                                          */ pIter!=this->end();
                 /*                                        */ ++pIter)
            {
                // The weight of the kk-th particle is the n2Weight() of its cell
                vWeights[kk]=pIter->second->pCell->n2Weight();
                ++kk;
            }
            
            // Initialize METIS vectors
            std::vector<int> vNeighbors;
            std::vector<int> vOffsets;
            
            // Populate METIS vectors
            int offSet(0);
            int pp(0);
            for (typename ParticleContainerType::const_iterator pIter =this->begin();
                 /*                                          */ pIter!=this->end();
                 /*                                        */ ++pIter)
            {
                vOffsets.push_back(offSet);
                vNeighbors.push_back(pp); // each particle is its only neighbor
                ++offSet;
                ++pp;
            }
            vOffsets.push_back(offSet); // push back an extra offset (required by Metis)
            
            
            // partition into the same number of MPI processes
            int numPart = this->mpiProcs();
            
            // Create a vector to store the output of the partitioner. Initialize at process=0
            std::vector<int> partitionerRankVector(this->size(),0); //
            
            // Call METIS
            if (numPart > 1)
            {
                int metisOptions[METIS_NOPTIONS];
                int numVert = this->size(); // number of nodes
                int ncon = 1;  // number of weights per node
                int objval;
                float ubvec = 1.001; // load imbalance tolerance
                METIS_SetDefaultOptions(metisOptions);
                METIS_PartGraphKway(&numVert,&ncon,&vOffsets[0],&vNeighbors[0],
                                    &vWeights[0],NULL,NULL,&numPart,NULL,
                                    &ubvec,metisOptions,&objval,&partitionerRankVector[0]);
            }
            
            //            int assignedParticleSize(0);
            
            for (unsigned int i=0; i<this->size(); ++i)
            {
                if (this->mpiRank()==partitionerRankVector[i])
                {
                    typename ParticleContainerType::iterator pIter(this->begin());
                    std::advance(pIter,i);
                    //                    assignedParticles.push_back(pIter->second);
                    assignedParticles.insert(std::make_pair(pIter->second->sID,pIter->second));
                    
                    //                    assignedParticleSize+=1;
                }
            }
            
            //          return assignedParticleSize;
        }
        
        
        
        
        
    public:
        
        
        
        
        /**********************************************************************/
        ParticleSystemParallel(int argc, char* argv[]) :
//        /* init list       */  ParticleSystemBaseType::ParticleSystemBase(cellSize),
        /* init list       */  ModelMPIbase(argc,argv),
        /* init list       */  partitionIsValid(false)
        {
//            particleSizeVector.resize(this->mpiProcs());
//            mpiIDoffsetVector.resize(this->mpiProcs());
        }
        
        /**********************************************************************/
        ParticleSystemParallel() :
        /* init list       */  partitionIsValid(false)
        {

        }
        
        /**********************************************************************/
        template <typename ...AdditionalConstructorTypes>
        ParticleType* addParticle(const typename ParticleSystemBaseType::PositionType& p, const AdditionalConstructorTypes&... args)
        {/*! Adds a ParticleType particle with position p and an arbitrary (variadic)
          * number of additional constructor parameters.
          */
            assert(mpiInitialized() && "CANNOT INSERT PARTICLES BEFORE CALLING ModelMPIbase::initMPI(argc,argv).");
            partitionIsValid=false;
            return ParticleSystemBaseType::addParticle(p,args...);
        }
        
        
        /**********************************************************************/
        float partionSystem()
        {/*! @param[in] useCellPartitioner Flag that determines whether the
          *  ParticleSystem is partitioned by particles or by cells.
          *
          *  Partitions particles among processors.
          */
            
            //          useCellPartitioner=useCellPartitioner_in;
            
            //            const float t0(clock());
            //            assignedCells.clear();
            
            //! -1 Clear the container of "assignedParticles"
            assignedParticles.clear();
            //            assignedParticles.reserve(this->size());
            
            //! -2 The number of ParticleType particles assigned to this MPI process
            //            int assignedParticleSize(0);
            
            if(this->useCellPartitioner)
            {
                partionByCells();
            }
            else
            {
                partionByParticles();
            }
            
            
            LongestProcessingTime<ParticleType> lpt;
            
            for (typename ParticleContainerType::iterator pIter =this->begin();
                 /*                                          */ pIter!=this->end();
                 /*                                        */ ++pIter)
            {
                lpt.insert(pIter->second->pCell->n2Weight(),pIter->second);
            }
            lpt.partition(this->mpiProcs());
            
            if (this->mpiRank()==0){lpt.show();}
            
 //           DON'T USE METIS !!!!
            
            /*! -3 Assign each particle in the "assignedParticles" container a unique ID (mpiID)
             *  All mpiID(s) of particles assigned to a given processor are sequential
             */
            
            //!     -3.1 Collect assignedParticleSize from each MPI process into particleSizeVector
            particleSizeVector.resize(this->mpiProcs());
            int assignedParticleSize(assignedParticles.size());
            MPI_Allgather(&assignedParticleSize,1,MPI_INT,&particleSizeVector[0],1,MPI_INT,MPI_COMM_WORLD);
            
            // Check that all particles have been partitioned
            unsigned int totalAssignedParticle(0);
            for (unsigned int k=0;k<particleSizeVector.size();++k)
            {
                totalAssignedParticle+=particleSizeVector[k];
            }
            model_removeAssert(totalAssignedParticle==this->size() && "SOMETHING WENT VERY WRONG");
            
            //!     -3.2 Compute mpiIDoffset by comulative summing assignedParticleSize across the MPI processes
            int mpiIDoffset(0);
            MPI_Exscan(&assignedParticleSize,&mpiIDoffset,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
            
            //!     -3.3 Collect mpiIDoffset from each MPI process into mpiIDoffsetVector
            mpiIDoffsetVector.resize(this->mpiProcs());
            MPI_Allgather(&mpiIDoffset,1,MPI_INT,&mpiIDoffsetVector[0],1,MPI_INT,MPI_COMM_WORLD);
            
            //!     -3.4 Loop over assignedParticles and assign an increasing ID with offset mpiIDoffsetVector[this->mpiRank()]
            std::vector<size_t> mpiIDvector(this->size(),0);
            size_t mpiLocalID(0);
            
            for (typename AssignedParticleContainerType::iterator pIter =assignedParticles.begin();
                 /*                                            */ pIter!=assignedParticles.end();
                 /*                                          */ ++pIter)
            {
                mpiIDvector[mpiLocalID+mpiIDoffsetVector[this->mpiRank()]]=pIter->second->sID;
                //                (*pIter)->mpiID=mpiLocalID+mpiIDoffsetVector[this->mpiRank()];
                ++mpiLocalID;
            }
            
            
            //            //! Syncronize InteractionType::resultVector among processors
            MPI_Allgatherv(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,
                           &mpiIDvector[0],&particleSizeVector[0],&mpiIDoffsetVector[0],MPI_UNSIGNED_LONG,MPI_COMM_WORLD);
            
            
            for (size_t k=0;k<this->size();++k)
            {
                const typename ParticleContainerType::iterator pIter(this->find(mpiIDvector[k]));
                pIter->second->mpiID=k;
            }
            
            //            size_t mpiLocalID(0);
            //            for (typename AssignedParticleContainerType::iterator pIter =assignedParticles.begin();
            //                 /*                                            */ pIter!=assignedParticles.end();
            //                 /*                                            */ pIter++)
            //            {
            //                (*pIter)->mpiID=mpiLocalID+mpiIDoffsetVector[this->mpiRank()];
            //                mpiLocalID++;
            //            }
            
            // Set flag "partitionIsValid" to true
            partitionIsValid=true;
            
            //            if (this->mpiRank() == 0)
            //            {
            //            std::cout<<std::scientific<<" ["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]."<<std::endl;
            //            }
            
            
            //            return HERE RETURN THE LOAD IMBALANCE;
            return 1.0;
            
        }
        
        /**********************************************************************/
        template <typename InteractionType>
        void computeNeighborInteraction()
        {/*! Compute nearest-neighbor particle interaction according to InteractionType
          */
            
            
            if (!partitionIsValid)
            {
                partionSystem();
            }
            model_removeAssert(partitionIsValid && "PARTITION IS NOT VALID");
            
            // Construct a vector to store the number of particles in each processor
            //            std::vector<int> particleSizeVector(this->mpiProcs(),0);
            //            particleSizeVector.resize(this->mpiProcs());
            
            
            //            assignMpiIDs(particleSizeVector);
            
            double t0(clock());
            
            
            //            InteractionType::resultVector.resize(3*this->size(),0.0);
            InteractionType::resize(3*this->size(),0.0);
            
            //           size_t Ninteractions(0);
            
            for (typename AssignedParticleContainerType::const_iterator pIter=assignedParticles.begin();
                 /*                                                  */ pIter!=assignedParticles.end();
                 /*                                                */ ++pIter)
            {
                // loop over neighbor cells of the current particle
                //                for (typename SpatialCellType::CellMapType::const_iterator cIter=(*pIter)->neighborCellsBegin();cIter!=(*pIter)->neighborCellsEnd();++cIter)
                for (typename SpatialCellType::CellMapType::const_iterator cIter=pIter->second->neighborCellsBegin();cIter!=pIter->second->neighborCellsEnd();++cIter)
                    
                {
                    // loop over particles in the neighbor cell
                    for(typename SpatialCellType::ParticleContainerType::const_iterator qIter=cIter->second->particleBegin();qIter!=cIter->second->particleEnd();++qIter) // loop over neighbor particles
                    {
                        //                        const typename ParticleContainerType::iterator nIter(this->find((*qIter)->sID));
                        //                        InteractionType(**pIter,*nIter->second);  // the constructor of InteractionType actually computes the binary interaction between *iterI and *iterJ
                        //                        InteractionType(**pIter,**qIter);  // the constructor of InteractionType actually computes the binary interaction between *iterI and *iterJ
                        InteractionType(*pIter->second,**qIter);  // the constructor of InteractionType actually computes the binary interaction between *iterI and *iterJ
                        
                        //                        ++Ninteractions;
                    }
                }
            }
            
            
            //mpiIDoffsetVector
            std::vector<int> interactionSizeVector(this->mpiProcs());
            std::vector<int> interactionRankOffsetVector(this->mpiProcs());
            
            for (unsigned int k=0;k<interactionSizeVector.size();++k)
            {
                // Store the number
                interactionSizeVector[k]=particleSizeVector[k]*InteractionType::DataPerParticle;
                interactionRankOffsetVector[k]=mpiIDoffsetVector[k]*InteractionType::DataPerParticle;
                
            }
            
            //! Syncronize InteractionType::resultVector among processors
            MPI_Allgatherv(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,
                           &InteractionType::resultVector[0],&interactionSizeVector[0],&interactionRankOffsetVector[0],MPI_DOUBLE,MPI_COMM_WORLD);
            
            if (this->mpiRank() == 0)
            {
                std::cout<<std::scientific<<" ["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]."<<std::endl;
            }
            
        }
        
        /**********************************************************************/
        template <typename InteractionType>
        typename InteractionType::ResultType getInteractionResult(const int& pID)
        {
            return InteractionType::getResult(*(this->find(pID)->second));
        }
        
        /**********************************************************************/
        template <typename InteractionType>
        void resetInteraction()
        {/*! Compute full particle interaction according to InteractionType
          */
            for (typename ParticleContainerType::iterator iterJ=this->begin();iterJ!=this->end();++iterJ) // FOR SYMMETRIC INTERACTION iterJ STARTS AT iterI
            {
                InteractionType::reset(*iterJ->second);  // the constructor of InteractionType actually computes the binary interaction between *iterI and *iterJ
            }
        }
        
        /**********************************************************************/
        template <class T>
        friend T& operator<< (T& os, const ParticleContainerType& pS)
        {/*! Operator << use ParticleType specific << operator
          */
//            int temp();
//            MPI_Comm_rank(MPI_COMM_WORLD,&mpiRank());
            //            std::cout<<"mpiRank()="<<mpiRank()<<std::endl;
            if (ModelMPIbase::mpiRank()==0)
            {
                for (typename ParticleContainerType::const_iterator pIter=pS.begin();pIter!=pS.end();++pIter)
                {
                    os<<(*pIter->second)<<std::endl;
                    }
                    }
                    MPI_Barrier(MPI_COMM_WORLD);
                    return os;
                    }
                    
                    };
                    
                    
                    } // end namespace
#endif
                    
                    
                    
                    
                    //        /**********************************************************************/
                    //        template <typename InteractionType>
                    //        void computeInteractionByCell()
                    //        {
                    //            for (typename AssignedCellContainerType::const_iterator cIter=assignedCells.begin();cIter!=assignedCells.end();++cIter)
                    //            { // loop over cells assigned to this process
                    //                for (typename SpatialCellType::ParticleContainerType::const_iterator pIter=(*cIter)->particleBegin();pIter!=(*cIter)->particleEnd();++pIter) // loop over particles in the current cell
                    //                {
                    //                    for (typename SpatialCellType::CellMapType::const_iterator nIter=(*cIter)->neighborCellsBegin();nIter!=(*cIter)->neighborCellsEnd();++nIter) // loop over neighbor cells
                    //                    {
                    //                        for(typename SpatialCellType::ParticleContainerType::const_iterator qIter=nIter->second->particleBegin();qIter!=nIter->second->particleEnd();++qIter) // loop over neighbor particles
                    //                        {
                    //                            //                            InteractionType(**pIter,**qIter);  // the constructor of InteractionType actually computes the binary interaction between *iterI and *iterJ
                    //                        }
                    //                    }
                    //                }
                    //            }
                    //        }
                    
                    //        /**********************************************************************/
                    //        template <typename InteractionType>
                    //        void computeInteractionByParticle()
                    //        {
                    //            // loop over particles assigned to this process
                    //            for (typename AssignedParticleContainerType::const_iterator pIter=assignedParticles.begin();pIter!=assignedParticles.end();++pIter)
                    //            {
                    //                // loop over neighbor cells of the current particle
                    //                for (typename SpatialCellType::CellMapType::const_iterator cIter=(*pIter)->neighborCellsBegin();cIter!=(*pIter)->neighborCellsEnd();++cIter)
                    //                {
                    //                    // loop over particles in the neighbor cell
                    //                    for(typename SpatialCellType::ParticleContainerType::const_iterator qIter=cIter->second->particleBegin();qIter!=cIter->second->particleEnd();++qIter) // loop over neighbor particles
                    //                    {
                    //                        //                            InteractionType(**pIter,**qIter);  // the constructor of InteractionType actually computes the binary interaction between *iterI and *iterJ
                    //                    }
                    //                }
                    //            }
                    //        }
                    
                    
                    
                    /* cells ****************************************/
                    //        SpatialCellObserverType cells()
                    //        {
                    //            return SpatialCellObserverType();
                    //        }
                    
                    //        /*****************************************/
                    //        template<char P='P',bool autodelete=true>
                    //        void MPIoutput()
                    //        {
                    //            //           int mpiRank();
                    //            //           MPI_Comm_rank(MPI_COMM_WORLD,&mpiRank());
                    //            if (this->mpiRank() == 0)
                    //                //if (PilMPI::mpiRank() == 0)
                    //            {
                    //                SequentialOutputFile<P,autodelete> pFile;
                    //                //                pFile<<*this<<std::endl;
                    //                pFile<<this->particles();
                    //
                    //
                    //                //                SpatialCellObserverType sco;
                    //                //                sco.end();
                    //
                    //                //                SpatialCellObserver<_ParticleType::SpatialCellType,_ParticleType::dim> sco;
                    //
                    //                SequentialOutputFile<'C',1> cFile;
                    //                //                cFile<<cells()<<std::endl;
                    //                //                cFile<<this->cells(); // this gives strange compilation error, so use the following
                    //                for (typename SpatialCellObserverType::CellMapType::const_iterator cIter =SpatialCellObserverType::cellBegin();
                    //                     /*                                                         */ cIter!=SpatialCellObserverType::cellEnd();
                    //                     /*                                                         */ cIter++)
                    //                {
                    //                    //os<<sC.cellID.transpose()<<" "<< sC.size()<<" "<<sC.neighborSize()<< " "<<sC.assignedRank;
                    //                    cFile<<cIter->second->cellID.transpose()<<" "<< cIter->second->size()<<" "<<cIter->second->neighborSize()<< " "<<cIter->second->assignedRank<<std::endl;
                    //                }
                    //
                    //
                    //
                    //                //                cFile<<sco<<std::endl;
                    //
                    //            }
                    //        }
