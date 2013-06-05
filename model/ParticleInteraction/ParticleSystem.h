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

#ifndef _MODEL_ParticleSystem_h_
#define _MODEL_ParticleSystem_h_


#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#endif

#ifdef _MODEL_MPI_
#include <model/MPI/ModelMPIbase.h>
#include <metis.h> // partitioner
#endif


#include <memory> // for auto_ptr
#include <utility> // for std::pair
#include <map> 
#include <vector>
#include <deque> 
#include <boost/ptr_container/ptr_map.hpp> // TO BE CHANGED WITH ACTUAL MPI IMPLEMENTATION
#include <time.h> 
#include <iomanip> // std::scientific

#include <Eigen/Core>

#include <model/ParticleInteraction/SystemProperties.h>
#include <model/SpaceDecomposition/SpatialCellObserver.h>
#include <model/Utilities/SequentialOutputFile.h>
#include <model/Utilities/CompareVectorsByComponent.h>


namespace model {
    
    
    // TO DO: REMOVE cellSize from Constructor and simplify constructors

    template <typename _ParticleType, typename UserSystemProperties>
    class ParticleSystemBase :
    /* inheritance */  protected boost::ptr_map<size_t,_ParticleType>,
    /* inheritance */  private UserSystemProperties
    {
    
        size_t particleCounter;
        
    public:
        
        typedef _ParticleType ParticleType; // make ParticleType available outside the class
        typedef typename _ParticleType::PositionType PositionType;
        typedef typename _ParticleType::SpatialCellType SpatialCellType;
        typedef typename _ParticleType::SpatialCellObserverType SpatialCellObserverType;
        typedef boost::ptr_map<size_t,_ParticleType> ParticleContainerType;

        /**********************************************************************/
        ParticleSystemBase(const double& cellSize=1.0) :
        /* init list */ particleCounter(0)
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
//            // 2- Obtain the sID (StaticID) of the new particle
//            const size_t sID(pN->sID);

            const size_t sID(particleCounter);
            particleCounter++;
            
            // 3- Insert the particle in the ParticleSystem using sID as the key. Assert succesful insertion
            assert(this->insert(sID , pN ).second && "CANNOT INSERT PARTICLE IN PARTICLE SYSTEM.");
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
        
    };
    
    /**************************************************************************/
    /**************************************************************************/
    /*! \brief Serial Version
     */
    template <bool useMPI, typename _ParticleType, typename UserSystemProperties = SystemProperties<> >
    class ParticleSystem :
    /* inheritance */  public ParticleSystemBase<_ParticleType,UserSystemProperties>
    {
        typedef ParticleSystemBase<_ParticleType,UserSystemProperties> ParticleSystemBaseType;
        typedef _ParticleType ParticleType; // make ParticleType available outside the class
//        typedef typename ParticleSystemBaseType::PositionType PositionType;
        typedef typename ParticleSystemBaseType::SpatialCellType SpatialCellType;
        typedef typename ParticleSystemBaseType::ParticleContainerType ParticleContainerType;
        typedef typename ParticleSystemBaseType::SpatialCellObserverType SpatialCellObserverType;

        
        
        
        void loopOverNeighCells(typename ParticleContainerType::iterator& pIter)
        {
            for (typename SpatialCellType::CellMapType::const_iterator cIter=(*pIter)->neighborCellsBegin();cIter!=(*pIter)->neighborCellsEnd();++cIter)
            {
                // loop over particles in the neighbor cell
                for(typename SpatialCellType::ParticleContainerType::const_iterator qIter=cIter->second->particleBegin();qIter!=cIter->second->particleEnd();++qIter) // loop over neighbor particles
                {
                    InteractionType(**pIter,**qIter);  // the constructor of InteractionType actually computes the binary interaction between *iterI and *iterJ
                }
            }
        }
        
        
    public:
        /*****************************************/
        ParticleSystem(int argc, char* argv[], const double& cellSize=1.0) :
        /* init list */  ParticleSystemBaseType::ParticleSystemBase(cellSize)
        {/*! 
          */

        }
        
        
        /*****************************************/
        template <typename InteractionType>
        void computeInteraction()
        {
            
#ifdef _OPENMP
#pragma omp parallel for
            for (unsigned int k=0;k<this->size();++k)
            {
                typename ParticleContainerType::iterator pIter(this->begin()); //  data within a parallel region is private to each thread
                std::advance(pIter,k);
                loopOverNeighCells(pIter);
            }
#else
            for (typename ParticleContainerType::const_iterator pIter=this->begin();pIter!=this->end();++pIter)
            {
                loopOverNeighCells(pIter);
            }
#endif
                

        }
        
    };
    
    
#ifdef _MODEL_MPI_
    /**************************************************************************/
    /**************************************************************************/
    /*! \brief Parallel Version (MPI)
     */
    template <typename _ParticleType, typename UserSystemProperties>
    class ParticleSystem<true,_ParticleType,UserSystemProperties> :
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
        typedef std::deque<SpatialCellType*> AssignedCellContainerType;
        typedef std::deque<_ParticleType*>   AssignedParticleContainerType;
        
        //! The SpatialCell cells assigned to this MPI process
        AssignedCellContainerType assignedCells;
        AssignedParticleContainerType assignedParticles;
        
        
        bool useCellPartitioner;

        bool partitionIsValid;
        std::vector<int> mpiIDoffsetVector;
        std::vector<int> particleSizeVector;
        
        /*****************************************/
        void assignMpiIDs()
        {/*! Assigns each particle in the system a unique ID (SpatialParticle::mpiID)
          * such that particles are sorted by MPI process
          */
            
//            // Find Offset
////            mpiIDoffset=0;
//            int mpiIDoffset(0);
//            // Compute mpiIDoffset by comulative summing assignedParticleSize across the MPI processes
//            MPI_Exscan(&assignedParticleSize,&mpiIDoffset,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
//
//            // Collect assignedParticleSize in the particleSizeVector the MPI processes
//            MPI_Allgather(&assignedParticleSize,1,MPI_INT,&particleSizeVector[0],1,MPI_INT,MPI_COMM_WORLD);

            
            
  //          std::cout<<"Processor "<<this->mpiRank<<" :  mpiIDoffset="<<mpiIDoffset<<std::endl;
            
            int mpiIDlocal(0);
            
            
            if (useCellPartitioner)
            {
                // loop over cells assigned to this process
                for (typename AssignedCellContainerType::iterator cIter=assignedCells.begin();cIter!=assignedCells.end();++cIter)
                {
                    // loop over particles in the current cell
                    for (typename SpatialCellType::ParticleContainerType::iterator pIter=(*cIter)->particleBegin();pIter!=(*cIter)->particleEnd();++pIter) 
                    {
//                                            (*pIter)->mpiID=mpiIDlocal+mpiIDoffset;
//                        (*pIter)->mpiID=mpiIDlocal+mpiIDoffsetVector[this->mpiRank];
                        mpiIDlocal++;
                    }
                }
            }
            else
            {
                // loop over cells assigned to this process
                for (typename AssignedParticleContainerType::iterator pIter=assignedParticles.begin();pIter!=assignedParticles.end();++pIter)
                {
//                    (*pIter)->mpiID=mpiIDlocal+mpiIDoffsetVector[this->mpiRank];
//                                            (*pIter)->mpiID=mpiIDlocal+mpiIDoffset;
                        mpiIDlocal++;
                }
            }
            
        
        }
        
        /**********************************************************************/
        int partionCells()
        {
            typedef typename SpatialCellObserverType::CellIdType CellIdType;
            typedef std::map<CellIdType,const int,CompareVectorsByComponent<int,_ParticleType::dim> > ijk2idMapType;
            //            typedef std::map<Eigen::Matrix<int,_ParticleType::dim,1>,const int,CompareVectorsByComponent<int,_ParticleType::dim> > ijk2idMapType;
            ijk2idMapType  ijk2idMap;;
            typedef std::map<int,const CellIdType> id2ijkMapType;
            //            typedef std::map<int,Eigen::Matrix<int,_ParticleType::dim,1> > id2ijkMapType;
            id2ijkMapType id2ijkMap;;
            
            
            std::vector<int> vWeights(SpatialCellObserverType::size(),1); // vector containing the weight of each cell
            
            
            int id(0);
            for (typename SpatialCellObserverType::CellMapType::const_iterator cIter=SpatialCellObserverType::cellBegin();cIter!=SpatialCellObserverType::cellEnd();++cIter)
            {
                ijk2idMap.insert(std::pair<Eigen::Matrix<int,_ParticleType::dim,1>,const int>(cIter->second->cellID,id));
                id2ijkMap.insert(std::pair<int,Eigen::Matrix<int,_ParticleType::dim,1> >(id,cIter->second->cellID));
                vWeights[id]=cIter->second->n2Weight();
                id++;
            }
            
            
            
            // - Populate vector of neighbors and offsets
            std::vector<int> vNeighbors;
            std::vector<int> vOffsets;
            
            int offSet(0);
            for (typename SpatialCellObserverType::CellMapType::const_iterator cIter=SpatialCellObserverType::cellBegin();cIter!=SpatialCellObserverType::cellEnd();++cIter)
            {
                vOffsets.push_back(offSet);
                for (typename SpatialCellObserverType::CellMapType::const_iterator nIter=cIter->second->neighborCellsBegin();nIter!=cIter->second->neighborCellsEnd();++nIter)
                {
                    if (cIter->second->cellID!=nIter->second->cellID)
                    {
                        typename ijk2idMapType::const_iterator mIter(ijk2idMap.find(nIter->second->cellID.transpose()));
                        vNeighbors.push_back(mIter->second);
                        offSet++;
                    }
                }
            }
            vOffsets.push_back(offSet); // push back an extra offset (required by Metis)
            
            
            // Create a vector to store the output of the partitioner. Initialize at process=0
            std::vector<int> partitionerRankVector(SpatialCellObserverType::size(),0);

            // partition into the same number of MPI processes
            int numPart = this->mpiProcs;

            if (numPart > 1)
            {
                int metisOptions[METIS_NOPTIONS];
                int numVert = SpatialCellObserverType::size(); // number of nodes
                int ncon = 1;  // number of weights per node
                int objval;
                float ubvec = 1.001; // load imbalance tolerance
                METIS_SetDefaultOptions(metisOptions);
                METIS_PartGraphKway(&numVert,&ncon,&vOffsets[0],&vNeighbors[0],
                                    &vWeights[0],NULL,NULL,&numPart,NULL,
                                    &ubvec,metisOptions,&objval,&partitionerRankVector[0]);
            }
            
            
            
            int assignedParticleSize(0);

            // Assign ranks to cells
            for (unsigned int i=0; i<SpatialCellObserverType::size(); ++i)
            {
                typename id2ijkMapType::const_iterator idIter(id2ijkMap.find(i));
                const CellIdType cellID(idIter->second);
                
                typename SpatialCellObserverType::isCellType isC(SpatialCellObserverType::isCell(cellID));
                assert(isC.first && "CELL MUST EXIST!");
                isC.second->assignedRank=partitionerRankVector[i];
                
                if (this->mpiRank==partitionerRankVector[i])
                {
                    assignedCells.push_back(isC.second);
                    assignedParticleSize+=isC.second->size();
                }
            }
            
            
            // Evaluate performance
            
            float totalWeight(0.0);
            for (unsigned int w=0;w<vWeights.size();++w)
            {
                totalWeight+=vWeights[w];
            }
//            double optimumLoad=totalWeight/vWeights.size();

            
            std::vector<int> procWeights(this->mpiProcs,0); // vector containing the weight of each processor

            for (unsigned int w=0;w<partitionerRankVector.size();++w)
            {
                procWeights[partitionerRankVector[w]]+=vWeights[w];
            }
            
            double speedup(this->mpiProcs);
            for (unsigned int w=0;w<procWeights.size();++w)
            {
                float temp(totalWeight/procWeights[w]);
                if(speedup>temp)
                {
                    speedup=temp;
                }
//                procWeights[partitionerRankVector[w]]+=vWeights[w];
            }
            std::cout<<"speedUp="<<speedup<<std::endl;
            
            return assignedParticleSize;
        }
        
        /**********************************************************************/
        int partionParticles()
        {
            // Initialize the vector containing the weight of each particle
            std::vector<int> vWeights(this->size(),1);
            
            // Populate the vector containing the weight of each particle
            int kk(0);
            for (typename ParticleContainerType::const_iterator pIter=this->begin();pIter!=this->end();++pIter)
            {
                // The weight of the kk-th particle is the n2Weight() of its cell
                vWeights[kk]=pIter->second->pCell->n2Weight();
                kk++;
            }

            // Initialize METIS vectors
            std::vector<int> vNeighbors;
            std::vector<int> vOffsets;
            
            // Populate METIS vectors
            int offSet(0);
            int pp(0);
            for (typename ParticleContainerType::const_iterator pIter=this->begin();pIter!=this->end();++pIter)
            {
                vOffsets.push_back(offSet);
                vNeighbors.push_back(pp); // each particle is its only neighbor
                offSet++;
                pp++;
            }
            vOffsets.push_back(offSet); // push back an extra offset (required by Metis)

            
            // partition into the same number of MPI processes
            int numPart = this->mpiProcs;
            
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
            
            int assignedParticleSize(0);
            
            for (unsigned int i=0; i<this->size(); ++i)
            {                
                if (this->mpiRank==partitionerRankVector[i])
                {
                    typename ParticleContainerType::iterator pIter(this->begin());
                    std::advance(pIter,i);
                    assignedParticles.push_back(pIter->second);
                    assignedParticleSize+=1;
                }
            }
            
            return assignedParticleSize;
        }
        
        /*****************************************/
        template <typename InteractionType>
        void computeInteractionByCell()
        {
            for (typename AssignedCellContainerType::const_iterator cIter=assignedCells.begin();cIter!=assignedCells.end();++cIter){ // loop over cells assigned to this process
                for (typename SpatialCellType::ParticleContainerType::const_iterator pIter=(*cIter)->particleBegin();pIter!=(*cIter)->particleEnd();++pIter) // loop over particles in the current cell
                {
                    for (typename SpatialCellType::CellMapType::const_iterator nIter=(*cIter)->neighborCellsBegin();nIter!=(*cIter)->neighborCellsEnd();++nIter) // loop over neighbor cells
                    {
                        for(typename SpatialCellType::ParticleContainerType::const_iterator qIter=nIter->second->particleBegin();qIter!=nIter->second->particleEnd();++qIter) // loop over neighbor particles
                        {
                            //                            InteractionType(**pIter,**qIter);  // the constructor of InteractionType actually computes the binary interaction between *iterI and *iterJ
                        }
                    }
                }
            }
        }
        
        /*****************************************/
        template <typename InteractionType>
        void computeInteractionByParticle()
        {
            // loop over particles assigned to this process
            for (typename AssignedParticleContainerType::const_iterator pIter=assignedParticles.begin();pIter!=assignedParticles.end();++pIter)
            {
                // loop over neighbor cells of the current particle
                for (typename SpatialCellType::CellMapType::const_iterator cIter=(*pIter)->neighborCellsBegin();cIter!=(*pIter)->neighborCellsEnd();++cIter)
                {
                    // loop over particles in the neighbor cell
                    for(typename SpatialCellType::ParticleContainerType::const_iterator qIter=cIter->second->particleBegin();qIter!=cIter->second->particleEnd();++qIter) // loop over neighbor particles
                    {
                        //                            InteractionType(**pIter,**qIter);  // the constructor of InteractionType actually computes the binary interaction between *iterI and *iterJ
                    }
                }
            }
        }
        

        
    public:
        
//        static bool useCellPartitioner;
        
        
        
        /*****************************************/
        ParticleSystem(int argc, char* argv[], const double& cellSize=1.0) :
        /* init list */  ParticleSystemBaseType::ParticleSystemBase(cellSize),
        /* init list */  ModelMPIbase(argc,argv),
        /* init list */  useCellPartitioner(true),
        /* init list */  partitionIsValid(false)
        {
            particleSizeVector.resize(this->mpiProcs);
            mpiIDoffsetVector.resize(this->mpiProcs);
        }
        
        /**********************************************************************/
        template <typename ...AdditionalConstructorTypes>
        int addParticle(const typename ParticleSystemBaseType::PositionType& p, const AdditionalConstructorTypes&... args)
        {/*! Adds a ParticleType particle with position p and an arbitrary (variadic)
          * number of additional constructor parameters.
          */
            partitionIsValid=false;
            return ParticleSystemBaseType::addParticle(p,args...);
        }
        

        /**********************************************************************/
        float partionSystem(const bool& useCellPartitioner_in)
        {
            
            useCellPartitioner=useCellPartitioner_in;
            
//            const float t0(clock());
            
            assignedCells.clear();
            assignedParticles.clear();

        
            //! The number of ParticleType particles assigned to this MPI process
            int assignedParticleSize(0);

            if(useCellPartitioner)
            {
//                if (this->mpiRank == 0)
//                {
//                    std::cout<<"Partitioning system by cells... "<<std::flush;
//                }
                assignedParticleSize=partionCells();
            }
            else
            {
//                std::cout<<"Partitioning system by particles... "<<std::flush;
                assignedParticleSize=partionParticles();
            }

            // Collect assignedParticleSize from each MPI process into particleSizeVector 
            MPI_Allgather(&assignedParticleSize,1,MPI_INT,&particleSizeVector[0],1,MPI_INT,MPI_COMM_WORLD);
            
            // Check that all particles have been partitioned
            unsigned int totalAssignedParticle(0);
            for (unsigned int k=0;k<particleSizeVector.size();++k)
            {
//                if (this->mpiRank == 0)
//                {
//                    std::cout<<particleSizeVector[k]<<std::endl;
//                }
                totalAssignedParticle+=particleSizeVector[k];
            }
//            if (this->mpiRank == 0)
//            {
//                std::cout<<totalAssignedParticle<<" vs "<<this->size()<<std::endl;
//            }
            assert(totalAssignedParticle==this->size() && "SOMETHING WENT VERY WRONG");

            // Compute mpiIDoffset by comulative summing assignedParticleSize across the MPI processes
            int mpiIDoffset(0);
            MPI_Exscan(&assignedParticleSize,&mpiIDoffset,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
            
            // Collect mpiIDoffset from each MPI process into mpiIDoffsetVector
            MPI_Allgather(&mpiIDoffset,1,MPI_INT,&mpiIDoffsetVector[0],1,MPI_INT,MPI_COMM_WORLD);
            
            // Update mpiID for each particle
            assignMpiIDs();
            
            // Set flag "partitionIsValid" to true
            partitionIsValid=true;
            
//            if (this->mpiRank == 0)
//            {
//                std::cout<<std::scientific<<" ["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]."<<std::endl;
//            }

            
//            return HERE RETURN THE LOAD IMBALANCE;
            return 1.0;

        }
        

        
        /**********************************************************************/
        template <typename InteractionType>
        void computeInteraction()
        {/*! Compute nearest-neighbor particle interaction according to InteractionType
          */
            
            assert(partitionIsValid && "PARTITION IS NOT VALID, TRY CALLING ParticleSystem::partionSystem()");
            
            // Construct a vector to store the number of particles in each processor
//            std::vector<int> particleSizeVector(this->mpiProcs,0);
//            particleSizeVector.resize(this->mpiProcs);

            
//            assignMpiIDs(particleSizeVector);
            
            
//            InteractionType::resultVector.resize(3*this->size(),0.0);
            InteractionType::resize(3*this->size(),0.0);

            
            if (useCellPartitioner)
            {
                computeInteractionByCell<InteractionType>();
            
            }
            else
            {
                computeInteractionByParticle<InteractionType>();
            }
            
            

            //InteractionType::resultVerctor;
            
            //mpiIDoffsetVector
            std::vector<int> interactionSizeVector(this->mpiProcs);            
            std::vector<int> interactionRankOffsetVector(this->mpiProcs);
            
            for (unsigned int k=0;k<interactionSizeVector.size();++k)
            {
                // Store the number 
                interactionSizeVector[k]=particleSizeVector[k]*InteractionType::DataPerParticle;
                interactionRankOffsetVector[k]=mpiIDoffsetVector[k]*InteractionType::DataPerParticle;
                
            }
            
            
            MPI_Allgatherv(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,
                           &InteractionType::resultVector[0],&interactionSizeVector[0],&interactionRankOffsetVector[0],MPI_DOUBLE,MPI_COMM_WORLD);
//            &InteractionType::resultVector[0],&interactionSizeVector[0],&interactionRankOffsetVector[0],InteractionType::ResultType,MPI_COMM_WORLD);

            
//            if (this->mpiRank==0)
//            {
//                for (unsigned int k=0;k<InteractionType::resultVector.size()/3;++k)
//                {
//                    std::cout<<InteractionType::resultVector[3*k]<<" "<<InteractionType::resultVector[3*k+1]<<" "<<InteractionType::resultVector[3*k+2]<<"\n";
//                }
//            
//            }
            
            
        }
        
        
        /*****************************************/
        template <typename InteractionType>
        typename InteractionType::ResultType getInteractionResult(const int& pID)
        {
            return InteractionType::getResult(*(this->find(pID)->second));
        }
        
        /*****************************************/
        template <typename InteractionType>
        void resetInteraction()
        {/*! Compute full particle interaction according to InteractionType
          */
            for (typename ParticleContainerType::iterator iterJ=this->begin();iterJ!=this->end();++iterJ){ // FOR SYMMETRIC INTERACTION iterJ STARTS AT iterI
                InteractionType::reset(*iterJ->second);  // the constructor of InteractionType actually computes the binary interaction between *iterI and *iterJ
            }
        }
        

        
        
        /* cells ****************************************/
//        SpatialCellObserverType cells()
//        {
//            return SpatialCellObserverType();
//        }
        
        /*****************************************/
        template<char P='P',bool autodelete=true>
        void MPIoutput()
        {
 //           int mpiRank;
 //           MPI_Comm_rank(MPI_COMM_WORLD,&mpiRank);
            if (this->mpiRank == 0)
            //if (PilMPI::mpiRank == 0)
            {
                SequentialOutputFile<P,autodelete> pFile;
                pFile<<*this<<std::endl;
                
                
//                SpatialCellObserverType sco;
//                sco.end();
                
//                SpatialCellObserver<_ParticleType::SpatialCellType,_ParticleType::dim> sco;
                
                SequentialOutputFile<'C',1> cFile;
//                cFile<<cells()<<std::endl;
                
                for (typename SpatialCellObserverType::CellMapType::const_iterator cIter=SpatialCellObserverType::cellBegin();cIter!=SpatialCellObserverType::cellEnd();++cIter)
                {
                    //os<<sC.cellID.transpose()<<" "<< sC.size()<<" "<<sC.neighborSize()<< " "<<sC.assignedRank;
                    cFile<<cIter->second->cellID.transpose()<<" "<< cIter->second->size()<<" "<<cIter->second->neighborSize()<< " "<<cIter->second->assignedRank<<std::endl;
                }
                

                
//                cFile<<sco<<std::endl;

            }
        }
        
        /*****************************************/
        template <class T>
        friend T& operator<< (T& os, const ParticleSystem<true,ParticleType,UserSystemProperties>& pS)
        {/*! Operator << use ParticleType specific << operator
          */
            for (typename ParticleContainerType::const_iterator pIter=pS.begin();pIter!=pS.end();++pIter)
            {
                os<<(*pIter->second)<<std::endl;
                
            }
            return os;
        }
                
    };
        
                
//    template <typename _ParticleType, typename UserSystemProperties>
//    bool ParticleSystem<true,_ParticleType,UserSystemProperties>::useCellPartitioner=true;
#endif // ParticleSystem<1,_ParticleType,UserSystemProperties>

                
} // end namespace
#endif
