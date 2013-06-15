/* This file is part of PIL, the Particle Interaction Library.
 *
 * Copyright (C) 2012 by Giacomo Po <gpo@ucla.edu>
 * Copyright (C) 2012 by Shao-Ching Huang <sch@ucla.edu>
 * Copyright (C) 2012 by Tajendra Singh <tvsingh@ucla.edu>
 * Copyright (C) 2012 by Tamer Crosby <tcrosby@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _MODELMPIBASE_h_
#define _MODELMPIBASE_h_

#include <assert.h>
#include <mpi.h>
#include <model/Utilities/TerminalColors.h>

namespace model {
    
    
    
    class ModelMPIbase
    {
        
//        /*****************************************/
//        static int getMPIrank()
//        {
//            //           assert(mpiInitialized && "MPI not initialized.");
//            int mpiRank_temp;
//            MPI_Comm_rank(MPI_COMM_WORLD,&mpiRank_temp);
//            return mpiRank_temp;
//        }
//        
//        /*****************************************/
//        static int getMPInProcs()
//        {
//            int nProcs_temp;
//            MPI_Comm_size(MPI_COMM_WORLD,&nProcs_temp);
//            return nProcs_temp;
//        }
        

        
        
        static int _mpiInitialized;
        static int _mpiRank;
        static int _mpiProcs;

        
    public:
        
        
        /* Constructor */
        ModelMPIbase(int argc, char* argv[])
//        :
//        /* init list */ mpiInitialized(initMPI(argc,argv)),
//        /* init list */ mpiRank(getMPIrank()),
//        /* init list */ mpiProcs(getMPInProcs())
        {
//            std::cout<<greenBoldColor<<"MPI process "<<mpiRank<<" of "<<mpiProcs
//            /*     */<<", mpiInitialized="<<mpiInitialized<<defaultColor<<std::endl;
            
//            assert(mpiInitialized && "MPI not initialized.");
            
            initMPI(argc,argv);
            
        }
        
        ModelMPIbase()
        {
        }
        
        static void initMPI(int argc, char* argv[])
        {
//            int temp(0);
            MPI_Initialized(&_mpiInitialized);
            if (!_mpiInitialized)
            {
                MPI_Init(&argc,&argv);

                // Verify
                MPI_Initialized(&_mpiInitialized);
                assert(_mpiInitialized && "MPI not initialized.");

                MPI_Comm_rank(MPI_COMM_WORLD,&_mpiRank);
                MPI_Comm_size(MPI_COMM_WORLD,&_mpiProcs);
                std::cout<<greenBoldColor<<"MPI process "<<_mpiRank<<" of "<<_mpiProcs
                /*     */<<", mpiInitialized="<<_mpiInitialized<<defaultColor<<std::endl;

            }
//            // Verify
//            MPI_Initialized(&_mpiInitialized);
//            assert(_mpiInitialized && "MPI not initialized.");
 //           return temp;
        }
        
        
        static const int& mpiRank() {return _mpiRank;}
        static const int& mpiProcs(){return _mpiProcs;}
        static const int& mpiInitialized(){return _mpiInitialized;}
        
        /* Destructor */
        ~ModelMPIbase()
        {
            MPI_Finalize();
        }
        
    };
    
    // static data
    int ModelMPIbase::_mpiInitialized=0;
    int ModelMPIbase::_mpiRank=0;
    int ModelMPIbase::_mpiProcs=1;
    
} // end namespace
#endif