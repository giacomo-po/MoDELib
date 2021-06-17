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

#ifndef _MODELMPIBASE_h_
#define _MODELMPIBASE_h_

#include <assert.h>
#include <string>
#include <mpi.h>
#include <TerminalColors.h>

namespace model {
    
    class ModelMPIbase
    {

        static int _mpiInitialized;
        static int _mpiRank;
        static int _mpiProcs;
        static std::string _processorName;

        
    public:
        
        /**********************************************************************/
        ModelMPIbase(int argc, char* argv[])
        {/*! Constructor 
          */
            initMPI(argc,argv);
        }
        
        /**********************************************************************/
        ModelMPIbase()
        {
            
        }
        
        /**********************************************************************/
        ~ModelMPIbase()
        {/*! Destructor
          */
            MPI_Finalize();
        }
        
        /**********************************************************************/
        static void initMPI(int argc, char* argv[])
        {
            MPI_Initialized(&_mpiInitialized);
            if (!_mpiInitialized)
            {
                MPI_Init(&argc,&argv);

                // Verify
                MPI_Initialized(&_mpiInitialized);
                assert(_mpiInitialized && "MPI not initialized.");
                
                
                MPI_Comm_rank(MPI_COMM_WORLD,&_mpiRank);
                MPI_Comm_size(MPI_COMM_WORLD,&_mpiProcs);
                
                int tempLen;
                char tempName[MPI_MAX_PROCESSOR_NAME];
                MPI_Get_processor_name(tempName, &tempLen);
                _processorName=tempName;
                
                std::cout<<greenBoldColor<<"MPI process "<<_mpiRank<<" of "<<_mpiProcs
                /*     */<<" (on node "<<_processorName<<"), mpiInitialized="<<_mpiInitialized<<defaultColor<<std::endl;

                MPI_Barrier(MPI_COMM_WORLD);

            }
        }
        
        
        static const int& mpiRank() {return _mpiRank;}
        static const int& mpiProcs(){return _mpiProcs;}
        static const int& mpiInitialized(){return _mpiInitialized;}
        static const std::string& processorName(){return _processorName;}
        
    };
    
    // static data
    int ModelMPIbase::_mpiInitialized=0;
    int ModelMPIbase::_mpiRank=0;
    int ModelMPIbase::_mpiProcs=1;
    std::string ModelMPIbase::_processorName="NO NAME ASSIGNED YET";
    
} // end namespace
#endif
