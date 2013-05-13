/* This file is part of PIL, the Particle Interaction Library.
 *
 * Copyright (C) 2012 by Giacomo Po <gpo@ucla.edu>
 * Copyright (C) 2012 by Shao-Ching Huang <sch@ucla.edu>
 * Copyright (C) 2012 by Tajendra Singh <tvsingh@ucla.edu>
 * Copyright (C) 2012 by Tamer Crosby <tcrosby@ucla.edu>
 *
 * PIL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _PilMPI_h
#define _PilMPI_h

#include <mpi.h>

//#include <metis.h> // partitioner

//#include <memory> // for auto_ptr
//#include <utility> // for std::pair
//#include <map> 
//#include <vector>
//#include <deque> 
//#include <boost/ptr_container/ptr_map.hpp> // TO BE CHANGED WITH ACTUAL MPI IMPLEMENTATION
//
//#include <Eigen/Core>
//
//#include <pil/SystemProperties.h>
//#include <pil/SpatialCells/SpatialCellObserver.h>
//
//
//#include <pil/Utilities/SequentialOutputFile.h>


namespace pil {
    
    //template <typename _ParticleType, typename UserSystemProperties = SystemProperties<> >
    class PilMPI
    {
        
        /*****************************************/
        static int getMPIrank()
        {
            int mpiRank_temp;
            MPI_Comm_rank(MPI_COMM_WORLD,&mpiRank_temp);
            return mpiRank_temp;
        }
        
        /*****************************************/
        static int getMPInProcs()
        {
            int nProcs_temp;
            MPI_Comm_size(MPI_COMM_WORLD,&nProcs_temp);
            return nProcs_temp;
        }
        
        
    public:
        const int mpiRank;
        const int nProcs;
        
        /* Constructor */
        PilMPI() :
        /* init list */ mpiRank(getMPIrank()),
        /* init list */ nProcs(getMPInProcs())
        {
        
        }
        
    };
    
    
//	const int PilMPI::mpiRank=PilMPI::getMPIrank();
//	const int PilMPI::nProcs=PilMPI::getMPInProcs();

                
} // end namespace
#endif
