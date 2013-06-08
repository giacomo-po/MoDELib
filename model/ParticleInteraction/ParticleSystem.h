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






#include <memory> // for auto_ptr
#include <utility> // for std::pair
#include <map>
#include <vector>
#include <deque>
#include <time.h>
#include <iomanip> // std::scientific

#include <Eigen/Core>

#include <model/ParticleInteraction/ParticleSystemBase.h>
#ifdef _MODEL_MPI_
#include <model/ParticleInteraction/ParticleSystemParallel.h>
#else
#include <model/ParticleInteraction/ParticleSystemSerial.h>
#endif

#include <model/ParticleInteraction/SystemProperties.h>

//#include <model/SpaceDecomposition/SpatialCellObserver.h>
#include <model/Utilities/SequentialOutputFile.h>
#include <model/Utilities/CompareVectorsByComponent.h>


namespace model {
    
    
    
    
    /**************************************************************************/
    /**************************************************************************/
    /*! \brief Serial Version
     */
    template <typename _ParticleType, typename UserSystemProperties = SystemProperties<> >
    struct ParticleSystem :
#ifdef _MODEL_MPI_
    /* inheritance */  public ParticleSystemParallel<_ParticleType,UserSystemProperties>
    {
        
        /*****************************************/
        ParticleSystem(int argc, char* argv[], const double& cellSize=1.0) :
        /* init list */  ParticleSystemParallel<_ParticleType,UserSystemProperties>(argc,argv,cellSize)
        {/*!
          */
        }
        
    };
#else
    /* inheritance */  public ParticleSystemSerial<_ParticleType,UserSystemProperties>
    {
        
        /*****************************************/
        ParticleSystem(int argc, char* argv[], const double& cellSize=1.0) :
        /* init list */  ParticleSystemSerial<_ParticleType,UserSystemProperties>(argc,argv,cellSize)
        {/*!
          */
        }
    };
#endif
    
    
    
    
    
    
    
    
    
    
} // end namespace
#endif

