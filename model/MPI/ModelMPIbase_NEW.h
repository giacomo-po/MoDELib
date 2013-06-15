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

#include <memory> // std::shared_ptr
#include <set>
#include <model/MPI/MPIinit.h>

namespace model {
    
    class ModelMPIbase
    {

        std::shared_ptr<MPIinit> get_SharedPtr(int argc, char* argv[])
        {
            if (thisSet.empty())
            {
                return std::shared_ptr<MPIinit>(new MPIinit(argc,argv));
            }
            else
            {
                return (*thisSet.begin())->p_MPIinit;
            }
            
        }
        
        static std::set<const ModelMPIbase*> thisSet;
        
    public:
        
        const std::shared_ptr<MPIinit> p_MPIinit;
        const int& mpiRank;
        const int& mpiProcs;
        
        /* Constructor */
        ModelMPIbase(int argc, char* argv[]) :
        /* init list */ p_MPIinit(get_SharedPtr(argc,argv)),
        /* init list */ mpiRank(p_MPIinit->mpiRank),
        /* init list */ mpiProcs(p_MPIinit->mpiProcs)
        {
            assert(thisSet.insert(this).second && "CANNOT INSERT THIS");   
        }
        
    };
    
    // static data
    std::set<const ModelMPIbase*> ModelMPIbase::thisSet;
    
} // end namespace
#endif
