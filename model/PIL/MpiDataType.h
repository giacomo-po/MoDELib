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

#ifndef _pil_MpiDataType_h
#define _pil_MpiDataType_h

#include <mpi.h>


namespace pil {
    
    
    template <typename DataType>
    struct MpiDataType
    {

    };

    template <>
    struct MpiDataType<double>
    {
        //enum{DataType=MPI_DOUBLE};
        //typedef MPI_DOUBLE DataType;
//#define DataType MPI_DOUBLE
        MPI::Datatype DataType = MPI::DOUBLE;
    };
    
    template <>
    struct MpiDataType<int>
    {
//        enum{DataType=MPI_INT};
//#define DataType MPI_INT

    };


    
} // end namespace
#endif
