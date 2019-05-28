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

#ifndef model_MPIcout_h_
#define model_MPIcout_h_

#include <iostream>

#ifdef _MODEL_MPI_
#include <mpi.h>
#include <ModelMPIbase.h>
#endif

namespace model
{
    
    struct MPIcout
    {
        typedef std::basic_ostream<char, std::char_traits<char> > StlEndl_IO;
        // define StlEndl as a pointer-to-function taking and returning a reference to StlEndl_IO
        typedef StlEndl_IO& (*StlEndl)(StlEndl_IO&);
        
        
        /**********************************************************************/
        MPIcout(const MPIcout&) :
        /* base init*/ os(std::cout)
        {/*! Copy constructor is private (no copying allowed)
          */
        }
        
        /**********************************************************************/
        MPIcout& operator=(const MPIcout&)
        {/*! Assignment operator is private (no assignment allowed)
          */
            return *this;
        }
        
        
        std::ostream& os;
        
    public:
        /**********************************************************************/
        MPIcout() :
        /* base init*/ os(std::cout)
        {/*!
          */
        }
        
        /**********************************************************************/
        template<typename T>
        MPIcout& operator<<(const T& t)
        {
#ifdef _MODEL_MPI_
            
            if (ModelMPIbase::mpiRank()==0)
            {
                os<<t;
            }
#else
            os<<t;
#endif
            return *this;
        }
        
        /**********************************************************************/
        MPIcout& operator<<(StlEndl manip)
        {/*! Overload << for Std::endl
          */
#ifdef _MODEL_MPI_
            if(ModelMPIbase::mpiRank()==0)
            {
                manip(os);
            }
#else
            manip(os);
            
#endif
            return *this;
        }
        
    };
    
    // declare object cout
    static MPIcout cout;
    
} // end namespace
#endif
