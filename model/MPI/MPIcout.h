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

#ifndef model_MPIcout_h_
#define model_MPIcout_h_

#include <assert.h>
#include <mpi.h>
#include <model/Utilities/TerminalColors.h>

namespace model {
    
    struct MPIcout : public std::ostream,
    public ModelMPIbase
    {
        typedef std::basic_ostream<char, std::char_traits<char> > StlEndl_IO;
        // define StlEndl as a pointer-to-function taking and returning a reference to StlEndl_IO
        typedef StlEndl_IO& (*StlEndl)(StlEndl_IO&);

#ifdef _MODEL_MPI_
//        template<typename T>
//        std::ostream& operator<<(const T& t)
//        {
//            if (this->mpiRank()==0)
//            {
//                return std::cout<<t;
//            }
//            else
//            {
//                return std::cout;
//            }
//        }
        template<typename T>
        MPIcout& operator<<(const T& t)
        {
            if (this->mpiRank()==0)
            {
                std::cout<<t;
                return *this;
            }
            else
            {
                return *this;
            }
        }
#else
//        template<typename T>
//        std::ostream& operator<<(const T& t)
//        {
//            return std::cout<<t;
//        }
        template<typename T>
        MPIcout& operator<<(const T& t)
        {
            std::cout<<t;
            return *this;
        }
#endif
//                /**********************************************************************/
//                std::ostream& operator<<(StlEndl manip)
//                {/*! Overload << for Std::endl
//                  */
//#ifdef _MODEL_MPI_
//                    //            int mpiRank;
//                    //            MPI_Comm_rank(MPI_COMM_WORLD,&mpiRank);
//                    if(mpiRank==0)
//                    {
//                        manip(std::cout);
//                    }
//#else
//                    manip(std::cout);
//                    
//#endif            
//                    return std::cout;
//                }

                /**********************************************************************/
                MPIcout& operator<<(StlEndl manip)
                {/*! Overload << for Std::endl
                  */
#ifdef _MODEL_MPI_
                    if(mpiRank==0)
                    {
                        manip(*this);
                    }
#else
                    manip(*this);
                    
#endif
                    return *this;
                }
                
                };
                
    MPIcout cout;
                
} // end namespace
#endif
