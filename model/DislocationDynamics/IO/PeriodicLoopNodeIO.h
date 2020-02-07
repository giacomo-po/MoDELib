/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PeriodicLoopNodeIO_H_
#define model_PeriodicLoopNodeIO_H_

#include <tuple>
#include <iomanip>

namespace model
{
    
    template<short unsigned int dim>
    struct PeriodicLoopNodeIO
    {
        
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef Eigen::Matrix<double,dim-1,1> VectorLowerDim;

        size_t periodicLoopID;          // sID
        size_t sID;          // sID
        VectorDim Pg;
        VectorLowerDim Pl;

        

        /**********************************************************************/
        template<typename PeriodicLoopNodeType>
        PeriodicLoopNodeIO(const PeriodicLoopNodeType& node) :
        /* init */ periodicLoopID(node.periodicLoop.sID)
        /* init */,sID(node.sID)
        /* init */,Pg(node.periodicLoop.periodicGlidePlane->getGlobalPosition(node.P))
        /* init */,Pl(node.P)
        {
        }
        
        /**********************************************************************/
        PeriodicLoopNodeIO(const size_t& pID_in,
                          const size_t& sID_in,
                           const VectorDim& Pg_in,
                           const VectorLowerDim& Pl_in) :
        /* init */ periodicLoopID(pID_in)
        /* init */,sID(sID_in)
        /* init */,Pg(Pg_in)
        /* init */,Pl(Pl_in)
        {
            
        }
        
        /**********************************************************************/
        PeriodicLoopNodeIO() :
        /* init */ periodicLoopID(0)
        /* init */,sID(0)
        /* init */,Pg(VectorDim::Zero())
        /* init */,Pl(VectorLowerDim::Zero())
        {
            
            
        }
        
        /**********************************************************************/
        PeriodicLoopNodeIO(std::stringstream& ss) :
        /* init */ periodicLoopID(0)
        /* init */,sID(0)
        /* init */,Pg(VectorDim::Zero())
        /* init */,Pl(VectorLowerDim::Zero())
        {
            ss>>periodicLoopID;
            ss>>sID;
            for(int d=0;d<dim;++d)
            {
                ss>>Pg(d);
            }
            for(int d=0;d<dim-1;++d)
            {
                ss>>Pl(d);
            }
        }

        
        /**********************************************************************/
        template <class T>
        friend T& operator << (T& os, const PeriodicLoopNodeIO<dim>& ds)
        {
            os  << ds.periodicLoopID<<" "
            /**/<< ds.sID<<" "
            /**/<< std::setprecision(15)<<std::scientific<<ds.Pg.transpose()<<" "
            /**/<< std::setprecision(15)<<std::scientific<<ds.Pl.transpose();
            return os;
        }
        
	};
	
}
#endif

