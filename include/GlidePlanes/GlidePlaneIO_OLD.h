/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_GlidePlaneIO_H_
#define model_GlidePlaneIO_H_

#include <tuple>
#include <iomanip>
#include <Eigen/Dense>
#include <GlidePlaneKey.h>

namespace model
{
    
    template<short unsigned int dim>
    struct GlidePlaneIO
    {
        
        typedef typename GlidePlaneKey<dim>::LongIntType LongIntType;
        
        typedef Eigen::Matrix<LongIntType,dim,1> VectorDimI;
        
        VectorDimI r;
        LongIntType h;
        LongIntType latticeID;
        
        /**********************************************************************/
        GlidePlaneIO(const GlidePlaneKey<dim>& key) :
        /* init */ r(key.reciprocalDirectionComponents()),
        /* init */ h(key.planeIndex()),
        /* init */ latticeID(key.latticeID())
        {
        }
                
        /**********************************************************************/
        GlidePlaneIO() :
        /* init */ r(VectorDimI::Zero()),
        /* init */ h(0),
        /* init */ latticeID(0)
        {
        }

        /**********************************************************************/
        GlidePlaneIO(std::stringstream& ss) :
        /* init */ r(VectorDimI::Zero()),
        /* init */ h(0),
        /* init */ latticeID(0)
        {
            for(int d=0;d<dim;++d)
            {
                ss>>r(d);
            }
            ss>>h;
            ss>>latticeID;
        }
        
        /**********************************************************************/
        template <class T>
        friend T& operator << (T& os, const GlidePlaneIO<dim>& ds)
        {
            os  << ds.r.transpose()<<"\t"
            /**/<<ds.h<<"\t"
            /**/<<ds.latticeID;
            return os;
        }
        
	};
	
}
#endif
