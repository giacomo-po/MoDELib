/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SimplexReader_H_
#define model_SimplexReader_H_

#include <iostream>
#include <Eigen/Dense>
#include <model/Mesh/SimplexTraits.h>
//#include <model/IO/VertexReader.h>
#include <model/IO/IDreader.h>

namespace model
{
    
    
    /**************************************************************************/
    /**************************************************************************/
    template<int dim>
    struct SimplexReader
    {
        
        //        static VertexReader<'N',dim+1,double> nodeReader;
        static IDreader<'N',1,dim,double> nodeReader;
        
        
        /**********************************************************************/
        static Eigen::Matrix<double,dim,1> get_P0(const typename SimplexTraits<dim,0>::SimplexIDType& xID)
        {
            //            const typename VertexReader<'N',dim+1,double>::const_iterator nIter(nodeReader.find((xID)(0)));
            const typename IDreader<'N',1,dim,double>::const_iterator nIter(nodeReader.find((xID)(0)));
            assert((nIter!=nodeReader.end()) && "MESH VERTEX NOT FOUND IN N/N_x.txt.");
            return Eigen::Map<const Eigen::Matrix<double,dim,1>>(nIter->second.data());
        }
        
    };
    
    template<int dim>
    //    VertexReader<'N',dim+1,double> SimplexReader<dim>::nodeReader;
    IDreader<'N',1,dim,double> SimplexReader<dim>::nodeReader;
    
    
}	// close namespace
#endif
