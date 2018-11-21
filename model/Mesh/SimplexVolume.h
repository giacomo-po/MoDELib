/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SimplexVolume_H_
#define model_SimplexVolume_H_

#include <SimplexTraits.h>
#include <CTM.h>

namespace model
{
	
	
	/**************************************************************************/
	/**************************************************************************/	
	template<short int dim, short int order>
	struct SimplexVolume
    {
        
        /**********************************************************************/
        static Eigen::Matrix<double,dim,SimplexTraits<dim,order>::nVertices-1>  edgeMatrix(const Eigen::Matrix<double,dim,SimplexTraits<dim,order>::nVertices>& posM)
        {            
            Eigen::Matrix<double,dim,SimplexTraits<dim,order>::nVertices-1> temp=posM.template block<dim,SimplexTraits<dim,order>::nVertices-1>(0,1);
            temp.colwise() -= posM.col(0);
            return temp;
        }
        
        /**********************************************************************/
        static double volume(const Eigen::Matrix<double,dim,SimplexTraits<dim,order>::nVertices>& posM)
        {
            const Eigen::Matrix<double,dim,SimplexTraits<dim,order>::nVertices-1> B=edgeMatrix(posM);
            return sqrt((B.transpose()*B).determinant())/CTM::factorial(SimplexTraits<dim,order>::nVertices-1);
        }
        
	};
    
    /**************************************************************************/
    /**************************************************************************/
    template<short int dim>
    struct SimplexVolume<dim,0>
    {
        /**********************************************************************/
        static double volume(const Eigen::Matrix<double,dim,SimplexTraits<dim,0>::nVertices>&)
        {
            return 0.0;
        }
        
    };
    
}	// close namespace
#endif
