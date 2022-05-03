/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DISLOCATIONLOCALREFERENCE_H
#define model_DISLOCATIONLOCALREFERENCE_H

#include <math.h>
#include <float.h>
#include <Eigen/Dense>

namespace model
{

	template <short unsigned int dim>
	struct DislocationLocalReference
    {
	
		typedef Eigen::Matrix<double,dim,1>   VectorDimD;
		typedef Eigen::Matrix<double,dim,dim> MatrixDimD;

		static MatrixDimD global2local(const VectorDimD& C, const VectorDimD& N);
	};

}
#endif


