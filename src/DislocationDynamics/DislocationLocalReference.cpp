/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DISLOCATIONLOCALREFERENCE_cpp_
#define model_DISLOCATIONLOCALREFERENCE_cpp_

#include <DislocationLocalReference.h>

namespace model
{


	
		template <short unsigned int dim>
		typename DislocationLocalReference<dim>::MatrixDimD DislocationLocalReference<dim>::global2local(const VectorDimD& C, const VectorDimD& N)
        {
			assert(C.norm()>FLT_EPSILON);
			assert(N.norm()>FLT_EPSILON);
			assert(std::fabs(C.dot(N))<FLT_EPSILON);
			MatrixDimD R(MatrixDimD::Zero());
			R.col(0)=C.normalized();
			R.col(2)=N.normalized();
			R.col(1)= R.col(2).cross(R.col(0));
//			assert(R.determinant())
			return R.transpose();
//			return R;
		}

		template struct DislocationLocalReference<3>;

}
#endif


