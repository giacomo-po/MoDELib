/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2011 by Benjamin Ramirez<ramirezbrf@gmail.com>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


	#include "model/Geometry/Splines/SplineNodeBase_corder0.h"

	private:
		
		//typedef Eigen::Matrix<double, dim, 1>	VectorDim;
		typedef Eigen::Matrix<double, dim, dim> MatrixDim;
		
		VectorDim T;
		


		public:

		MatrixDim prjM;	//! the projection matrix. THIS SHOULD BE PRIVATE


		const VectorDim & get_T() const {
			return T;
		}
		
