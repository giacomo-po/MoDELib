/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2011 by Benjamin Ramirez<ramirezbrf@gmail.com>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

	#include "model/Geometry/Splines/SplineNodeBase_corder1.h"
		
	protected:
		
		//typedef Eigen::Matrix<double, dim, 1> VectorDim;
		VectorDim K;
		
	public:
		
		virtual	const VectorDim & get_K() const{
			return K;
		}
		
	
