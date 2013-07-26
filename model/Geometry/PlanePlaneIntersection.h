/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2011 by Benjamin Ramirez<ramirezbrf@gmail.com>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */



#ifndef model_PLANEPLANEINTERSECTION_H_
#define model_PLANEPLANEINTERSECTION_H_

#include <Eigen/Dense>
#include <float.h> // defines FLT_EPSILON

namespace model {
    
    template <short unsigned int dim>
    struct PlanePlaneIntersection
    {
        typedef Eigen::Matrix<double,dim,1>   VectorDim;
        enum{parallelPlanes=0,coincidentPlanes=1,incidentPlanes=2};
        
        /**********************************************************************/
		static int planePlaneType(const VectorDim & normal1,
                                  const VectorDim & normal2,
                                  const VectorDim & pt1,
                                  const VectorDim & pt2)
        {
            
			bool areParallelNormals(normal1.cross(normal2).norm()<FLT_EPSILON);
			bool areCoincidentPoints(std::fabs((pt1-pt2).dot(normal1))<FLT_EPSILON);
			
			int i;
			if(areParallelNormals && !areCoincidentPoints)
            {
				/*parallel planes: no intersection*/
				i= 0;
			}
            else if(areParallelNormals && areCoincidentPoints)
            {
				/*unique planes: planar intersection intersection*/
				i= 1;
			}
            else
            {
				/*angle planes: angular intersection*/
				i= 2;
			}
			
			return i;
		}
        
        
    };
	
	
	//////////////////////////////////////////////////////////////s
} // namespace model
#endif

