/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_CRYSTALSTRUCTURES_H_
#define model_CRYSTALSTRUCTURES_H_

#include <vector>
#include <Eigen/Dense>

namespace model {
    
    /**************************************************************************/
    /**************************************************************************/
    struct FCC {
        
        template <int dim>
        static std::vector<Eigen::Matrix<double,dim,1> > getPlaneNormals(){
            typedef Eigen::Matrix<double,dim,1> VectorDim;
            typedef std::vector<VectorDim> PlaneNormalContainerType;
            
            VectorDim  alpha(-1.0, 1.0,-1.0);
			VectorDim  beta( 1.0,-1.0,-1.0);
			VectorDim gamma(-1.0,-1.0, 1.0);
			VectorDim delta( 1.0, 1.0, 1.0);
            
            PlaneNormalContainerType temp;
            temp.push_back(alpha);
            temp.push_back(beta);
            temp.push_back(gamma);
            temp.push_back(delta);
            
            return temp;
        }

        
    };
    
    /**************************************************************************/
    /**************************************************************************/
    
//    template <>
//    class FCC<3> {
//        
//        enum{dim=3};
//        
//        
//        typedef Eigen::Matrix<double,dim,1> VectorDim;
//        typedef std::vector<VectorDim> PlaneNormalContainerType;
//        
//
//        
//    public:
//        
// //       static const PlaneNormalContainerType referencePlaneNormals;
//   
// //       static const PlaneNormalContainerType referencePlaneNormals;
//        
//        
//        static PlaneNormalContainerType getPlaneNormals(){
//            VectorDim  alpha(-1.0, 1.0,-1.0);
//			VectorDim  beta( 1.0,-1.0,-1.0);
//			VectorDim gamma(-1.0,-1.0, 1.0);
//			VectorDim delta( 1.0, 1.0, 1.0);
//            
//            PlaneNormalContainerType temp;
//            temp.push_back(alpha);
//            temp.push_back(beta);
//            temp.push_back(gamma);
//            temp.push_back(delta);
//            
//            return temp;
//        }
//        
//        
//    };
    
//    template <int dim>
//    const std::vector<Eigen::Matrix<double,dim,1> > FCC<dim>::referencePlaneNormals=FCC<dim>::getPlaneNormals();
//    const typename FCC<dim>::PlaneNormalContainerType FCC<dim>::referencePlaneNormals=FCC<dim>::getPlaneNormals();
    
    
    /**************************************************************************/
} // namespace model
#endif
