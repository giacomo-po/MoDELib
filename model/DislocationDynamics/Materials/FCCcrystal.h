/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_FCCcrystal_H_
#define model_FCCcrystal_H_

#include <vector>
#include <Eigen/Dense>

#include <model/LatticeMath/LatticeMath.h>
#include <model/DislocationDynamics/Materials/SlipSystem.h>

namespace model
{
    
    struct FCC
    {
        static constexpr bool enable111planes=true;
        static constexpr bool enable110planes=false; // {110} planes don't work in this version since PlanePlane intersections can be off-lattice
        
        /**********************************************************************/
        template <int dim,typename MaterialType>
        static Eigen::Matrix<double,dim,dim> getLatticeBasis()
        {/*!\returns The matrix of lattice vectors (cartesian cooridinates in columns),
          * in units of the crystallographic Burgers vector.
          */
            
            Eigen::Matrix<double,dim,dim> temp;
            temp << 0.0, 1.0, 1.0,
            /*   */ 1.0, 0.0, 1.0,
            /*   */ 1.0, 1.0, 0.0;
            
            return temp/sqrt(2.0);
        }
        
        /**********************************************************************/
        template <int dim>
        static std::vector<LatticePlaneBase> reciprocalPlaneNormals()
        {/*!\returns a std::vector of ReciprocalLatticeDirection(s) corresponding
          * the slip plane normals of the FCC lattice
          */
            typedef Eigen::Matrix<long int,dim,1> VectorDimI;
            
            typedef LatticeVector<dim> LatticeVectorType;
            LatticeVectorType a1(VectorDimI(1,0,0));
            LatticeVectorType a2(VectorDimI(0,1,0));
            LatticeVectorType a3(VectorDimI(0,0,1));
            
            std::vector<LatticePlaneBase> temp;
            
            if(enable111planes)
            {
                temp.emplace_back(a1,a3);           // is (-1, 1,-1) in cartesian
                temp.emplace_back(a3,a2);           // is ( 1,-1,-1) in cartesian
                temp.emplace_back(a2,a1);           // is (-1,-1, 1) in cartesian
                temp.emplace_back(a1-a3,a2-a3);     // is ( 1, 1, 1) in cartesian
            }
            
            if(enable110planes)
            {
                temp.emplace_back(a1+a2-a3,a3);
                temp.emplace_back(a1+a2-a3,a1-a2);
                temp.emplace_back(a1+a3-a2,a2);
                temp.emplace_back(a1+a3-a2,a1-a3);
                temp.emplace_back(a2+a3-a1,a1);
                temp.emplace_back(a2+a3-a1,a2-a3);
                
            }
            
            return temp;
        }
        
        /**********************************************************************/
        static std::vector<SlipSystem> slipSystems()
        {/*!\returns a std::vector of ReciprocalLatticeDirection(s) corresponding
          * the slip plane normals of the FCC lattice
          */
            typedef Eigen::Matrix<long int,3,1> VectorDimI;
            
            typedef LatticeVector<3> LatticeVectorType;
            LatticeVectorType a1(VectorDimI(1,0,0));
            LatticeVectorType a2(VectorDimI(0,1,0));
            LatticeVectorType a3(VectorDimI(0,0,1));
            
            std::vector<SlipSystem> temp;
            
            if(enable111planes)
            {// <110>{111}
                
                temp.emplace_back(a1,a3, a1);           // is (-1, 1,-1) in cartesian
                temp.emplace_back(a1,a3,-a1);           // is (-1, 1,-1) in cartesian
                temp.emplace_back(a1,a3, a3);           // is (-1, 1,-1) in cartesian
                temp.emplace_back(a1,a3,-a3);           // is (-1, 1,-1) in cartesian
                temp.emplace_back(a1,a3,a1-a3);           // is (-1, 1,-1) in cartesian
                temp.emplace_back(a1,a3,a3-a1);           // is (-1, 1,-1) in cartesian
                
                temp.emplace_back(a3,a2, a3);           // is ( 1,-1,-1) in cartesian
                temp.emplace_back(a3,a2,-a3);           // is ( 1,-1,-1) in cartesian
                temp.emplace_back(a3,a2, a2);           // is ( 1,-1,-1) in cartesian
                temp.emplace_back(a3,a2,-a2);           // is ( 1,-1,-1) in cartesian
                temp.emplace_back(a3,a2,a3-a2);           // is ( 1,-1,-1) in cartesian
                temp.emplace_back(a3,a2,a2-a3);           // is ( 1,-1,-1) in cartesian
                
                temp.emplace_back(a2,a1, a2);           // is (-1,-1, 1) in cartesian
                temp.emplace_back(a2,a1,-a2);           // is (-1,-1, 1) in cartesian
                temp.emplace_back(a2,a1, a1);           // is (-1,-1, 1) in cartesian
                temp.emplace_back(a2,a1,-a1);           // is (-1,-1, 1) in cartesian
                temp.emplace_back(a2,a1,a2-a1);           // is (-1,-1, 1) in cartesian
                temp.emplace_back(a2,a1,a1-a2);           // is (-1,-1, 1) in cartesian
                
                temp.emplace_back(a1-a3,a2-a3, a1-a3);     // is ( 1, 1, 1) in cartesian
                temp.emplace_back(a1-a3,a2-a3, a3-a1);     // is ( 1, 1, 1) in cartesian
                temp.emplace_back(a1-a3,a2-a3,a2-a3);     // is ( 1, 1, 1) in cartesian
                temp.emplace_back(a1-a3,a2-a3,a3-a2);     // is ( 1, 1, 1) in cartesian
                temp.emplace_back(a1-a3,a2-a3, a1-a2);     // is ( 1, 1, 1) in cartesian
                temp.emplace_back(a1-a3,a2-a3, a2-a1);     // is ( 1, 1, 1) in cartesian
            }
            
            if(enable110planes)
            {// <110>{111}
                temp.emplace_back(a1+a2-a3,a3, a3);
                temp.emplace_back(a1+a2-a3,a3,-a3);
                
                temp.emplace_back(a1+a2-a3,a1-a2,a1-a2);
                temp.emplace_back(a1+a2-a3,a1-a2,a2-a1);
                
                temp.emplace_back(a1+a3-a2,a2, a2);
                temp.emplace_back(a1+a3-a2,a2,-a2);
                
                temp.emplace_back(a1+a3-a2,a1-a3,a1-a3);
                temp.emplace_back(a1+a3-a2,a1-a3,a3-a1);
                
                temp.emplace_back(a2+a3-a1,a1, a1);
                temp.emplace_back(a2+a3-a1,a1,-a1);
                
                temp.emplace_back(a2+a3-a1,a2-a3,a2-a3);
                temp.emplace_back(a2+a3-a1,a2-a3,a3-a2);
            }
            
            return temp;
        }
        
        
        
    };
    
    
} // namespace model
#endif

