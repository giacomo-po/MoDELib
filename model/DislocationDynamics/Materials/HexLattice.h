/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_HexLattice_H_
#define model_HexLattice_H_

#include <vector>
#include <Eigen/Dense>

#include <model/LatticeMath/LatticeMath.h>
#include <model/DislocationDynamics/Materials/SlipSystem.h>

namespace model
{
    
    struct Hexagonal
    {

        /**********************************************************************/
        template <int dim,typename MaterialType>
        static Eigen::Matrix<double,dim,dim> getLatticeBasis()
        {/*!\returns The matrix of lattice vectors (cartesian cooridinates in columns),
          * in units of the crystallographic Burgers vector.
          */
            
            Eigen::Matrix<double,dim,dim> temp;
            temp << 1.0, 0.5,           0.0,
            /*   */ 0.0, 0.5*sqrt(3.0), 0.0,
            /*   */ 0.0, 0.0,           MaterialType::c2a;
            
            return temp;
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
            LatticeVectorType a3(a2-a1);
            LatticeVectorType c(VectorDimI(0,0,1));
            
            std::vector<LatticePlaneBase> temp;
            temp.emplace_back(a1,a2);           // basal plane
            
            temp.emplace_back(a1,c);           // prismatic plane
            temp.emplace_back(a2,c);           // prismatic plane
            temp.emplace_back(a3,c);           // prismatic plane
            
            temp.emplace_back(a1,a2+c);         // pyramidal plane
            temp.emplace_back(a2,a3+c);         // pyramidal plane
            temp.emplace_back(a3,-a1+c);        // pyramidal plane
            temp.emplace_back(-a1,-a2+c);       // pyramidal plane
            temp.emplace_back(-a2,-a3+c);       // pyramidal plane
            temp.emplace_back(-a3,a1+c);        // pyramidal plane
            
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
            LatticeVectorType a3(a2-a1);
            LatticeVectorType c(VectorDimI(0,0,1));
            
            std::vector<SlipSystem> temp;
            
            // <a> type slip
            temp.emplace_back(a1,a2, a1);           // basal plane
            temp.emplace_back(a1,a2,-a1);           // basal plane
            temp.emplace_back(a1,a2, a2);           // basal plane
            temp.emplace_back(a1,a2,-a2);           // basal plane
            temp.emplace_back(a1,a2, a3);           // basal plane
            temp.emplace_back(a1,a2,-a3);           // basal plane
            
            temp.emplace_back(a1,c, a1);           // prismatic plane
            temp.emplace_back(a1,c,-a1);           // prismatic plane
            temp.emplace_back(a2,c, a2);           // prismatic plane
            temp.emplace_back(a2,c,-a2);           // prismatic plane
            temp.emplace_back(a3,c, a3);           // prismatic plane
            temp.emplace_back(a3,c,-a3);           // prismatic plane
            
            temp.emplace_back(a1,a2+c, a1);         // pyramidal plane
            temp.emplace_back(a1,a2+c,-a1);         // pyramidal plane
            temp.emplace_back(a2,a3+c, a2);         // pyramidal plane
            temp.emplace_back(a2,a3+c,-a2);         // pyramidal plane
            temp.emplace_back(a3,-a1+c, a3);        // pyramidal plane
            temp.emplace_back(a3,-a1+c,-a3);        // pyramidal plane
            temp.emplace_back(-a1,-a2+c, a1);       // pyramidal plane
            temp.emplace_back(-a1,-a2+c,-a1);       // pyramidal plane
            temp.emplace_back(-a2,-a3+c, a2);       // pyramidal plane
            temp.emplace_back(-a2,-a3+c,-a2);       // pyramidal plane
            temp.emplace_back(-a3,a1+c, a3);        // pyramidal plane
            temp.emplace_back(-a3,a1+c,-a3);        // pyramidal plane
            
            // <a+c> type slip
            // TO BE COMPLETED
            
            return temp;
        }
        
        
        
    };
    
    
} // namespace model
#endif

