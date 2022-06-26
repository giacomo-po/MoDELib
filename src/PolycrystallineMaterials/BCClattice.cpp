/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_BCClattice_cpp_
#define model_BCClattice_cpp_

#include <vector>
#include <Eigen/Dense>

#include <LatticeModule.h>
#include <SlipSystem.h>
#include <PolycrystallineMaterialBase.h>
#include <BCClattice.h>

namespace model
{

        Eigen::Matrix<double,3,3> BCClattice<3>::getLatticeBasis()
        {/*!\returns The matrix of lattice vectors (in columns), in units of the
          * crystallographic Burgers vector.
          */
            
            Eigen::Matrix<double,dim,dim> temp;
            temp << -1.0,  1.0,  1.0,
            /*   */  1.0, -1.0,  1.0,
            /*   */  1.0,  1.0, -1.0;
            
            return temp/sqrt(3.0);
        }
        
        std::vector<LatticePlaneBase> BCClattice<3>::reciprocalPlaneNormals(const Lattice<dim>& lat)
        {/*!\returns a std::vector of ReciprocalLatticeDirection(s) corresponding
          * the slip plane normals of the FCC lattice
          */
            typedef Eigen::Matrix<long int,dim,1> VectorDimI;
            
            typedef LatticeVector<dim> LatticeVectorType;
            LatticeVectorType a1(VectorDimI(1,0,0),lat);
            LatticeVectorType a2(VectorDimI(0,1,0),lat);
            LatticeVectorType a3(VectorDimI(0,0,1),lat);
            LatticeVectorType  y(VectorDimI(1,1,1),lat);
            
            std::vector<LatticePlaneBase> temp;
            temp.emplace_back(a3,a1); // is ( 1, 0, 1) in cartesian
            temp.emplace_back( y,a2); // is ( 1, 0,-1) in cartesian
            temp.emplace_back(a2,a3); // is ( 0, 1, 1) in cartesian
            temp.emplace_back( y,a1); // is ( 0,-1, 1) in cartesian
            temp.emplace_back(a1,a2); // is ( 1, 1, 0) in cartesian
            temp.emplace_back( y,a3); // is (-1, 1, 0) in cartesian
            
            return temp;
        }

        std::vector<std::shared_ptr<SlipSystem>> BCClattice<3>::slipSystems(const std::map<std::string,std::shared_ptr<DislocationMobilityBase>>& mobilities,
                                                                    const Lattice<dim>& lat,
                                                                    const PolycrystallineMaterialBase& )
        {/*!\returns a std::vector of ReciprocalLatticeDirection(s) corresponding
          * the slip systems of the BCC lattice
          */
            typedef Eigen::Matrix<long int,3,1> VectorDimI;
            
            typedef LatticeVector<3> LatticeVectorType;
            LatticeVectorType a1(VectorDimI(1,0,0),lat);
            LatticeVectorType a2(VectorDimI(0,1,0),lat);
            LatticeVectorType a3(VectorDimI(0,0,1),lat);
            LatticeVectorType  y(VectorDimI(1,1,1),lat);
            
//            std::shared_ptr<DislocationMobilityBase> bccMobility(new DislocationMobilityBCC(materialBase));
            const std::shared_ptr<DislocationMobilityBase>& bccMobility(mobilities.at("bcc"));

            
            std::vector<std::shared_ptr<SlipSystem>> temp;
            
            temp.emplace_back(new SlipSystem(a3,a1, a3,bccMobility,nullptr)); // is ( 1, 0, 1) in cartesian
            temp.emplace_back(new SlipSystem(a3,a1,a3*(-1),bccMobility,nullptr)); // is ( 1, 0, 1) in cartesian
            temp.emplace_back(new SlipSystem(a3,a1, a1,bccMobility,nullptr)); // is ( 1, 0, 1) in cartesian
            temp.emplace_back(new SlipSystem(a3,a1,a1*(-1),bccMobility,nullptr)); // is ( 1, 0, 1) in cartesian
            
            temp.emplace_back(new SlipSystem( y,a2, y,bccMobility,nullptr)); // is ( 1, 0,-1) in cartesian
            temp.emplace_back(new SlipSystem( y,a2,y*(-1),bccMobility,nullptr)); // is ( 1, 0,-1) in cartesian
            temp.emplace_back(new SlipSystem( y,a2, a2,bccMobility,nullptr)); // is ( 1, 0,-1) in cartesian
            temp.emplace_back(new SlipSystem( y,a2,a2*(-1),bccMobility,nullptr)); // is ( 1, 0,-1) in cartesian
            
            temp.emplace_back(new SlipSystem(a2,a3, a2,bccMobility,nullptr)); // is ( 0, 1, 1) in cartesian
            temp.emplace_back(new SlipSystem(a2,a3,a2*(-1),bccMobility,nullptr)); // is ( 0, 1, 1) in cartesian
            temp.emplace_back(new SlipSystem(a2,a3, a3,bccMobility,nullptr)); // is ( 0, 1, 1) in cartesian
            temp.emplace_back(new SlipSystem(a2,a3,a3*(-1),bccMobility,nullptr)); // is ( 0, 1, 1) in cartesian
            
            temp.emplace_back(new SlipSystem( y,a1, y,bccMobility,nullptr)); // is ( 0,-1, 1) in cartesian
            temp.emplace_back(new SlipSystem( y,a1,y*(-1),bccMobility,nullptr)); // is ( 0,-1, 1) in cartesian
            temp.emplace_back(new SlipSystem( y,a1, a1,bccMobility,nullptr)); // is ( 0,-1, 1) in cartesian
            temp.emplace_back(new SlipSystem( y,a1,a1*(-1),bccMobility,nullptr)); // is ( 0,-1, 1) in cartesian
            
            temp.emplace_back(new SlipSystem(a1,a2, a1,bccMobility,nullptr)); // is ( 1, 1, 0) in cartesian
            temp.emplace_back(new SlipSystem(a1,a2,a1*(-1),bccMobility,nullptr)); // is ( 1, 1, 0) in cartesian
            temp.emplace_back(new SlipSystem(a1,a2, a2,bccMobility,nullptr)); // is ( 1, 1, 0) in cartesian
            temp.emplace_back(new SlipSystem(a1,a2,a2*(-1),bccMobility,nullptr)); // is ( 1, 1, 0) in cartesian
            
            temp.emplace_back(new SlipSystem( y,a3, y,bccMobility,nullptr)); // is (-1, 1, 0) in cartesian
            temp.emplace_back(new SlipSystem( y,a3,y*(-1),bccMobility,nullptr)); // is (-1, 1, 0) in cartesian
            temp.emplace_back(new SlipSystem( y,a3, a3,bccMobility,nullptr)); // is (-1, 1, 0) in cartesian
            temp.emplace_back(new SlipSystem( y,a3,a3*(-1),bccMobility,nullptr)); // is (-1, 1, 0) in cartesian
            
            return temp;
        }
  
//template struct BCClattice<3>;


} // namespace model
#endif

