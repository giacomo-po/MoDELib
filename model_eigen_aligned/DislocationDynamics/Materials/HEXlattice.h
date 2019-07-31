/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_HEXlattice_H_
#define model_HEXlattice_H_

#include <vector>
#include <Eigen/Dense>

#include <LatticeMath.h>
#include <SlipSystem.h>
#include <DislocatedMaterial.h>


namespace model
{
    
    template<int dim>
    struct HEXlattice
    {
        
    };
    
    template<>
    struct HEXlattice<3> : public Lattice<3>
    {

        static constexpr int dim=3;
        static constexpr auto name="HEX";
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;

        /**********************************************************************/
        HEXlattice(const MatrixDim& Q) :
        /* init */ Lattice<dim>(getLatticeBasis(),Q)
        {
            
        }
        
        /**********************************************************************/
        static Eigen::Matrix<double,dim,dim> getLatticeBasis()
        {/*!\returns The matrix of lattice vectors (cartesian cooridinates in columns),
          * in units of the crystallographic Burgers vector.
          */
            
            Eigen::Matrix<double,dim,dim> temp;
            temp << 1.0, 0.5,           0.0,
            /*   */ 0.0, 0.5*sqrt(3.0), 0.0,
            /*   */ 0.0, 0.0,           sqrt(8.0/3.0);
            
            return temp;
        }
        
        /**********************************************************************/
        static std::vector<LatticePlaneBase> reciprocalPlaneNormals(const Lattice<dim>& lat)
        {/*!\returns a std::vector of ReciprocalLatticeDirection(s) corresponding
          * the slip plane normals of the FCC lattice
          */
            typedef Eigen::Matrix<long int,dim,1> VectorDimI;
            
            typedef LatticeVector<dim> LatticeVectorType;
            LatticeVectorType a1(VectorDimI(1,0,0),lat);
            LatticeVectorType a2(VectorDimI(0,1,0),lat);
            LatticeVectorType a3(a2-a1);
            LatticeVectorType c(VectorDimI(0,0,1),lat);
            
            std::vector<LatticePlaneBase> temp;
            temp.emplace_back(a1,a2);           // basal plane
            
            temp.emplace_back(a1,c);           // prismatic plane
            temp.emplace_back(a2,c);           // prismatic plane
            temp.emplace_back(a3,c);           // prismatic plane
            
            temp.emplace_back(a1,a2+c);         // pyramidal plane
            temp.emplace_back(a2,a3+c);         // pyramidal plane
            temp.emplace_back(a3,c-a1);        // pyramidal plane
            temp.emplace_back(a1*(-1),c-a2);       // pyramidal plane
            temp.emplace_back(a2*(-1),c-a3);       // pyramidal plane
            temp.emplace_back(a3*(-1),a1+c);        // pyramidal plane
            
            return temp;
        }
        
        /**********************************************************************/
        static std::vector<std::shared_ptr<SlipSystem>> slipSystems(const DislocatedMaterial<dim,Isotropic>& material,const Lattice<dim>& lat)
        {/*!\returns a std::vector of ReciprocalLatticeDirection(s) corresponding
          * the slip plane normals of the FCC lattice
          */
            typedef Eigen::Matrix<long int,3,1> VectorDimI;
            
            typedef LatticeVector<3> LatticeVectorType;
            LatticeVectorType a1(VectorDimI(1,0,0),lat);
            LatticeVectorType a2(VectorDimI(0,1,0),lat);
            LatticeVectorType a3(a2-a1);
            LatticeVectorType c(VectorDimI(0,0,1),lat);
            

            std::vector<std::shared_ptr<SlipSystem>> temp;
            
            // <a> type slip
            temp.emplace_back(new SlipSystem(a1,a2,a1,material.dislocationMobilities.at("hexBasal")));           // basal plane
            temp.emplace_back(new SlipSystem(a1,a2,a1*(-1),material.dislocationMobilities.at("hexBasal")));           // basal plane
            temp.emplace_back(new SlipSystem(a1,a2,a2,material.dislocationMobilities.at("hexBasal")));           // basal plane
            temp.emplace_back(new SlipSystem(a1,a2,a2*(-1),material.dislocationMobilities.at("hexBasal")));           // basal plane
            temp.emplace_back(new SlipSystem(a1,a2,a3,material.dislocationMobilities.at("hexBasal")));           // basal plane
            temp.emplace_back(new SlipSystem(a1,a2,a3*(-1),material.dislocationMobilities.at("hexBasal")));           // basal plane
            
            temp.emplace_back(new SlipSystem(a1,c,a1,material.dislocationMobilities.at("hexPrismatic")));           // prismatic plane
            temp.emplace_back(new SlipSystem(a1,c,a1*(-1),material.dislocationMobilities.at("hexPrismatic")));           // prismatic plane
            temp.emplace_back(new SlipSystem(a2,c,a2,material.dislocationMobilities.at("hexPrismatic")));           // prismatic plane
            temp.emplace_back(new SlipSystem(a2,c,a2*(-1),material.dislocationMobilities.at("hexPrismatic")));           // prismatic plane
            temp.emplace_back(new SlipSystem(a3,c,a3,material.dislocationMobilities.at("hexPrismatic")));           // prismatic plane
            temp.emplace_back(new SlipSystem(a3,c,a3*(-1),material.dislocationMobilities.at("hexPrismatic")));           // prismatic plane
            
            temp.emplace_back(new SlipSystem(a1,a2+c,a1,material.dislocationMobilities.at("hexPyramidal")));         // pyramidal plane
            temp.emplace_back(new SlipSystem(a1,a2+c,a1*(-1),material.dislocationMobilities.at("hexPyramidal")));         // pyramidal plane
            temp.emplace_back(new SlipSystem(a2,a3+c,a2,material.dislocationMobilities.at("hexPyramidal")));         // pyramidal plane
            temp.emplace_back(new SlipSystem(a2,a3+c,a2*(-1),material.dislocationMobilities.at("hexPyramidal")));         // pyramidal plane
            temp.emplace_back(new SlipSystem(a3,c-a1,a3,material.dislocationMobilities.at("hexPyramidal")));        // pyramidal plane
            temp.emplace_back(new SlipSystem(a3,c-a1,a3*(-1),material.dislocationMobilities.at("hexPyramidal")));        // pyramidal plane
            temp.emplace_back(new SlipSystem(a1*(-1),c-a2,a1,material.dislocationMobilities.at("hexPyramidal")));       // pyramidal plane
            temp.emplace_back(new SlipSystem(a1*(-1),c-a2,a1*(-1),material.dislocationMobilities.at("hexPyramidal")));       // pyramidal plane
            temp.emplace_back(new SlipSystem(a2*(-1),c-a3,a2,material.dislocationMobilities.at("hexPyramidal")));       // pyramidal plane
            temp.emplace_back(new SlipSystem(a2*(-1),c-a3,a2*(-1),material.dislocationMobilities.at("hexPyramidal")));       // pyramidal plane
            temp.emplace_back(new SlipSystem(a3*(-1),a1+c,a3,material.dislocationMobilities.at("hexPyramidal")));        // pyramidal plane
            temp.emplace_back(new SlipSystem(a3*(-1),a1+c,a3*(-1),material.dislocationMobilities.at("hexPyramidal")));        // pyramidal plane
            
            // <a+c> type slip
            // TO BE COMPLETED
            
            return temp;
        }
        
        
        
    };
    
    
} // namespace model
#endif

