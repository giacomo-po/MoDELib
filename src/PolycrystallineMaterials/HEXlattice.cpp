/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_HEXlattice_cpp_
#define model_HEXlattice_cpp_

#include <HEXlattice.h>
#include <DislocationMobilityHEX.h>

namespace model
{
    
        
        HEXlattice<3>::HEXlattice(const MatrixDim& Q,const PolycrystallineMaterialBase& material,const std::string& polyFile) :
        /* init */ SingleCrystalBase<dim>(getLatticeBasis(),Q)
        /* init */,PlaneNormalContainerType(getPlaneNormals())
//        /* init */,DislocationMobilityContainerType(getMobilities())
        /* init */,SlipSystemContainerType(getSlipSystems(material,polyFile,*this))
        /* init */,SecondPhaseContainerType(getSecondPhases(material,*this))
        {
            
        }
        
        Eigen::Matrix<double,3,3> HEXlattice<3>::getLatticeBasis()
        {/*!\returns The matrix of lattice vectors (cartesian cooridinates in columns),
          * in units of the crystallographic Burgers vector.
          */
            
            Eigen::Matrix<double,dim,dim> temp;
            temp << 1.0, 0.5,           0.0,
            /*   */ 0.0, 0.5*sqrt(3.0), 0.0,
            /*   */ 0.0, 0.0,           sqrt(8.0/3.0);
            
            return temp;
        }

    const typename HEXlattice<3>::PlaneNormalContainerType& HEXlattice<3>::planeNormals() const
    {
        return *this;
    }

    const typename HEXlattice<3>::SlipSystemContainerType& HEXlattice<3>::slipSystems() const
    {
        return *this;
    }

    const typename HEXlattice<3>::SecondPhaseContainerType& HEXlattice<3>::secondPhases() const
    {
        return *this;
    }


//    const typename HEXlattice<3>::DislocationMobilityContainerType& dislocationMobilities() const
//    {
//        return *this;
//    }


        
        std::vector<std::shared_ptr<LatticePlaneBase>> HEXlattice<3>::getPlaneNormals() const
        {/*!\returns a std::vector of ReciprocalLatticeDirection(s) corresponding
          * the slip plane normals of the HEX lattice
          */
            
            const bool enableBasalSlipSystems(TextFileParser(material.materialFile).readScalar<int>("enableBasalSlipSystems",false));
            const bool enablePrismaticSlipSystems(TextFileParser(material.materialFile).readScalar<int>("enablePrismaticSlipSystems",false));
            const bool enablePyramidalSlipSystems(TextFileParser(material.materialFile).readScalar<int>("enablePyramidalSlipSystems",false));

            
            typedef Eigen::Matrix<long int,dim,1> VectorDimI;
            
            typedef LatticeVector<dim> LatticeVectorType;
            LatticeVectorType a1((VectorDimI()<<1,0,0).finished(),*this);
            LatticeVectorType a2((VectorDimI()<<0,1,0).finished(),*this);
            LatticeVectorType a3(a2-a1);
            LatticeVectorType  c((VectorDimI()<<0,0,1).finished(),*this);

            std::vector<std::shared_ptr<LatticePlaneBase>> temp;
            if(enableBasalSlipSystems)
            {
                temp.emplace_back(new LatticePlaneBase(a1,a2));           // basal plane
            }

            if(enablePrismaticSlipSystems)
            {
                temp.emplace_back(new LatticePlaneBase(a1,c));           // prismatic plane
                temp.emplace_back(new LatticePlaneBase(a2,c));           // prismatic plane
                temp.emplace_back(new LatticePlaneBase(a3,c));           // prismatic plane
            }
            
            if(enablePyramidalSlipSystems)
            {
                temp.emplace_back(new LatticePlaneBase(a1,a2+c));         // pyramidal plane
                temp.emplace_back(new LatticePlaneBase(a2,a3+c));         // pyramidal plane
                temp.emplace_back(new LatticePlaneBase(a3,c-a1));        // pyramidal plane
                temp.emplace_back(new LatticePlaneBase(a1*(-1),c-a2));       // pyramidal plane
                temp.emplace_back(new LatticePlaneBase(a2*(-1),c-a3));       // pyramidal plane
                temp.emplace_back(new LatticePlaneBase(a3*(-1),a1+c));        // pyramidal plane
            }
            
            return temp;
        }
        
        std::vector<std::shared_ptr<SlipSystem>> HEXlattice<3>::getSlipSystems(const PolycrystallineMaterialBase& material,
                                                                               const std::string& polyFile,
                                                                               const PlaneNormalContainerType& plN)
        {/*!\returns a std::vector of ReciprocalLatticeDirection(s) corresponding
          * the slip systems of the Hexagonal lattice
          */
            
            
            
            std::vector<std::shared_ptr<SlipSystem>> temp;
            
                FINISH THIS
            
            return temp;
        }
        
        

    std::vector<std::shared_ptr<SecondPhase<3>>> HEXlattice<3>::getSecondPhases(const PolycrystallineMaterialBase& material,
                                                                                const PlaneNormalContainerType& plN)
    {
        
        const std::vector<std::string> spNames(TextFileParser(material.materialFile).readArray<std::string>("secondPhases",true));
        std::vector<std::shared_ptr<SecondPhase<3>>> temp;

        for(const std::string& sp : spNames)
        {
                throw std::runtime_error("Unnown SecondPhase "+sp+" in HEX crystal.");
        }
        return temp;
    }
        
} // namespace model
#endif

