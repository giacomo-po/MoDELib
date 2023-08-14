/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_BCClattice_cpp_
#define model_BCClattice_cpp_

#include <BCClattice.h>
#include <DislocationMobilityBCC.h>

namespace model
{
    
        
        BCClattice<3>::BCClattice(const MatrixDim& Q,const PolycrystallineMaterialBase& material,const std::string& polyFile) :
        /* init */ SingleCrystalBase<dim>(getLatticeBasis(),Q)
        /* init */,PlaneNormalContainerType(getPlaneNormals())
        /* init */,SlipSystemContainerType(getSlipSystems(material,polyFile,*this))
        /* init */,SecondPhaseContainerType(getSecondPhases(material,*this))
        {
            
        }
        
        Eigen::Matrix<double,3,3> BCClattice<3>::getLatticeBasis()
        {/*!\returns The matrix of lattice vectors (cartesian cooridinates in columns),
          * in units of the crystallographic Burgers vector.
          */
            
            Eigen::Matrix<double,dim,dim> temp;
            temp << -1.0,  1.0,  1.0,
            /*   */  1.0, -1.0,  1.0,
            /*   */  1.0,  1.0, -1.0;
            
            return temp/sqrt(3.0);
        }

    const typename BCClattice<3>::PlaneNormalContainerType& BCClattice<3>::planeNormals() const
    {
        return *this;
    }

    const typename BCClattice<3>::SlipSystemContainerType& BCClattice<3>::slipSystems() const
    {
        return *this;
    }

    const typename BCClattice<3>::SecondPhaseContainerType& BCClattice<3>::secondPhases() const
    {
        return *this;
    }


//    const typename BCClattice<3>::DislocationMobilityContainerType& dislocationMobilities() const
//    {
//        return *this;
//    }


        
        std::vector<std::shared_ptr<LatticePlaneBase>> BCClattice<3>::getPlaneNormals() const
        {/*!\returns a std::vector of ReciprocalLatticeDirection(s) corresponding
          * the slip plane normals of the BCC lattice
          */
            
            typedef Eigen::Matrix<long int,dim,1> VectorDimI;
            typedef LatticeVector<dim> LatticeVectorType;

            LatticeVectorType a1((VectorDimI()<<1,0,0).finished(),*this);
            LatticeVectorType a2((VectorDimI()<<0,1,0).finished(),*this);
            LatticeVectorType a3((VectorDimI()<<0,0,1).finished(),*this);
            LatticeVectorType  y((VectorDimI()<<1,1,1).finished(),*this);

            std::vector<std::shared_ptr<LatticePlaneBase>> temp;
            temp.emplace_back(new LatticePlaneBase(a3,a1));
            temp.emplace_back(new LatticePlaneBase( y,a2));
            temp.emplace_back(new LatticePlaneBase(a2,a3));
            temp.emplace_back(new LatticePlaneBase( y,a1));
            temp.emplace_back(new LatticePlaneBase(a1,a2));
            temp.emplace_back(new LatticePlaneBase( y,a3));
            
            return temp;

        }
        
        std::vector<std::shared_ptr<SlipSystem>> BCClattice<3>::getSlipSystems(const PolycrystallineMaterialBase& material,
                                                                               const std::string& polyFile,
                                                                               const PlaneNormalContainerType& plN)
        {/*!\returns a std::vector of ReciprocalLatticeDirection(s) corresponding
          * the slip systems of the Hexagonal lattice
          */
            
            const std::string dislocationMobilityType(TextFileParser(polyFile).readString("dislocationMobilityType",true));
            DislocationMobilitySelector mobilitySelector("BCC");
            const std::shared_ptr<DislocationMobilityBase> mobility110(mobilitySelector.getMobility(dislocationMobilityType,material));

            
//            const std::shared_ptr<DislocationMobilityBase> mobility110(new DislocationMobilityBCC(material));

            
            const int solidSolutionNoiseMode(TextFileParser(polyFile).readScalar<int>("solidSolutionNoiseMode",true));
            const int stackingFaultNoiseMode(TextFileParser(polyFile).readScalar<int>("stackingFaultNoiseMode",true));
            std::shared_ptr<GlidePlaneNoise> planeNoise((solidSolutionNoiseMode||stackingFaultNoiseMode)? new GlidePlaneNoise(polyFile,material) : nullptr);

            
            typedef Eigen::Matrix<double,dim,1> VectorDimD;
            const double d110(this->reciprocalLatticeDirection(this->C2G*(VectorDimD()<<1.0,1.0,0.0).finished()).planeSpacing());

            std::vector<std::shared_ptr<SlipSystem>> temp;
            for(const auto& planeBase : plN)
            {
                if(std::fabs(planeBase->planeSpacing()-d110)<FLT_EPSILON)
                {// a {110} plane
                    const auto& a1(planeBase->primitiveVectors.first);
                    const auto& a3(planeBase->primitiveVectors.second);
                    temp.emplace_back(new SlipSystem(*planeBase, a1,mobility110,nullptr,planeNoise));
                    temp.emplace_back(new SlipSystem(*planeBase,a1*(-1),mobility110,nullptr,planeNoise));
                    temp.emplace_back(new SlipSystem(*planeBase, a3,mobility110,nullptr,planeNoise));
                    temp.emplace_back(new SlipSystem(*planeBase,a3*(-1),mobility110,nullptr,planeNoise));
                }
            }

            return temp;
        }
        
        

    std::vector<std::shared_ptr<SecondPhase<3>>> BCClattice<3>::getSecondPhases(const PolycrystallineMaterialBase& material,
                                                                                const PlaneNormalContainerType&)
    {
        
        const std::vector<std::string> spNames(TextFileParser(material.materialFile).readArray<std::string>("secondPhases",true));
        std::vector<std::shared_ptr<SecondPhase<3>>> temp;

        for(const std::string& sp : spNames)
        {
                throw std::runtime_error("Unnown SecondPhase "+sp+" in BCC crystal.");
        }
        return temp;
    }
        
} // namespace model
#endif
