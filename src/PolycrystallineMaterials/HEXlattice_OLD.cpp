/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_HEXlattice_cpp_
#define model_HEXlattice_cpp_

#include <vector>
#include <Eigen/Dense>

#include <LatticeModule.h>
#include <SlipSystem.h>
#include <PolycrystallineMaterialBase.h>
#include <DislocationMobilityHEXbasal.h>
#include <DislocationMobilityHEXprismatic.h>
#include <DislocationMobilityHEXpyramidal.h>
#include <HEXlattice.h>

namespace model
{
    
        
        HEXlattice<3>::HEXlattice(const MatrixDim& Q) :
        /* init */ Lattice<dim>(getLatticeBasis(),Q)
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
        
        std::vector<std::shared_ptr<LatticePlaneBase>> HEXlattice<3>::reciprocalPlaneNormals(const PolycrystallineMaterialBase& material,const Lattice<dim>& lat)
        {/*!\returns a std::vector of ReciprocalLatticeDirection(s) corresponding
          * the slip plane normals of the FCC lattice
          */
            
            const bool enableBasalSlipSystems(TextFileParser(material.materialFile).readScalar<int>("enableBasalSlipSystems",false));
            const bool enablePrismaticSlipSystems(TextFileParser(material.materialFile).readScalar<int>("enablePrismaticSlipSystems",false));
            const bool enablePyramidalSlipSystems(TextFileParser(material.materialFile).readScalar<int>("enablePyramidalSlipSystems",false));

            
            typedef Eigen::Matrix<long int,dim,1> VectorDimI;
            
            typedef LatticeVector<dim> LatticeVectorType;
            LatticeVectorType a1(VectorDimI(1,0,0),lat);
            LatticeVectorType a2(VectorDimI(0,1,0),lat);
            LatticeVectorType a3(a2-a1);
            LatticeVectorType c(VectorDimI(0,0,1),lat);
            
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
        
        std::vector<std::shared_ptr<SlipSystem>> HEXlattice<3>::slipSystems(const std::map<std::string,std::shared_ptr<DislocationMobilityBase>>& mobilities,
                                                                    const Lattice<dim>& lat,
                                                                    const PolycrystallineMaterialBase& material,
                                                                    const bool& enablePartials)
        {/*!\returns a std::vector of ReciprocalLatticeDirection(s) corresponding
          * the slip systems of the FCC lattice
          */
            
            const bool enableBasalSlipSystems(TextFileParser(material.materialFile).readScalar<int>("enableBasalSlipSystems",true));
            const bool enablePrismaticSlipSystems(TextFileParser(material.materialFile).readScalar<int>("enablePrismaticSlipSystems",true));
            const bool enablePyramidalSlipSystems(TextFileParser(material.materialFile).readScalar<int>("enablePyramidalSlipSystems",true));
            
            
            typedef Eigen::Matrix<long int,3,1> VectorDimI;
            
            typedef LatticeVector<3> LatticeVectorType;
            LatticeVectorType a1(VectorDimI(1,0,0),lat);
            LatticeVectorType a2(VectorDimI(0,1,0),lat);
            LatticeVectorType a3(a2-a1);
            LatticeVectorType c(VectorDimI(0,0,1),lat);
            
            const std::shared_ptr<DislocationMobilityBase>& basalMobility(mobilities.at("hexBasal"));
            const std::shared_ptr<DislocationMobilityBase>& prismaticMobility(mobilities.at("hexPrismatic"));
            const std::shared_ptr<DislocationMobilityBase>& pyramidalMobility(mobilities.at("hexPyramidal"));
            
            std::vector<std::shared_ptr<SlipSystem>> temp;
            
            if(enablePartials)
            {
                if(enableBasalSlipSystems)
                {
                    // Paritials for basal plane
                    const double ISF(TextFileParser(material.materialFile).readScalar<double>("ISF_SI",true)/(material.mu_SI*material.b_SI));
                    const double USF(TextFileParser(material.materialFile).readScalar<double>("USF_SI",true)/(material.mu_SI*material.b_SI));
                    const double MSF(TextFileParser(material.materialFile).readScalar<double>("MSF_SI",true)/(material.mu_SI*material.b_SI));
                    
                    const Eigen::Matrix<double,3,2> waveVectorsBasal((Eigen::Matrix<double,3,2>()<<0.0, 0.0, // value at origin
                                                                      /*                        */ 0.0, 1.0,
                                                                      /*                        */ 1.0,-1.0).finished());
                    
                    const Eigen::Matrix<double,4,3> fBasal((Eigen::Matrix<double,4,3>()<<0.00,0.0, 0.0, // value at origin
                                                            /*                        */ 0.50,sqrt(3.0)/6.0, ISF,
                                                            /*                        */ 0.25,sqrt(3.0)/12.0,USF,
                                                            /*                        */ 1.00,sqrt(3.0)/3.0, MSF).finished());
                    
                    const int rotSymmBasal(3);
                    const std::vector<Eigen::Matrix<double,2,1>> mirSymmBasal;
                    std::shared_ptr<GammaSurface> gammaSurface0(new GammaSurface(LatticePlaneBase(a1,a2),waveVectorsBasal,fBasal,rotSymmBasal,mirSymmBasal));
                    temp.emplace_back(new SlipSystem(a1,a2, RationalLatticeDirection<3>(Rational(1,3),(a1+a2)*(+1)),basalMobility,gammaSurface0));
                    temp.emplace_back(new SlipSystem(a1,a2, RationalLatticeDirection<3>(Rational(1,3),(a1+a2)*(-1)),basalMobility,gammaSurface0));
                    temp.emplace_back(new SlipSystem(a1,a2, RationalLatticeDirection<3>(Rational(1,3),(a3-a1)*(+1)),basalMobility,gammaSurface0));
                    temp.emplace_back(new SlipSystem(a1,a2, RationalLatticeDirection<3>(Rational(1,3),(a3-a1)*(-1)),basalMobility,gammaSurface0));
                    temp.emplace_back(new SlipSystem(a1,a2, RationalLatticeDirection<3>(Rational(1,3),(a2+a3)*(+1)),basalMobility,gammaSurface0));
                    temp.emplace_back(new SlipSystem(a1,a2, RationalLatticeDirection<3>(Rational(1,3),(a2+a3)*(-1)),basalMobility,gammaSurface0));
                }
                
                if(enablePrismaticSlipSystems)
                {
                    // Paritials for prism plane
                    const double PSF0(TextFileParser(material.materialFile).readScalar<double>("PSF0_SI",true)/(material.mu_SI*material.b_SI));
                    const double PSF1(TextFileParser(material.materialFile).readScalar<double>("PSF1_SI",true)/(material.mu_SI*material.b_SI));
                    const double PSF2(TextFileParser(material.materialFile).readScalar<double>("PSF2_SI",true)/(material.mu_SI*material.b_SI));
                    const double PSF3(TextFileParser(material.materialFile).readScalar<double>("PSF3_SI",true)/(material.mu_SI*material.b_SI));
                    
                    const Eigen::Matrix<double,5,2> waveVectorsPrism((Eigen::Matrix<double,5,2>()<<0.0, 0.0,
                                                                      /*                        */ 1.0, 0.0,
                                                                      /*                        */ 0.0, 1.0,
                                                                      /*                        */ 1.0, 1.0,
                                                                      /*                        */ 2.0, 0.0).finished());
                    
                    const Eigen::Matrix<double,5,3> fPrism((Eigen::Matrix<double,5,3>()<<0.00,0.0, 0.0, // value at origin
                                                            /*                        */ 0.50,              0.0,PSF0,
                                                            /*                        */ 0.00,sqrt(8.0/3.0)/2.0,PSF1,
                                                            /*                        */ 0.25,sqrt(8.0/3.0)/4.0,PSF2,
                                                            /*                        */ 0.50,sqrt(8.0/3.0)/2.0,PSF3).finished());
                    
                    const int rotSymmPrism(1);
                    std::vector<Eigen::Matrix<double,2,1>> mirSymmPrism;
                    mirSymmPrism.push_back((Eigen::Matrix<double,2,1>()<<1.0,0.0).finished()); // symm with respect to local y-axis
                    mirSymmPrism.push_back((Eigen::Matrix<double,2,1>()<<0.0,1.0).finished()); // symm with respect to local x-axis
                    
                    std::shared_ptr<GammaSurface> gammaSurface1(new GammaSurface(LatticePlaneBase(a1,c),waveVectorsPrism,fPrism,rotSymmPrism,mirSymmPrism));
                    temp.emplace_back(new SlipSystem(a1,c, RationalLatticeDirection<3>(Rational(1,2),a1*(+1)),prismaticMobility,gammaSurface1));
                    temp.emplace_back(new SlipSystem(a1,c, RationalLatticeDirection<3>(Rational(1,2),a1*(-1)),prismaticMobility,gammaSurface1));
                    
                    std::shared_ptr<GammaSurface> gammaSurface2(new GammaSurface(LatticePlaneBase(a2,c),waveVectorsPrism,fPrism,rotSymmPrism,mirSymmPrism));
                    temp.emplace_back(new SlipSystem(a2,c, RationalLatticeDirection<3>(Rational(1,2),a2*(+1)),prismaticMobility,gammaSurface2));
                    temp.emplace_back(new SlipSystem(a2,c, RationalLatticeDirection<3>(Rational(1,2),a2*(-1)),prismaticMobility,gammaSurface2));
                    
                    std::shared_ptr<GammaSurface> gammaSurface3(new GammaSurface(LatticePlaneBase(a3,c),waveVectorsPrism,fPrism,rotSymmPrism,mirSymmPrism));
                    temp.emplace_back(new SlipSystem(a3,c, RationalLatticeDirection<3>(Rational(1,2),a3*(+1)),prismaticMobility,gammaSurface3));
                    temp.emplace_back(new SlipSystem(a3,c, RationalLatticeDirection<3>(Rational(1,2),a3*(-1)),prismaticMobility,gammaSurface3));
                }
                
                
            }
            else
            {
                if(enableBasalSlipSystems)
                {
                    // <a> type slip
                    temp.emplace_back(new SlipSystem(a1,a2,a1,basalMobility,nullptr));           // basal plane
                    temp.emplace_back(new SlipSystem(a1,a2,a1*(-1),basalMobility,nullptr));           // basal plane
                    temp.emplace_back(new SlipSystem(a1,a2,a2,basalMobility,nullptr));           // basal plane
                    temp.emplace_back(new SlipSystem(a1,a2,a2*(-1),basalMobility,nullptr));           // basal plane
                    temp.emplace_back(new SlipSystem(a1,a2,a3,basalMobility,nullptr));           // basal plane
                    temp.emplace_back(new SlipSystem(a1,a2,a3*(-1),basalMobility,nullptr));           // basal plane
                }
                
                if(enablePrismaticSlipSystems)
                {
                    temp.emplace_back(new SlipSystem(a1,c,a1,prismaticMobility,nullptr));           // prismatic plane
                    temp.emplace_back(new SlipSystem(a1,c,a1*(-1),prismaticMobility,nullptr));           // prismatic plane
                    temp.emplace_back(new SlipSystem(a2,c,a2,prismaticMobility,nullptr));           // prismatic plane
                    temp.emplace_back(new SlipSystem(a2,c,a2*(-1),prismaticMobility,nullptr));           // prismatic plane
                    temp.emplace_back(new SlipSystem(a3,c,a3,prismaticMobility,nullptr));           // prismatic plane
                    temp.emplace_back(new SlipSystem(a3,c,a3*(-1),prismaticMobility,nullptr));           // prismatic plane
                }
                
                
                if(enablePyramidalSlipSystems)
                {
                    temp.emplace_back(new SlipSystem(a1,a2+c,a1,pyramidalMobility,nullptr));         // pyramidal plane
                    temp.emplace_back(new SlipSystem(a1,a2+c,a1*(-1),pyramidalMobility,nullptr));         // pyramidal plane
                    temp.emplace_back(new SlipSystem(a2,a3+c,a2,pyramidalMobility,nullptr));         // pyramidal plane
                    temp.emplace_back(new SlipSystem(a2,a3+c,a2*(-1),pyramidalMobility,nullptr));         // pyramidal plane
                    temp.emplace_back(new SlipSystem(a3,c-a1,a3,pyramidalMobility,nullptr));        // pyramidal plane
                    temp.emplace_back(new SlipSystem(a3,c-a1,a3*(-1),pyramidalMobility,nullptr));        // pyramidal plane
                    temp.emplace_back(new SlipSystem(a1*(-1),c-a2,a1,pyramidalMobility,nullptr));       // pyramidal plane
                    temp.emplace_back(new SlipSystem(a1*(-1),c-a2,a1*(-1),pyramidalMobility,nullptr));       // pyramidal plane
                    temp.emplace_back(new SlipSystem(a2*(-1),c-a3,a2,pyramidalMobility,nullptr));       // pyramidal plane
                    temp.emplace_back(new SlipSystem(a2*(-1),c-a3,a2*(-1),pyramidalMobility,nullptr));       // pyramidal plane
                    temp.emplace_back(new SlipSystem(a3*(-1),a1+c,a3,pyramidalMobility,nullptr));        // pyramidal plane
                    temp.emplace_back(new SlipSystem(a3*(-1),a1+c,a3*(-1),pyramidalMobility,nullptr));        // pyramidal plane
                }
                
                // <a+c> type slip
                // TO BE COMPLETED
                
            }
            
            return temp;
        }

//    template struct HEXlattice<3>;

}
#endif

