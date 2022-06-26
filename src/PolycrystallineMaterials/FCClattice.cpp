/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_FCClattice_cpp_
#define model_FCClattice_cpp_

#include <FCClattice.h>

namespace model
{
    
        
        FCClattice<3>::FCClattice(const MatrixDim& Q) :
        /* init */ Lattice<dim>(getLatticeBasis(),Q)
        {
            
        }
        
        Eigen::Matrix<double,3,3> FCClattice<3>::getLatticeBasis()
        {/*!\returns The matrix of lattice vectors (cartesian cooridinates in columns),
          * in units of the crystallographic Burgers vector.
          */
            
            Eigen::Matrix<double,dim,dim> temp;
            temp << 0.0, 1.0, 1.0,
            /*   */ 1.0, 0.0, 1.0,
            /*   */ 1.0, 1.0, 0.0;
            
            return temp/sqrt(2.0);
        }
        
        std::vector<LatticePlaneBase> FCClattice<3>::reciprocalPlaneNormals(const Lattice<dim>& lat)
        {/*!\returns a std::vector of ReciprocalLatticeDirection(s) corresponding
          * the slip plane normals of the FCC lattice
          */
            
            
            typedef Eigen::Matrix<long int,dim,1> VectorDimI;
            
            typedef LatticeVector<dim> LatticeVectorType;
            LatticeVectorType a1(VectorDimI(1,0,0),lat); // [011]
            LatticeVectorType a2(VectorDimI(0,1,0),lat); // [101]
            LatticeVectorType a3(VectorDimI(0,0,1),lat); // [110]
            
            std::vector<LatticePlaneBase> temp;
            
            if(enable111planes)
            {// {111} planes
                temp.emplace_back(a1,a3);           // is (-1, 1,-1) in cartesian
                temp.emplace_back(a3,a2);           // is ( 1,-1,-1) in cartesian
                temp.emplace_back(a2,a1);           // is (-1,-1, 1) in cartesian
                temp.emplace_back(a1-a3,a2-a3);     // is ( 1, 1, 1) in cartesian
            }
            
            if(enable110planes)
            {// {110} planes
                temp.emplace_back(a1+a2-a3,a3);
                temp.emplace_back(a1+a2-a3,a1-a2);
                temp.emplace_back(a1+a3-a2,a2);
                temp.emplace_back(a1+a3-a2,a1-a3);
                temp.emplace_back(a2+a3-a1,a1);
                temp.emplace_back(a2+a3-a1,a2-a3);
            }
            
            
            return temp;
        }
        
        std::vector<std::shared_ptr<SlipSystem>> FCClattice<3>::slipSystems(const std::map<std::string,std::shared_ptr<DislocationMobilityBase>>& mobilities,
                                                                    const Lattice<dim>& lat,
                                                                    const PolycrystallineMaterialBase& material,
                                                                    const bool& enablePartials)
        {/*!\returns a std::vector of ReciprocalLatticeDirection(s) corresponding
          * the slip systems of the Hexagonal lattice
          */
            typedef Eigen::Matrix<long int,3,1> VectorDimI;
            
            
            typedef LatticeVector<3> LatticeVectorType;
            LatticeVectorType a1(VectorDimI(1,0,0),lat);
            LatticeVectorType a2(VectorDimI(0,1,0),lat);
            LatticeVectorType a3(VectorDimI(0,0,1),lat);
            
            //            std::shared_ptr<DislocationMobilityBase> fccMobility(new DislocationMobilityFCC(materialBase));
            const std::shared_ptr<DislocationMobilityBase>& fccMobility(mobilities.at("fcc"));
            
            std::vector<std::shared_ptr<SlipSystem>> temp;
            
            if(enable111planes)
            {// <110>{111}
                if(enablePartials)
                {
                    
                    const double ISF(TextFileParser(material.materialFile).readScalar<double>("ISF_SI",true)/(material.mu_SI*material.b_SI));
                    const double USF(TextFileParser(material.materialFile).readScalar<double>("USF_SI",true)/(material.mu_SI*material.b_SI));
                    const double MSF(TextFileParser(material.materialFile).readScalar<double>("MSF_SI",true)/(material.mu_SI*material.b_SI));
                    
                    const Eigen::Matrix<double,3,2> waveVectors((Eigen::Matrix<double,3,2>()<<0.0, 0.0,
                                                                 /*                        */ 0.0, 1.0,
                                                                 /*                        */ 1.0,-1.0
                                                                 ).finished());
                    
                    const Eigen::Matrix<double,4,3> f((Eigen::Matrix<double,4,3>()<<0.00,0.0, 0.0, 
                                                       /*                        */ 0.50,sqrt(3.0)/6.0, ISF,
                                                       /*                        */ 0.25,sqrt(3.0)/12.0,USF,
                                                       /*                        */ 1.00,sqrt(3.0)/3.0, MSF).finished());
                    
                    const int rotSymm(3);
                    const std::vector<Eigen::Matrix<double,2,1>> mirSymm;
                    std::shared_ptr<GammaSurface> gammaSurface0(new GammaSurface(LatticePlaneBase(a1,a3),waveVectors,f,rotSymm,mirSymm));
                    temp.emplace_back(new SlipSystem(a1,a3, RationalLatticeDirection<3>(Rational(1,3),(a1+a3)*(+1)),fccMobility,gammaSurface0));               // is (-1, 1,-1) in cartesian
                    temp.emplace_back(new SlipSystem(a1,a3, RationalLatticeDirection<3>(Rational(1,3),(a1+a3)*(-1)),fccMobility,gammaSurface0));               // is (-1, 1,-1) in cartesian
                    temp.emplace_back(new SlipSystem(a1,a3, RationalLatticeDirection<3>(Rational(1,3),(a1*2-a3)*(+1)),fccMobility,gammaSurface0));               // is (-1, 1,-1) in cartesian
                    temp.emplace_back(new SlipSystem(a1,a3, RationalLatticeDirection<3>(Rational(1,3),(a1*2-a3)*(-1)),fccMobility,gammaSurface0));               // is (-1, 1,-1) in cartesian
                    temp.emplace_back(new SlipSystem(a1,a3, RationalLatticeDirection<3>(Rational(1,3),(a3*2-a1)*(+1)),fccMobility,gammaSurface0));               // is (-1, 1,-1) in cartesian
                    temp.emplace_back(new SlipSystem(a1,a3, RationalLatticeDirection<3>(Rational(1,3),(a3*2-a1)*(-1)),fccMobility,gammaSurface0));               // is (-1, 1,-1) in cartesian
                    
                    std::shared_ptr<GammaSurface> gammaSurface1(new GammaSurface(LatticePlaneBase(a3,a2),waveVectors,f,rotSymm,mirSymm));
                    temp.emplace_back(new SlipSystem(a3,a2, RationalLatticeDirection<3>(Rational(1,3),(a3+a2)*(+1)),fccMobility,gammaSurface1));               // is (-1, 1,-1) in cartesian
                    temp.emplace_back(new SlipSystem(a3,a2, RationalLatticeDirection<3>(Rational(1,3),(a3+a2)*(-1)),fccMobility,gammaSurface1));               // is (-1, 1,-1) in cartesian
                    temp.emplace_back(new SlipSystem(a3,a2, RationalLatticeDirection<3>(Rational(1,3),(a3*2-a2)*(+1)),fccMobility,gammaSurface1));               // is (-1, 1,-1) in cartesian
                    temp.emplace_back(new SlipSystem(a3,a2, RationalLatticeDirection<3>(Rational(1,3),(a3*2-a2)*(-1)),fccMobility,gammaSurface1));               // is (-1, 1,-1) in cartesian
                    temp.emplace_back(new SlipSystem(a3,a2, RationalLatticeDirection<3>(Rational(1,3),(a2*2-a3)*(+1)),fccMobility,gammaSurface1));               // is (-1, 1,-1) in cartesian
                    temp.emplace_back(new SlipSystem(a3,a2, RationalLatticeDirection<3>(Rational(1,3),(a2*2-a3)*(-1)),fccMobility,gammaSurface1));               // is (-1, 1,-1) in cartesian
                    
                    std::shared_ptr<GammaSurface> gammaSurface2(new GammaSurface(LatticePlaneBase(a2,a1),waveVectors,f,rotSymm,mirSymm));
                    temp.emplace_back(new SlipSystem(a2,a1, RationalLatticeDirection<3>(Rational(1,3),(a2+a1)*(+1)),fccMobility,gammaSurface2));               // is (-1, 1,-1) in cartesian
                    temp.emplace_back(new SlipSystem(a2,a1, RationalLatticeDirection<3>(Rational(1,3),(a2+a1)*(-1)),fccMobility,gammaSurface2));               // is (-1, 1,-1) in cartesian
                    temp.emplace_back(new SlipSystem(a2,a1, RationalLatticeDirection<3>(Rational(1,3),(a2*2-a1)*(+1)),fccMobility,gammaSurface2));               // is (-1, 1,-1) in cartesian
                    temp.emplace_back(new SlipSystem(a2,a1, RationalLatticeDirection<3>(Rational(1,3),(a2*2-a1)*(-1)),fccMobility,gammaSurface2));               // is (-1, 1,-1) in cartesian
                    temp.emplace_back(new SlipSystem(a2,a1, RationalLatticeDirection<3>(Rational(1,3),(a1*2-a2)*(+1)),fccMobility,gammaSurface2));               // is (-1, 1,-1) in cartesian
                    temp.emplace_back(new SlipSystem(a2,a1, RationalLatticeDirection<3>(Rational(1,3),(a1*2-a2)*(-1)),fccMobility,gammaSurface2));               // is (-1, 1,-1) in cartesian
                    
                    std::shared_ptr<GammaSurface> gammaSurface3(new GammaSurface(LatticePlaneBase(a1-a3,a2-a3),waveVectors,f,rotSymm,mirSymm));
                    temp.emplace_back(new SlipSystem(a1-a3,a2-a3, RationalLatticeDirection<3>(Rational(1,3),(a1+a2-a3*2)*(+1)),fccMobility,gammaSurface3));               // is (-1, 1,-1) in cartesian
                    temp.emplace_back(new SlipSystem(a1-a3,a2-a3, RationalLatticeDirection<3>(Rational(1,3),(a1+a2-a3*2)*(-1)),fccMobility,gammaSurface3));               // is (-1, 1,-1) in cartesian
                    temp.emplace_back(new SlipSystem(a1-a3,a2-a3, RationalLatticeDirection<3>(Rational(1,3),(a1*2-a2-a3)*(+1)),fccMobility,gammaSurface3));               // is (-1, 1,-1) in cartesian
                    temp.emplace_back(new SlipSystem(a1-a3,a2-a3, RationalLatticeDirection<3>(Rational(1,3),(a1*2-a2-a3)*(-1)),fccMobility,gammaSurface3));               // is (-1, 1,-1) in cartesian
                    temp.emplace_back(new SlipSystem(a1-a3,a2-a3, RationalLatticeDirection<3>(Rational(1,3),(a2*2-a1-a3)*(+1)),fccMobility,gammaSurface3));               // is (-1, 1,-1) in cartesian
                    temp.emplace_back(new SlipSystem(a1-a3,a2-a3, RationalLatticeDirection<3>(Rational(1,3),(a2*2-a1-a3)*(-1)),fccMobility,gammaSurface3));                 // is (-1, 1,-1) in cartesian           // is (-1, 1,-1) in cartesian
                    
                }
                else
                {
                    temp.emplace_back(new SlipSystem(a1,a3, a1,fccMobility,nullptr));               // is (-1, 1,-1) in cartesian
                    temp.emplace_back(new SlipSystem(a1,a3,a1*(-1),fccMobility,nullptr));           // is (-1, 1,-1) in cartesian
                    temp.emplace_back(new SlipSystem(a1,a3, a3,fccMobility,nullptr));               // is (-1, 1,-1) in cartesian
                    temp.emplace_back(new SlipSystem(a1,a3,a3*(-1),fccMobility,nullptr));           // is (-1, 1,-1) in cartesian
                    temp.emplace_back(new SlipSystem(a1,a3,a1-a3,fccMobility,nullptr));             // is (-1, 1,-1) in cartesian
                    temp.emplace_back(new SlipSystem(a1,a3,a3-a1,fccMobility,nullptr));             // is (-1, 1,-1) in cartesian
                    
                    temp.emplace_back(new SlipSystem(a3,a2, a3,fccMobility,nullptr));               // is ( 1,-1,-1) in cartesian
                    temp.emplace_back(new SlipSystem(a3,a2,a3*(-1),fccMobility,nullptr));           // is ( 1,-1,-1) in cartesian
                    temp.emplace_back(new SlipSystem(a3,a2, a2,fccMobility,nullptr));               // is ( 1,-1,-1) in cartesian
                    temp.emplace_back(new SlipSystem(a3,a2,a2*(-1),fccMobility,nullptr));           // is ( 1,-1,-1) in cartesian
                    temp.emplace_back(new SlipSystem(a3,a2,a3-a2,fccMobility,nullptr));             // is ( 1,-1,-1) in cartesian
                    temp.emplace_back(new SlipSystem(a3,a2,a2-a3,fccMobility,nullptr));             // is ( 1,-1,-1) in cartesian
                    
                    temp.emplace_back(new SlipSystem(a2,a1, a2,fccMobility,nullptr));               // is (-1,-1, 1) in cartesian
                    temp.emplace_back(new SlipSystem(a2,a1,a2*(-1),fccMobility,nullptr));           // is (-1,-1, 1) in cartesian
                    temp.emplace_back(new SlipSystem(a2,a1, a1,fccMobility,nullptr));               // is (-1,-1, 1) in cartesian
                    temp.emplace_back(new SlipSystem(a2,a1,a1*(-1),fccMobility,nullptr));           // is (-1,-1, 1) in cartesian
                    temp.emplace_back(new SlipSystem(a2,a1,a2-a1,fccMobility,nullptr));             // is (-1,-1, 1) in cartesian
                    temp.emplace_back(new SlipSystem(a2,a1,a1-a2,fccMobility,nullptr));             // is (-1,-1, 1) in cartesian
                    
                    temp.emplace_back(new SlipSystem(a1-a3,a2-a3, a1-a3,fccMobility,nullptr));      // is ( 1, 1, 1) in cartesian
                    temp.emplace_back(new SlipSystem(a1-a3,a2-a3, a3-a1,fccMobility,nullptr));      // is ( 1, 1, 1) in cartesian
                    temp.emplace_back(new SlipSystem(a1-a3,a2-a3,a2-a3,fccMobility,nullptr));       // is ( 1, 1, 1) in cartesian
                    temp.emplace_back(new SlipSystem(a1-a3,a2-a3,a3-a2,fccMobility,nullptr));       // is ( 1, 1, 1) in cartesian
                    temp.emplace_back(new SlipSystem(a1-a3,a2-a3, a1-a2,fccMobility,nullptr));      // is ( 1, 1, 1) in cartesian
                    temp.emplace_back(new SlipSystem(a1-a3,a2-a3, a2-a1,fccMobility,nullptr));      // is ( 1, 1, 1) in cartesian
                }
            }
            
            
            
            if(enable110planes)
            {// <110>{110}
                temp.emplace_back(new SlipSystem(a1+a2-a3,a3, a3,fccMobility,nullptr));
                temp.emplace_back(new SlipSystem(a1+a2-a3,a3, a3*(-1),fccMobility,nullptr));
                
                temp.emplace_back(new SlipSystem(a1+a2-a3,a1-a2,a1-a2,fccMobility,nullptr));
                temp.emplace_back(new SlipSystem(a1+a2-a3,a1-a2,a2-a1,fccMobility,nullptr));
                
                temp.emplace_back(new SlipSystem(a1+a3-a2,a2, a2,fccMobility,nullptr));
                temp.emplace_back(new SlipSystem(a1+a3-a2,a2, a2*(-1),fccMobility,nullptr));
                
                temp.emplace_back(new SlipSystem(a1+a3-a2,a1-a3,a1-a3,fccMobility,nullptr));
                temp.emplace_back(new SlipSystem(a1+a3-a2,a1-a3,a3-a1,fccMobility,nullptr));
                
                temp.emplace_back(new SlipSystem(a2+a3-a1,a1, a1,fccMobility,nullptr));
                temp.emplace_back(new SlipSystem(a2+a3-a1,a1, a1*(-1),fccMobility,nullptr));
                
                temp.emplace_back(new SlipSystem(a2+a3-a1,a2-a3,a2-a3,fccMobility,nullptr));
                temp.emplace_back(new SlipSystem(a2+a3-a1,a2-a3,a3-a2,fccMobility,nullptr));
            }
            
            return temp;
        }
        
        
//template struct FCClattice<3>;

    
    
} // namespace model
#endif

