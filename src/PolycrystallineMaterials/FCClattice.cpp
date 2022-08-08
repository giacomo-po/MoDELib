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
#include <DislocationMobilityFCC.h>

namespace model
{
    
        
        FCClattice<3>::FCClattice(const MatrixDim& Q,const PolycrystallineMaterialBase& material,const std::string& polyFile) :
        /* init */ SingleCrystalBase<dim>(getLatticeBasis(),Q)
        /* init */,PlaneNormalContainerType(getPlaneNormals())
//        /* init */,DislocationMobilityContainerType(getMobilities())
        /* init */,SlipSystemContainerType(getSlipSystems(material,polyFile,*this))
        /* init */,SecondPhaseContainerType(getSecondPhases(material,*this))
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

    const typename FCClattice<3>::PlaneNormalContainerType& FCClattice<3>::planeNormals() const
    {
        return *this;
    }

    const typename FCClattice<3>::SlipSystemContainerType& FCClattice<3>::slipSystems() const
    {
        return *this;
    }

    const typename FCClattice<3>::SecondPhaseContainerType& FCClattice<3>::secondPhases() const
    {
        return *this;
    }


//    const typename FCClattice<3>::DislocationMobilityContainerType& dislocationMobilities() const
//    {
//        return *this;
//    }


        
        std::vector<std::shared_ptr<LatticePlaneBase>> FCClattice<3>::getPlaneNormals() const
        {/*!\returns a std::vector of ReciprocalLatticeDirection(s) corresponding
          * the slip plane normals of the FCC lattice
          */
            
            
            typedef Eigen::Matrix<long int,dim,1> VectorDimI;

            typedef LatticeVector<dim> LatticeVectorType;
            LatticeVectorType a1((VectorDimI()<<1,0,0).finished(),*this);
            LatticeVectorType a2((VectorDimI()<<0,1,0).finished(),*this);
            LatticeVectorType a3((VectorDimI()<<0,0,1).finished(),*this);
            
            std::vector<std::shared_ptr<LatticePlaneBase>> temp;
            
            if(enable111planes)
            {// {111} planes
                temp.emplace_back(new LatticePlaneBase(a1,a3));           // is (-1, 1,-1) in cartesian
                temp.emplace_back(new LatticePlaneBase(a3,a2));           // is ( 1,-1,-1) in cartesian
                temp.emplace_back(new LatticePlaneBase(a2,a1));           // is (-1,-1, 1) in cartesian
                temp.emplace_back(new LatticePlaneBase(a1-a3,a2-a3));     // is ( 1, 1, 1) in cartesian
            }
            
//            if(enable110planes)
//            {// {110} planes
//                temp.emplace_back(new LatticePlaneBase(a1+a2-a3,a3));
//                temp.emplace_back(new LatticePlaneBase(a1+a2-a3,a1-a2));
//                temp.emplace_back(new LatticePlaneBase(a1+a3-a2,a2));
//                temp.emplace_back(new LatticePlaneBase(a1+a3-a2,a1-a3));
//                temp.emplace_back(new LatticePlaneBase(a2+a3-a1,a1));
//                temp.emplace_back(new LatticePlaneBase(a2+a3-a1,a2-a3));
//            }
            
            
            
            
            return temp;
        }
        
        std::vector<std::shared_ptr<SlipSystem>> FCClattice<3>::getSlipSystems(const PolycrystallineMaterialBase& material,
                                                                               const std::string& polyFile,
                                                                               const PlaneNormalContainerType& plN)
        {/*!\returns a std::vector of ReciprocalLatticeDirection(s) corresponding
          * the slip systems of the Hexagonal lattice
          */
            const std::shared_ptr<DislocationMobilityBase> fccMobility(new DislocationMobilityFCC(material));
            
            std::vector<std::shared_ptr<SlipSystem>> temp;
            
            const bool enablePartials(TextFileParser(polyFile).readScalar<int>("enablePartials",true));
            
//            typedef Eigen::Matrix<long int,dim,1> VectorDimI;
            typedef Eigen::Matrix<double,dim,1> VectorDimD;

            
//            const double d111(ReciprocalLatticeVector<3>(VectorDimI(1,1,1), *this).planeSpacing());
            const double d111(this->reciprocalLatticeDirection((VectorDimD()<<1.0,1.0,1.0).finished()).planeSpacing());

//            for(const auto& planeBase : plN)
//            {
//                const auto& a1(planeBase->primitiveVectors.first);
//                const auto& a3(planeBase->primitiveVectors.second);
//                if(std::fabs(planeBase->planeSpacing()-d111)<FLT_EPSILON)
//                {// a {111} plane
//
//                }
//            }
            
            
//            if(enable111planes)
//            {// <110>{111}
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
                    for(const auto& planeBase : plN)
                    {
                        if(std::fabs(planeBase->planeSpacing()-d111)<FLT_EPSILON)
                        {// a {111} plane
                            std::shared_ptr<GammaSurface> gammaSurface(new GammaSurface(*planeBase,waveVectors,f,rotSymm,mirSymm));
                            const auto& a1(planeBase->primitiveVectors.first);
                            const auto& a3(planeBase->primitiveVectors.second);
                            temp.emplace_back(new SlipSystem(*planeBase, RationalLatticeDirection<3>(Rational(1,3),(a1+a3)*(+1)),fccMobility,gammaSurface));               // is (-1, 1,-1) in cartesian
                            temp.emplace_back(new SlipSystem(*planeBase, RationalLatticeDirection<3>(Rational(1,3),(a1+a3)*(-1)),fccMobility,gammaSurface));               // is (-1, 1,-1) in cartesian
                            temp.emplace_back(new SlipSystem(*planeBase, RationalLatticeDirection<3>(Rational(1,3),(a1*2-a3)*(+1)),fccMobility,gammaSurface));               // is (-1, 1,-1) in cartesian
                            temp.emplace_back(new SlipSystem(*planeBase, RationalLatticeDirection<3>(Rational(1,3),(a1*2-a3)*(-1)),fccMobility,gammaSurface));               // is (-1, 1,-1) in cartesian
                            temp.emplace_back(new SlipSystem(*planeBase, RationalLatticeDirection<3>(Rational(1,3),(a3*2-a1)*(+1)),fccMobility,gammaSurface));               // is (-1, 1,-1) in cartesian
                            temp.emplace_back(new SlipSystem(*planeBase, RationalLatticeDirection<3>(Rational(1,3),(a3*2-a1)*(-1)),fccMobility,gammaSurface));               // is (-1, 1,-1) in cartesian
                        }
                    }
                    
                    
//                    std::shared_ptr<GammaSurface> gammaSurface0(new GammaSurface(LatticePlaneBase(a1,a3),waveVectors,f,rotSymm,mirSymm));
//                    temp.emplace_back(new SlipSystem(a1,a3, RationalLatticeDirection<3>(Rational(1,3),(a1+a3)*(+1)),fccMobility,gammaSurface0));               // is (-1, 1,-1) in cartesian
//                    temp.emplace_back(new SlipSystem(a1,a3, RationalLatticeDirection<3>(Rational(1,3),(a1+a3)*(-1)),fccMobility,gammaSurface0));               // is (-1, 1,-1) in cartesian
//                    temp.emplace_back(new SlipSystem(a1,a3, RationalLatticeDirection<3>(Rational(1,3),(a1*2-a3)*(+1)),fccMobility,gammaSurface0));               // is (-1, 1,-1) in cartesian
//                    temp.emplace_back(new SlipSystem(a1,a3, RationalLatticeDirection<3>(Rational(1,3),(a1*2-a3)*(-1)),fccMobility,gammaSurface0));               // is (-1, 1,-1) in cartesian
//                    temp.emplace_back(new SlipSystem(a1,a3, RationalLatticeDirection<3>(Rational(1,3),(a3*2-a1)*(+1)),fccMobility,gammaSurface0));               // is (-1, 1,-1) in cartesian
//                    temp.emplace_back(new SlipSystem(a1,a3, RationalLatticeDirection<3>(Rational(1,3),(a3*2-a1)*(-1)),fccMobility,gammaSurface0));               // is (-1, 1,-1) in cartesian
//
//                    std::shared_ptr<GammaSurface> gammaSurface1(new GammaSurface(LatticePlaneBase(a3,a2),waveVectors,f,rotSymm,mirSymm));
//                    temp.emplace_back(new SlipSystem(a3,a2, RationalLatticeDirection<3>(Rational(1,3),(a3+a2)*(+1)),fccMobility,gammaSurface1));               // is (-1, 1,-1) in cartesian
//                    temp.emplace_back(new SlipSystem(a3,a2, RationalLatticeDirection<3>(Rational(1,3),(a3+a2)*(-1)),fccMobility,gammaSurface1));               // is (-1, 1,-1) in cartesian
//                    temp.emplace_back(new SlipSystem(a3,a2, RationalLatticeDirection<3>(Rational(1,3),(a3*2-a2)*(+1)),fccMobility,gammaSurface1));               // is (-1, 1,-1) in cartesian
//                    temp.emplace_back(new SlipSystem(a3,a2, RationalLatticeDirection<3>(Rational(1,3),(a3*2-a2)*(-1)),fccMobility,gammaSurface1));               // is (-1, 1,-1) in cartesian
//                    temp.emplace_back(new SlipSystem(a3,a2, RationalLatticeDirection<3>(Rational(1,3),(a2*2-a3)*(+1)),fccMobility,gammaSurface1));               // is (-1, 1,-1) in cartesian
//                    temp.emplace_back(new SlipSystem(a3,a2, RationalLatticeDirection<3>(Rational(1,3),(a2*2-a3)*(-1)),fccMobility,gammaSurface1));               // is (-1, 1,-1) in cartesian
//
//                    std::shared_ptr<GammaSurface> gammaSurface2(new GammaSurface(LatticePlaneBase(a2,a1),waveVectors,f,rotSymm,mirSymm));
//                    temp.emplace_back(new SlipSystem(a2,a1, RationalLatticeDirection<3>(Rational(1,3),(a2+a1)*(+1)),fccMobility,gammaSurface2));               // is (-1, 1,-1) in cartesian
//                    temp.emplace_back(new SlipSystem(a2,a1, RationalLatticeDirection<3>(Rational(1,3),(a2+a1)*(-1)),fccMobility,gammaSurface2));               // is (-1, 1,-1) in cartesian
//                    temp.emplace_back(new SlipSystem(a2,a1, RationalLatticeDirection<3>(Rational(1,3),(a2*2-a1)*(+1)),fccMobility,gammaSurface2));               // is (-1, 1,-1) in cartesian
//                    temp.emplace_back(new SlipSystem(a2,a1, RationalLatticeDirection<3>(Rational(1,3),(a2*2-a1)*(-1)),fccMobility,gammaSurface2));               // is (-1, 1,-1) in cartesian
//                    temp.emplace_back(new SlipSystem(a2,a1, RationalLatticeDirection<3>(Rational(1,3),(a1*2-a2)*(+1)),fccMobility,gammaSurface2));               // is (-1, 1,-1) in cartesian
//                    temp.emplace_back(new SlipSystem(a2,a1, RationalLatticeDirection<3>(Rational(1,3),(a1*2-a2)*(-1)),fccMobility,gammaSurface2));               // is (-1, 1,-1) in cartesian
//
//                    std::shared_ptr<GammaSurface> gammaSurface3(new GammaSurface(LatticePlaneBase(a1-a3,a2-a3),waveVectors,f,rotSymm,mirSymm));
//                    temp.emplace_back(new SlipSystem(a1-a3,a2-a3, RationalLatticeDirection<3>(Rational(1,3),(a1+a2-a3*2)*(+1)),fccMobility,gammaSurface3));               // is (-1, 1,-1) in cartesian
//                    temp.emplace_back(new SlipSystem(a1-a3,a2-a3, RationalLatticeDirection<3>(Rational(1,3),(a1+a2-a3*2)*(-1)),fccMobility,gammaSurface3));               // is (-1, 1,-1) in cartesian
//                    temp.emplace_back(new SlipSystem(a1-a3,a2-a3, RationalLatticeDirection<3>(Rational(1,3),(a1*2-a2-a3)*(+1)),fccMobility,gammaSurface3));               // is (-1, 1,-1) in cartesian
//                    temp.emplace_back(new SlipSystem(a1-a3,a2-a3, RationalLatticeDirection<3>(Rational(1,3),(a1*2-a2-a3)*(-1)),fccMobility,gammaSurface3));               // is (-1, 1,-1) in cartesian
//                    temp.emplace_back(new SlipSystem(a1-a3,a2-a3, RationalLatticeDirection<3>(Rational(1,3),(a2*2-a1-a3)*(+1)),fccMobility,gammaSurface3));               // is (-1, 1,-1) in cartesian
//                    temp.emplace_back(new SlipSystem(a1-a3,a2-a3, RationalLatticeDirection<3>(Rational(1,3),(a2*2-a1-a3)*(-1)),fccMobility,gammaSurface3));                 // is (-1, 1,-1) in cartesian           // is (-1, 1,-1) in cartesian
                    
                }
                else
                {
                    
                    for(const auto& planeBase : plN)
                    {
                        if(std::fabs(planeBase->planeSpacing()-d111)<FLT_EPSILON)
                        {// a {111} plane
                        const auto& a1(planeBase->primitiveVectors.first);
                        const auto& a3(planeBase->primitiveVectors.second);
                        temp.emplace_back(new SlipSystem(*planeBase, a1,fccMobility,nullptr));               // is (-1, 1,-1) in cartesian
                        temp.emplace_back(new SlipSystem(*planeBase,a1*(-1),fccMobility,nullptr));           // is (-1, 1,-1) in cartesian
                        temp.emplace_back(new SlipSystem(*planeBase, a3,fccMobility,nullptr));               // is (-1, 1,-1) in cartesian
                        temp.emplace_back(new SlipSystem(*planeBase,a3*(-1),fccMobility,nullptr));           // is (-1, 1,-1) in cartesian
                        temp.emplace_back(new SlipSystem(*planeBase,a1-a3,fccMobility,nullptr));             // is (-1, 1,-1) in cartesian
                        temp.emplace_back(new SlipSystem(*planeBase,a3-a1,fccMobility,nullptr));             // is (-1, 1,-1) in cartesian
                        }
                    }
                    
//                    temp.emplace_back(new SlipSystem(a1,a3, a1,fccMobility,nullptr));               // is (-1, 1,-1) in cartesian
//                    temp.emplace_back(new SlipSystem(a1,a3,a1*(-1),fccMobility,nullptr));           // is (-1, 1,-1) in cartesian
//                    temp.emplace_back(new SlipSystem(a1,a3, a3,fccMobility,nullptr));               // is (-1, 1,-1) in cartesian
//                    temp.emplace_back(new SlipSystem(a1,a3,a3*(-1),fccMobility,nullptr));           // is (-1, 1,-1) in cartesian
//                    temp.emplace_back(new SlipSystem(a1,a3,a1-a3,fccMobility,nullptr));             // is (-1, 1,-1) in cartesian
//                    temp.emplace_back(new SlipSystem(a1,a3,a3-a1,fccMobility,nullptr));             // is (-1, 1,-1) in cartesian
//
//                    temp.emplace_back(new SlipSystem(a3,a2, a3,fccMobility,nullptr));               // is ( 1,-1,-1) in cartesian
//                    temp.emplace_back(new SlipSystem(a3,a2,a3*(-1),fccMobility,nullptr));           // is ( 1,-1,-1) in cartesian
//                    temp.emplace_back(new SlipSystem(a3,a2, a2,fccMobility,nullptr));               // is ( 1,-1,-1) in cartesian
//                    temp.emplace_back(new SlipSystem(a3,a2,a2*(-1),fccMobility,nullptr));           // is ( 1,-1,-1) in cartesian
//                    temp.emplace_back(new SlipSystem(a3,a2,a3-a2,fccMobility,nullptr));             // is ( 1,-1,-1) in cartesian
//                    temp.emplace_back(new SlipSystem(a3,a2,a2-a3,fccMobility,nullptr));             // is ( 1,-1,-1) in cartesian
//
//                    temp.emplace_back(new SlipSystem(a2,a1, a2,fccMobility,nullptr));               // is (-1,-1, 1) in cartesian
//                    temp.emplace_back(new SlipSystem(a2,a1,a2*(-1),fccMobility,nullptr));           // is (-1,-1, 1) in cartesian
//                    temp.emplace_back(new SlipSystem(a2,a1, a1,fccMobility,nullptr));               // is (-1,-1, 1) in cartesian
//                    temp.emplace_back(new SlipSystem(a2,a1,a1*(-1),fccMobility,nullptr));           // is (-1,-1, 1) in cartesian
//                    temp.emplace_back(new SlipSystem(a2,a1,a2-a1,fccMobility,nullptr));             // is (-1,-1, 1) in cartesian
//                    temp.emplace_back(new SlipSystem(a2,a1,a1-a2,fccMobility,nullptr));             // is (-1,-1, 1) in cartesian
//
//                    temp.emplace_back(new SlipSystem(a1-a3,a2-a3, a1-a3,fccMobility,nullptr));      // is ( 1, 1, 1) in cartesian
//                    temp.emplace_back(new SlipSystem(a1-a3,a2-a3, a3-a1,fccMobility,nullptr));      // is ( 1, 1, 1) in cartesian
//                    temp.emplace_back(new SlipSystem(a1-a3,a2-a3,a2-a3,fccMobility,nullptr));       // is ( 1, 1, 1) in cartesian
//                    temp.emplace_back(new SlipSystem(a1-a3,a2-a3,a3-a2,fccMobility,nullptr));       // is ( 1, 1, 1) in cartesian
//                    temp.emplace_back(new SlipSystem(a1-a3,a2-a3, a1-a2,fccMobility,nullptr));      // is ( 1, 1, 1) in cartesian
//                    temp.emplace_back(new SlipSystem(a1-a3,a2-a3, a2-a1,fccMobility,nullptr));      // is ( 1, 1, 1) in cartesian
                }
//            }
            
            
            
//            if(enable110planes)
//            {// <110>{110}
//                temp.emplace_back(new SlipSystem(a1+a2-a3,a3, a3,fccMobility,nullptr));
//                temp.emplace_back(new SlipSystem(a1+a2-a3,a3, a3*(-1),fccMobility,nullptr));
//
//                temp.emplace_back(new SlipSystem(a1+a2-a3,a1-a2,a1-a2,fccMobility,nullptr));
//                temp.emplace_back(new SlipSystem(a1+a2-a3,a1-a2,a2-a1,fccMobility,nullptr));
//
//                temp.emplace_back(new SlipSystem(a1+a3-a2,a2, a2,fccMobility,nullptr));
//                temp.emplace_back(new SlipSystem(a1+a3-a2,a2, a2*(-1),fccMobility,nullptr));
//
//                temp.emplace_back(new SlipSystem(a1+a3-a2,a1-a3,a1-a3,fccMobility,nullptr));
//                temp.emplace_back(new SlipSystem(a1+a3-a2,a1-a3,a3-a1,fccMobility,nullptr));
//
//                temp.emplace_back(new SlipSystem(a2+a3-a1,a1, a1,fccMobility,nullptr));
//                temp.emplace_back(new SlipSystem(a2+a3-a1,a1, a1*(-1),fccMobility,nullptr));
//
//                temp.emplace_back(new SlipSystem(a2+a3-a1,a2-a3,a2-a3,fccMobility,nullptr));
//                temp.emplace_back(new SlipSystem(a2+a3-a1,a2-a3,a3-a2,fccMobility,nullptr));
//            }
            
            
            return temp;
        }
        
        

    std::vector<std::shared_ptr<SecondPhase<3>>> FCClattice<3>::getSecondPhases(const PolycrystallineMaterialBase& material,
                                                                                const PlaneNormalContainerType& plN)
    {
        
        const std::vector<std::string> spNames(TextFileParser(material.materialFile).readArray<std::string>("secondPhases",true));
        std::vector<std::shared_ptr<SecondPhase<3>>> temp;

        for(const std::string& sp : spNames)
        {
            if(sp=="L12")
            {
                const double APB(TextFileParser(material.materialFile).readScalar<double>("APB_SI",true)/(material.mu_SI*material.b_SI));
                const double SISF(TextFileParser(material.materialFile).readScalar<double>("SISF_SI",true)/(material.mu_SI*material.b_SI));
                const double CESF(TextFileParser(material.materialFile).readScalar<double>("CESF_SI",true)/(material.mu_SI*material.b_SI));
                const double CISF(TextFileParser(material.materialFile).readScalar<double>("CISF_SI",true)/(material.mu_SI*material.b_SI));
                const double SESF(TextFileParser(material.materialFile).readScalar<double>("SESF_SI",true)/(material.mu_SI*material.b_SI));

                const Eigen::Matrix<double,4,2> waveVectors111(0.5*(Eigen::Matrix<double,4,2>()<<0.0, 0.0,
                                                             /*                        */ 0.0, 1.0,
                                                             /*                        */ 1.0,-1.0,
                                                             /*                        */ 2.0, 0.0).finished()); // the factor 0.5 accounts for the fact that A_L12=2*A_fcc
                const Eigen::Matrix<double,6,3> f111((Eigen::Matrix<double,6,3>()<<0.00,0.0, 0.0,
                                                   /*                        */ 0.50,sqrt(3.0)/6.0, CESF,
                                                   /*                        */ 1.00,sqrt(3.0)/3.0,SISF,
                                                   /*                        */ 1.50,sqrt(3.0)/2.0, APB,
                                                   /*                        */ 2.00,2.0*sqrt(3.0)/3.0, SESF,
                                                   /*                        */ 2.50,5.0*sqrt(3.0)/6.0, CISF).finished());
                const int rotSymm111(3);
                const std::vector<Eigen::Matrix<double,2,1>> mirSymm111;

                std::map<const LatticePlaneBase*,std::shared_ptr<GammaSurface>> gsMap;
                for(const auto& n : plN)
                {
                    if(std::abs(n->planeSpacing()-sqrt(6.0)/3)<FLT_EPSILON)
                    {// a 111 plane
                        std::shared_ptr<GammaSurface> gammaSurface(new GammaSurface(*n,waveVectors111,f111,rotSymm111,mirSymm111));
                        gsMap.emplace(n.get(),gammaSurface);
                    }
                }
                temp.emplace_back(new SecondPhase<3>("L12",gsMap));
            }
            else
            {
                throw std::runtime_error("Unnown SecondPhase "+sp+" in FCC crystal.");
            }
        }
        return temp;
    }
        
} // namespace model
#endif

