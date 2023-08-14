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
#include <GlidePlaneNoise.h>


namespace model
{
    FCClattice<3>::FCClattice(const MatrixDim& Q,const PolycrystallineMaterialBase& material,const std::string& polyFile) :
    /* init */ SingleCrystalBase<dim>(getLatticeBasis(),Q)
    /* init */,PlaneNormalContainerType(getPlaneNormals())
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
        
        return temp;
    }

    std::vector<std::shared_ptr<SlipSystem>> FCClattice<3>::getSlipSystems(const PolycrystallineMaterialBase& material,
                                                                           const std::string& polyFile,
                                                                           const PlaneNormalContainerType& plN) const
    {/*!\returns a std::vector of ReciprocalLatticeDirection(s) corresponding
      * the slip systems of the Hexagonal lattice
      */
        
        const std::string dislocationMobilityType(TextFileParser(polyFile).readString("dislocationMobilityType",true));
        DislocationMobilitySelector mobilitySelector("FCC");
        const std::shared_ptr<DislocationMobilityBase> fccMobility(mobilitySelector.getMobility(dislocationMobilityType,material));
//
//        const std::shared_ptr<DislocationMobilityBase> fccMobility(new DislocationMobilityFCC(material));
        
        std::vector<std::shared_ptr<SlipSystem>> temp;
        
        const bool enablePartials(TextFileParser(polyFile).readScalar<int>("enablePartials",true));
        
        typedef Eigen::Matrix<double,dim,1> VectorDimD;
        
        
        const int solidSolutionNoiseMode(TextFileParser(polyFile).readScalar<int>("solidSolutionNoiseMode",true));
        const int stackingFaultNoiseMode(TextFileParser(polyFile).readScalar<int>("stackingFaultNoiseMode",true));

        std::shared_ptr<GlidePlaneNoise> planeNoise((solidSolutionNoiseMode||stackingFaultNoiseMode)? new GlidePlaneNoise(polyFile,material) : nullptr);
        
//        /* init */,planeNoise((TextFileParser(simulationParameters.traitsIO.noiseFile).readScalar<int>("solidSolutionNoiseMode") || TextFileParser(simulationParameters.traitsIO.noiseFile).readScalar<int>("stackingFaultNoiseMode"))? new GlidePlaneNoise(simulationParameters.traitsIO,poly) : nullptr)

        
        
        const double d111(this->reciprocalLatticeDirection(this->C2G*(VectorDimD()<<1.0,1.0,1.0).finished()).planeSpacing());
        
        if(enablePartials)
        {
            
            const double ISF(TextFileParser(material.materialFile).readScalar<double>("ISF_SI",true)/(material.mu_SI*material.b_SI));
            const double USF(TextFileParser(material.materialFile).readScalar<double>("USF_SI",true)/(material.mu_SI*material.b_SI));
            const double MSF(TextFileParser(material.materialFile).readScalar<double>("MSF_SI",true)/(material.mu_SI*material.b_SI));
            
            const Eigen::Matrix<double,3,2> waveVectors((Eigen::Matrix<double,3,2>()<<0.0, 0.0,
                                                         /*                        */ 0.0, 1.0,
//                                                         /*                        */ 1.0,-1.0
                                                         /*                        */ 1.0,1.0
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
                    const auto& a1(planeBase->primitiveVectors.first);
                    const auto& a3(planeBase->primitiveVectors.second);

                    const auto b1(a1);
                    const auto b2(a3-a1);
                    const auto b3(a3*(-1));
//                    std::shared_ptr<GammaSurface> gammaSurface(new GammaSurface(*planeBase,waveVectors,f,rotSymm,mirSymm));
                    std::shared_ptr<GammaSurface> gammaSurface(new GammaSurface(b1,b2,waveVectors,f,rotSymm,mirSymm));
                    temp.emplace_back(new SlipSystem(LatticePlaneBase(b1,b2), RationalLatticeDirection<3>(Rational(1,3),b1-b3),fccMobility,gammaSurface,planeNoise));               // is (-1, 1,-1) in cartesian
                    temp.emplace_back(new SlipSystem(LatticePlaneBase(b1,b2), RationalLatticeDirection<3>(Rational(1,3),b1-b2),fccMobility,gammaSurface,planeNoise));               // is (-1, 1,-1) in cartesian
                    temp.emplace_back(new SlipSystem(LatticePlaneBase(b2,b3), RationalLatticeDirection<3>(Rational(1,3),b2-b1),fccMobility,gammaSurface,planeNoise));               // is (-1, 1,-1) in cartesian
                    temp.emplace_back(new SlipSystem(LatticePlaneBase(b2,b3), RationalLatticeDirection<3>(Rational(1,3),b2-b3),fccMobility,gammaSurface,planeNoise));               // is (-1, 1,-1) in cartesian
                    temp.emplace_back(new SlipSystem(LatticePlaneBase(b3,b1), RationalLatticeDirection<3>(Rational(1,3),b3-b2),fccMobility,gammaSurface,planeNoise));               // is (-1, 1,-1) in cartesian
                    temp.emplace_back(new SlipSystem(LatticePlaneBase(b3,b1), RationalLatticeDirection<3>(Rational(1,3),b3-b1),fccMobility,gammaSurface,planeNoise));               // is (-1, 1,-1) in cartesian
                }
            }
        }
        else
        {
            
            for(const auto& planeBase : plN)
            {
                if(std::fabs(planeBase->planeSpacing()-d111)<FLT_EPSILON)
                {// a {111} plane
                    const auto& a1(planeBase->primitiveVectors.first);
                    const auto& a3(planeBase->primitiveVectors.second);
                    
                    const auto b1(a1);
                    const auto b2(a3-a1);
                    const auto b3(a3*(-1));

                    temp.emplace_back(new SlipSystem(LatticePlaneBase(b1,b2), b1,fccMobility,nullptr,planeNoise));
                    temp.emplace_back(new SlipSystem(LatticePlaneBase(b1,b2),b1*(-1),fccMobility,nullptr,planeNoise));
                    temp.emplace_back(new SlipSystem(LatticePlaneBase(b2,b3), b2,fccMobility,nullptr,planeNoise));
                    temp.emplace_back(new SlipSystem(LatticePlaneBase(b2,b3),b2*(-1),fccMobility,nullptr,planeNoise));
                    temp.emplace_back(new SlipSystem(LatticePlaneBase(b3,b1),b3,fccMobility,nullptr,planeNoise));
                    temp.emplace_back(new SlipSystem(LatticePlaneBase(b3,b1),b3*(-1),fccMobility,nullptr,planeNoise));

                }
            }
        }
        
        return temp;
    }



    std::vector<std::shared_ptr<SecondPhase<3>>> FCClattice<3>::getSecondPhases(const PolycrystallineMaterialBase& material,
                                                                                const SlipSystemContainerType& slipSystems) const
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
//                                                                    /*                        */ 1.0,-1.0,
                                                                    /*                        */ 1.0,1.0,
                                                                    /*                        */ 2.0, 0.0).finished()); // the factor 0.5 accounts for the fact that A_L12=2*A_fcc
                const Eigen::Matrix<double,6,3> f111((Eigen::Matrix<double,6,3>()<<0.00,0.0, 0.0,
                                                      /*                        */ 0.50,sqrt(3.0)/6.0, CISF,
                                                      /*                        */ 1.00,sqrt(3.0)/3.0,SESF,
                                                      /*                        */ 1.50,sqrt(3.0)/2.0, APB,
                                                      /*                        */ 2.00,2.0*sqrt(3.0)/3.0, SISF,
                                                      /*                        */ 2.50,5.0*sqrt(3.0)/6.0, CESF).finished());
                const int rotSymm111(3);
                const std::vector<Eigen::Matrix<double,2,1>> mirSymm111;
                
                std::map<std::shared_ptr<SlipSystem>,std::shared_ptr<GammaSurface>> gsMap;
                for(const auto& ss : slipSystems)
                {
                    if(std::abs(ss->n.planeSpacing()-sqrt(6.0)/3)<FLT_EPSILON)
                    {// a 111 plane
                        const auto& b1(ss->n.primitiveVectors.first);
                        const auto& b2(ss->n.primitiveVectors.second);
//                        std::shared_ptr<GammaSurface> gammaSurface(new GammaSurface(ss->n,waveVectors111,f111,rotSymm111,mirSymm111));
                        std::shared_ptr<GammaSurface> gammaSurface(new GammaSurface(b1,b2,waveVectors111,f111,rotSymm111,mirSymm111));

                        gsMap.emplace(ss,gammaSurface);
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

