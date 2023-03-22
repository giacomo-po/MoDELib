/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_SingleCrystalBase_cpp_
#define model_SingleCrystalBase_cpp_

#include <cmath>
#include <string>
#include <vector>
#include <tuple>

// defines mode::cout
#include <TerminalColors.h> // defines mode::cout
//#include <MaterialSymmetry.h>
//#include <PolycrystallineMaterial.h>
//#include <BCClattice.h>
//#include <FCClattice.h>
//#include <HEXlattice.h>
//#include <LatticeModule.h>
//#include <SlipSystem.h>
//#include <TextFileParser.h>
#include <SingleCrystalBase.h>

namespace model
{

//    template<int dim>
//    Eigen::Matrix<double,dim,dim> SingleCrystalBase<dim>::getLatticeBasis(const PolycrystallineMaterial<dim,Isotropic>& material)
//    {
//        if(material.crystalStructure=="BCC")
//        {
//            return BCClattice<dim>::getLatticeBasis();
//        }
//        else if(material.crystalStructure=="FCC")
//        {
//            return FCClattice<dim>::getLatticeBasis();
//        }
//        else if(material.crystalStructure=="HEX")
//        {
//            return HEXlattice<dim>::getLatticeBasis();
//        }
//        else
//        {
//            std::cout<<"Unknown crystal structure '"<<material.crystalStructure<<"'. Exiting."<<std::endl;
//            exit(EXIT_FAILURE);
//        }
//    }
//
//    template<int dim>
//    typename SingleCrystalBase<dim>::PlaneNormalContainerType SingleCrystalBase<dim>::getPlaneNormals(const PolycrystallineMaterial<dim,Isotropic>& material,
//                                                                                              const LatticeType& lat)
//    {
//        if(material.crystalStructure=="BCC")
//        {
//            return BCClattice<dim>::reciprocalPlaneNormals(lat);
//        }
//        else if(material.crystalStructure=="FCC")
//        {
//            return FCClattice<dim>::reciprocalPlaneNormals(lat);
//        }
//        else if(material.crystalStructure=="HEX")
//        {
//            return HEXlattice<dim>::reciprocalPlaneNormals(material,lat);
//        }
//        else
//        {
//            std::cout<<"Unknown crystal structure '"<<material.crystalStructure<<"'. Exiting."<<std::endl;
//            exit(EXIT_FAILURE);
//        }
//    }
//
//    template<int dim>
//    typename SingleCrystalBase<dim>::SlipSystemContainerType SingleCrystalBase<dim>::getSlipSystems(const PolycrystallineMaterial<dim,Isotropic>& material,
//                                                                                            const LatticeType& lat,
//                                                                                            const std::string& polyFile)
//    {
//        if(material.crystalStructure=="BCC")
//        {
//            return BCClattice<dim>::slipSystems(material.dislocationMobilities,lat,material);
//        }
//        else if(material.crystalStructure=="FCC")
//        {
//            const bool enablePartials(TextFileParser(polyFile).readScalar<int>("enablePartials",true));
//            return FCClattice<dim>::slipSystems(material.dislocationMobilities,lat,material,enablePartials);
//        }
//        else if(material.crystalStructure=="HEX")
//        {
//            const bool enablePartials(TextFileParser(polyFile).readScalar<int>("enablePartials",true));
//            return HEXlattice<dim>::slipSystems(material.dislocationMobilities,lat,material,enablePartials);
//        }
//        else
//        {
//            std::cout<<"Unknown crystal structure '"<<material.crystalStructure<<"'. Exiting."<<std::endl;
//            exit(EXIT_FAILURE);
//        }
//    }
//
//    template<int dim>
//    typename SingleCrystalBase<dim>::SecondPhaseContainerType SingleCrystalBase<dim>::getSecondPhases(const PolycrystallineMaterial<dim,Isotropic>& material,
//                                                                                            const SlipSystemContainerType& slipSystems)
//    {
//        if(material.crystalStructure=="BCC")
//        {
////            return BCClattice<dim>::secondPhases(material.dislocationMobilities,lat,material);
//        }
//        else if(material.crystalStructure=="FCC")
//        {
////            const bool enablePartials(TextFileParser(polyFile).readScalar<int>("enablePartials",true));
//            return FCClattice<dim>::secondPhases(material,slipSystems);
//        }
//        else if(material.crystalStructure=="HEX")
//        {
////            const bool enablePartials(TextFileParser(polyFile).readScalar<int>("enablePartials",true));
////            return HEXlattice<dim>::secondPhases(material.dislocationMobilities,lat,material,enablePartials);
//        }
//        else
//        {
//            std::cout<<"Unknown crystal structure '"<<material.crystalStructure<<"'. Exiting."<<std::endl;
//            exit(EXIT_FAILURE);
//        }
//    }

    template<int dim>
    SingleCrystalBase<dim>::SingleCrystalBase(const MatrixDim& A,
                                      const MatrixDim& C2G) :
    /* init */ LatticeType(A,C2G)
//    /* init */,bccLattice(this->crystalStructure=="BCC"? new BCClattice<3>(*this) : nullptr)
//    /* init */,fccLattice(this->crystalStructure=="FCC"? new FCClattice<3>(*this) : nullptr)
//    /* init */,hcpLattice(this->crystalStructure=="HCP"? new HCPlattice<3>(*this) : nullptr)

//    /* init */,planeNormals(getPlaneNormals(material,*this))
//    /* init */,slipSystems(getSlipSystems(material,*this,polyFile))
//    /* init */,secondPhases(getSecondPhases(material,*this))
    {
        
    }

template<int dim>
SingleCrystalBase<dim>::~SingleCrystalBase(){}
//    template<int dim>
//    const typename SingleCrystalBase<dim>::LatticeType& SingleCrystalBase<dim>::lattice() const
//    {
//        return *this;
//    }

//    template<int dim>
//    const typename SingleCrystalBase<dim>::PlaneNormalContainerType& SingleCrystalBase<dim>::planeNormals() const
//    {
//        if(bccLattice)
//        {
//            return bccLattice->planeNormals;
//        }
//        else if(fccLattice)
//        {
//            return fccLattice->planeNormals;
//        }
//        else if(hcpLattice)
//        {
//            return hcpLattice->planeNormals;
//        }
//        else
//        {
//            throw std::runtime_error("SingleCrystal::planeNormals() no available lattice.");
//            return typename SingleCrystalBase<dim>::PlaneNormalContainerType();
//        }
//    }
//
//    template<int dim>
//    const typename SingleCrystalBase<dim>::SlipSystemContainerType& SingleCrystalBase<dim>::slipSystems() const
//    {
//        if(bccLattice)
//        {
//            return bccLattice->slipSystems;
//        }
//        else if(fccLattice)
//        {
//            return fccLattice->slipSystems;
//        }
//        else if(hcpLattice)
//        {
//            return hcpLattice->slipSystems;
//        }
//        else
//        {
//            throw std::runtime_error("SingleCrystal::planeNormals() no available lattice.");
//            return typename SingleCrystalBase<dim>::SlipSystemContainerType();
//        }
//    }
//
//    template<int dim>
//    const typename SingleCrystalBase<dim>::SecondPhaseContainerType& SingleCrystalBase<dim>::secondPhases() const
//    {
//        if(bccLattice)
//        {
//            return bccLattice->secondPhases;
//        }
//        else if(fccLattice)
//        {
//            return fccLattice->secondPhases;
//        }
//        else if(hcpLattice)
//        {
//            return hcpLattice->secondPhases;
//        }
//        else
//        {
//            throw std::runtime_error("SingleCrystal::planeNormals() no available lattice.");
//            return typename SingleCrystalBase<dim>::SecondPhaseContainerType();
//        }
//    }

template struct SingleCrystalBase<3>;

}
#endif
