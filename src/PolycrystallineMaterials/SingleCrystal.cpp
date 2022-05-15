/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_SingleCrystal_cpp_
#define model_SingleCrystal_cpp_

#include <cmath>
#include <string>
#include <vector>
#include <tuple>

// defines mode::cout
#include <TerminalColors.h> // defines mode::cout
#include <MaterialSymmetry.h>
#include <PolycrystallineMaterial.h>
#include <BCClattice.h>
#include <FCClattice.h>
#include <HEXlattice.h>
#include <LatticeModule.h>
#include <SlipSystem.h>
#include <TextFileParser.h>
#include <SingleCrystal.h>

namespace model
{

    template<int dim>
    Eigen::Matrix<double,dim,dim> SingleCrystal<dim>::getLatticeBasis(const PolycrystallineMaterial<dim,Isotropic>& material)
    {
        if(material.crystalStructure=="BCC")
        {
            return BCClattice<dim>::getLatticeBasis();
        }
        else if(material.crystalStructure=="FCC")
        {
            return FCClattice<dim>::getLatticeBasis();
        }
        else if(material.crystalStructure=="HEX")
        {
            return HEXlattice<dim>::getLatticeBasis();
        }
        else
        {
            std::cout<<"Unknown crystal structure '"<<material.crystalStructure<<"'. Exiting."<<std::endl;
            exit(EXIT_FAILURE);
        }
    }

    template<int dim>
    typename SingleCrystal<dim>::PlaneNormalContainerType SingleCrystal<dim>::getPlaneNormals(const PolycrystallineMaterial<dim,Isotropic>& material,
                                                                                              const LatticeType& lat)
    {
        if(material.crystalStructure=="BCC")
        {
            return BCClattice<dim>::reciprocalPlaneNormals(lat);
        }
        else if(material.crystalStructure=="FCC")
        {
            return FCClattice<dim>::reciprocalPlaneNormals(lat);
        }
        else if(material.crystalStructure=="HEX")
        {
            return HEXlattice<dim>::reciprocalPlaneNormals(material,lat);
        }
        else
        {
            std::cout<<"Unknown crystal structure '"<<material.crystalStructure<<"'. Exiting."<<std::endl;
            exit(EXIT_FAILURE);
        }
    }

    template<int dim>
    typename SingleCrystal<dim>::SlipSystemContainerType SingleCrystal<dim>::getSlipSystems(const PolycrystallineMaterial<dim,Isotropic>& material,
                                                                                            const LatticeType& lat,
                                                                                            const std::string& polyFile)
    {
        if(material.crystalStructure=="BCC")
        {
            return BCClattice<dim>::slipSystems(material.dislocationMobilities,lat,material);
        }
        else if(material.crystalStructure=="FCC")
        {
            const bool enablePartials(TextFileParser(polyFile).readScalar<int>("enablePartials",true));
            return FCClattice<dim>::slipSystems(material.dislocationMobilities,lat,material,enablePartials);
        }
        else if(material.crystalStructure=="HEX")
        {
            const bool enablePartials(TextFileParser(polyFile).readScalar<int>("enablePartials",true));
            return HEXlattice<dim>::slipSystems(material.dislocationMobilities,lat,material,enablePartials);
        }
        else
        {
            std::cout<<"Unknown crystal structure '"<<material.crystalStructure<<"'. Exiting."<<std::endl;
            exit(EXIT_FAILURE);
        }
    }

    template<int dim>
    SingleCrystal<dim>::SingleCrystal(const PolycrystallineMaterial<dim,Isotropic>& material,const MatrixDim& C2G,const std::string& polyFile) :
    /* init */ LatticeType(getLatticeBasis(material),C2G)
    /* init */,PlaneNormalContainerType(getPlaneNormals(material,*this))
    /* init */,SlipSystemContainerType(getSlipSystems(material,*this,polyFile))
    {
        
    }

    template<int dim>
    const typename SingleCrystal<dim>::LatticeType& SingleCrystal<dim>::lattice() const
    {
        return *this;
    }

    template<int dim>
    const typename SingleCrystal<dim>::PlaneNormalContainerType& SingleCrystal<dim>::planeNormals() const
    {
        return *this;
    }

    template<int dim>
    const typename SingleCrystal<dim>::SlipSystemContainerType& SingleCrystal<dim>::slipSystems() const
    {
        return *this;
    }

template class SingleCrystal<3>;

}
#endif
