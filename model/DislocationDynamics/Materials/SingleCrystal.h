/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_SingleCrystal_H_
#define model_SingleCrystal_H_

#include <cmath>
#include <string>
#include <vector>
#include <tuple>

#include <MPIcout.h> // defines mode::cout
#include <TerminalColors.h> // defines mode::cout
#include <MaterialSymmetry.h>
#include <DislocatedMaterial.h>
#include <BCClattice.h>
#include <FCClattice.h>
#include <HEXlattice.h>
#include <LatticeMath.h>
#include <SlipSystem.h>
#include <TextFileParser.h>


namespace model
{
    
    template<int dim>
    class SingleCrystal : public Lattice<dim>
    /*                 */,private std::vector<LatticePlaneBase>
    /*                 */,private std::vector<std::shared_ptr<SlipSystem>>
    {
        
        typedef Lattice<dim> LatticeType;
        typedef std::vector<LatticePlaneBase> PlaneNormalContainerType;
        typedef std::vector<std::shared_ptr<SlipSystem>> SlipSystemContainerType;
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        
        /**********************************************************************/
        static Eigen::Matrix<double,dim,dim> getLatticeBasis(const DislocatedMaterial<dim,Isotropic>& material)
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
        
        /**********************************************************************/
        static PlaneNormalContainerType getPlaneNormals(const DislocatedMaterial<dim,Isotropic>& material,
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
        
        /**********************************************************************/
        static SlipSystemContainerType getSlipSystems(const DislocatedMaterial<dim,Isotropic>& material,
                                                      const LatticeType& lat)
        {
            if(material.crystalStructure=="BCC")
            {
                return BCClattice<dim>::slipSystems(material.dislocationMobilities,lat,material);
            }
            else if(material.crystalStructure=="FCC")
            {
                const bool enablePartials(TextFileParser("inputFiles/DD.txt").readScalar<int>("enablePartials",true));
                return FCClattice<dim>::slipSystems(material.dislocationMobilities,lat,material,enablePartials);
            }
            else if(material.crystalStructure=="HEX")
            {
                const bool enablePartials(TextFileParser("inputFiles/DD.txt").readScalar<int>("enablePartials",true));
                return HEXlattice<dim>::slipSystems(material.dislocationMobilities,lat,material,enablePartials);
            }
            else
            {
                std::cout<<"Unknown crystal structure '"<<material.crystalStructure<<"'. Exiting."<<std::endl;
                exit(EXIT_FAILURE);
            }
        }
        
    public:
        

        /**********************************************************************/
        SingleCrystal(const DislocatedMaterial<dim,Isotropic>& material,const MatrixDim& C2G) :
        /* init */ LatticeType(getLatticeBasis(material),C2G)
        /* init */,PlaneNormalContainerType(getPlaneNormals(material,*this))
        /* init */,SlipSystemContainerType(getSlipSystems(material,*this))
        {
                        
        }
        
        /**********************************************************************/
        const LatticeType& lattice() const
        {
            return *this;
        }
        
        /**********************************************************************/
        const PlaneNormalContainerType& planeNormals() const
        {
            return *this;
        }
        
        /**********************************************************************/
        const SlipSystemContainerType& slipSystems() const
        {
            return *this;
        }
        
    };
    
    
}
#endif
