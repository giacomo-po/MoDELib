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

//#include <PeriodicElement.h>
//#include <CrystalOrientation.h>
#include <MPIcout.h> // defines mode::cout
#include <TerminalColors.h> // defines mode::cout
//#include <EigenDataReader.h> // defines mode::cout
#include <MaterialSymmetry.h>
//#include <MaterialBase.h>
#include <BCClattice.h>
#include <FCClattice.h>
#include <LatticeMath.h>
#include <SlipSystem.h>
#include <TextFileParser.h>
#include <DislocationMobility.h>
#include <Material.h>


namespace model
{
    
    template<int dim>
    class SingleCrystal : //public SingleCrystalBase
    /*                 */ public Lattice<dim>
    /*                 */,private std::vector<LatticePlaneBase>
    /*                 */,private std::vector<SlipSystem>
    {
        
        typedef Lattice<dim> LatticeType;
        typedef std::vector<LatticePlaneBase> PlaneNormalContainerType;
        typedef std::vector<SlipSystem> SlipSystemContainerType;
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        
        /**********************************************************************/
        static Lattice<dim> getLattice(const std::string& crystalStructure,
                                       const MatrixDim& C2G)
        {
            if(crystalStructure=="BCC")
            {
                return BCClattice<dim>(C2G);
            }
            else if(crystalStructure=="FCC")
            {
                return FCClattice<dim>(C2G);
            }
            else
            {
                std::cout<<"Unknown crystal structure '"<<crystalStructure<<"'. Exiting."<<std::endl;
                exit(EXIT_FAILURE);
            }
        }
        
        /**********************************************************************/
        static PlaneNormalContainerType getPlaneNormals(const std::string& crystalStructure,
                                                        const LatticeType& lat)
        {
            if(crystalStructure=="BCC")
            {
                return BCClattice<dim>::reciprocalPlaneNormals(lat);
            }
            else if(crystalStructure=="FCC")
            {
                return FCClattice<dim>::reciprocalPlaneNormals(lat);
            }
            else
            {
                std::cout<<"Unknown crystal structure '"<<crystalStructure<<"'. Exiting."<<std::endl;
                exit(EXIT_FAILURE);
            }
        }
        
        /**********************************************************************/
        static SlipSystemContainerType getSlipSystems(const std::string& crystalStructure,
                                                      const LatticeType& lat)
        {
            if(crystalStructure=="BCC")
            {
                return BCClattice<dim>::slipSystems(lat);
            }
            else if(crystalStructure=="FCC")
            {
                return FCClattice<dim>::slipSystems(lat);
            }
            else
            {
                std::cout<<"Unknown crystal structure '"<<crystalStructure<<"'. Exiting."<<std::endl;
                exit(EXIT_FAILURE);
            }
        }
        
    public:
        

        /**********************************************************************/
        SingleCrystal(const Material<dim,Isotropic>& material,const MatrixDim& C2G) :
        /* init */ LatticeType(getLattice(material.crystalStructure,C2G))
        /* init */,PlaneNormalContainerType(getPlaneNormals(material.crystalStructure,*this))
        /* init */,SlipSystemContainerType(getSlipSystems(material.crystalStructure,*this))
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
