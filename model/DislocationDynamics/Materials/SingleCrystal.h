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

//#include <model/DislocationDynamics/Materials/PeriodicElement.h>
//#include <model/DislocationDynamics/Materials/CrystalOrientation.h>
#include <model/MPI/MPIcout.h> // defines mode::cout
#include <model/Utilities/TerminalColors.h> // defines mode::cout
//#include <model/IO/EigenDataReader.h> // defines mode::cout
#include <model/DislocationDynamics/Materials/MaterialSymmetry.h>
//#include <model/DislocationDynamics/Materials/MaterialBase.h>
#include <model/DislocationDynamics/Materials/BCClattice.h>
#include <model/DislocationDynamics/Materials/FCClattice.h>
#include <model/LatticeMath/LatticeMath.h>
#include <model/DislocationDynamics/Materials/SlipSystem.h>
#include <model/IO/TextFileParser.h>



//#include <model/IO/TextFileParser.h>



namespace model
{
    
    
    class SingleCrystalBase : private std::tuple<std::string
    /*                                        */,std::string>
    {
        
        typedef std::tuple<std::string
        /*              */,std::string> BaseType;
        
        static BaseType init(const std::string& fileName)
        {
            TextFileParser parser(fileName);
            return std::make_tuple(fileName
                                   ,parser.readString("crystalStructure")
                                   );
            
        }
        
    public:
        const std::string& materialFile;
        const std::string& crystalStructure;
        
        SingleCrystalBase(const std::string& fileName) :
        /* init */ BaseType(init(fileName))
        /* init */,materialFile(std::get<0>(*this))
        /* init */,crystalStructure(std::get<1>(*this))
        {
            
        }
        
    };
    
    template<int dim>
    class SingleCrystal : public SingleCrystalBase
    /*                 */,public Lattice<dim>
    /*                 */,private std::vector<LatticePlaneBase>
    /*                 */,private std::vector<SlipSystem>
    {
        
        typedef Lattice<dim> LatticeType;
        typedef std::vector<LatticePlaneBase> PlaneNormalContainerType;
        typedef std::vector<SlipSystem> SlipSystemContainerType;
        
        /**********************************************************************/
        static Lattice<dim> getLattice(const std::string& crystalStructure)
        {
            if(crystalStructure=="BCC")
            {
                return BCClattice<dim>();
            }
            else if(crystalStructure=="FCC")
            {
                return FCClattice<dim>();
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
        SingleCrystal(const std::string& materialName) :
        /* init */ SingleCrystalBase(materialName)
        /* init */,LatticeType(getLattice(this->crystalStructure))
        /* init */,PlaneNormalContainerType(getPlaneNormals(this->crystalStructure,*this))
        /* init */,SlipSystemContainerType(getSlipSystems(this->crystalStructure,*this))
        {
            
            model::cout<<greenColor<<"Creating SingleCrystal:"<<defaultColor<<std::endl;
            model::cout<<"  material="<<this->materialFile<<std::endl;
            model::cout<<"  crystal structure="<<this->crystalStructure<<std::endl;
            model::cout<<"  # plane normals="<<planeNormals().size()<<std::endl;
            model::cout<<"  # slip systems="<<slipSystems().size()<<std::endl;
            
            //            PROBLEM IS COPYING LATTICE VECTORS IN SLIP SYSTEM, BECAUSE THEY HOLD A REFERENCE TO A LATTICE
            
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
