/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_SingleCrystalBase_H_
#define model_SingleCrystalBase_H_

#include <cmath>
#include <string>
#include <vector>
#include <tuple>

#include <TerminalColors.h> // defines mode::cout
//#include <BCClattice.h>
//#include <FCClattice.h>
//#include <HEXlattice.h>
#include <LatticeModule.h>
#include <SlipSystem.h>
#include <SecondPhase.h>


#include <TextFileParser.h>


namespace model
{
    
    template<int dim>
    struct SingleCrystalBase : public Lattice<dim>
//    /*                 */,private std::vector<std::shared_ptr<LatticePlaneBase>>
//    /*                 */,private std::vector<std::shared_ptr<SlipSystem>>
//    /*                 */,private std::vector<std::shared_ptr<SecondPhase<dim>>>
    {
        
        typedef Lattice<dim> LatticeType;
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        typedef std::vector<std::shared_ptr<LatticePlaneBase>> PlaneNormalContainerType;
        typedef std::vector<std::shared_ptr<SlipSystem>> SlipSystemContainerType;
        typedef std::vector<std::shared_ptr<SecondPhase<dim>>> SecondPhaseContainerType;
//        typedef std::map<std::string,std::shared_ptr<DislocationMobilityBase>> DislocationMobilityContainerType;
       

//        static Eigen::Matrix<double,dim,dim> getLatticeBasis(const PolycrystallineMaterial<dim,Isotropic>& material);
//        static PlaneNormalContainerType getPlaneNormals(const PolycrystallineMaterial<dim,Isotropic>& material,const LatticeType& lat);
//        static SlipSystemContainerType getSlipSystems(const PolycrystallineMaterial<dim,Isotropic>& material,
//                                                      const LatticeType& lat,
//                                                      const std::string& polyFile);
//
//        static SecondPhaseContainerType getSecondPhases(const PolycrystallineMaterial<dim,Isotropic>& material,
//                                                        const SlipSystemContainerType& slipSystems);
        
        
//        const std::shared_ptr<BCClattice<dim>> bccLattice;
//        const std::shared_ptr<FCClattice<dim>> fccLattice;
//        const std::shared_ptr<HEXlattice<dim>> hexLattice;
        

        SingleCrystalBase(const MatrixDim& A,
                          const MatrixDim& C2G);
        
//        const LatticeType& lattice() const;
        
        virtual const PlaneNormalContainerType& planeNormals() const =0;
        virtual const SlipSystemContainerType& slipSystems() const =0;
        virtual const SecondPhaseContainerType& secondPhases() const =0;
//        virtual const DislocationMobilityContainerType& dislocationMobilities() const =0;
        
        
    };
    
    
}
#endif
