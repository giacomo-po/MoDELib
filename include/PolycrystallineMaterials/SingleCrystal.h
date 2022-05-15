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

#include <TerminalColors.h> // defines mode::cout
#include <MaterialSymmetry.h>
#include <PolycrystallineMaterial.h>
#include <BCClattice.h>
#include <FCClattice.h>
#include <HEXlattice.h>
#include <LatticeModule.h>
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

        static Eigen::Matrix<double,dim,dim> getLatticeBasis(const PolycrystallineMaterial<dim,Isotropic>& material);
        static PlaneNormalContainerType getPlaneNormals(const PolycrystallineMaterial<dim,Isotropic>& material,const LatticeType& lat);
        static SlipSystemContainerType getSlipSystems(const PolycrystallineMaterial<dim,Isotropic>& material,
                                                      const LatticeType& lat,
                                                      const std::string& polyFile);
        
    public:

        SingleCrystal(const PolycrystallineMaterial<dim,Isotropic>& material,const MatrixDim& C2G,const std::string& polyFile);
        const LatticeType& lattice() const;
        const PlaneNormalContainerType& planeNormals() const;
        const SlipSystemContainerType& slipSystems() const;
        
    };
    
    
}
#endif
