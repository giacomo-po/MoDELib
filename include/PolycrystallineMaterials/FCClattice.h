/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_FCClattice_H_
#define model_FCClattice_H_

#include <memory>
#include <vector>
#include <Eigen/Dense>

#include <LatticeModule.h>
#include <SlipSystem.h>
#include <PolycrystallineMaterialBase.h>
#include <DislocationMobilityFCC.h>
#include <RationalLatticeDirection.h>

namespace model
{
    
    template<int dim>
    struct FCClattice
    {
        
    };
    
    template<>
    struct FCClattice<3> : public Lattice<3>
    {
        static constexpr int dim=3;
        static constexpr auto name="FCC";
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        
        static constexpr bool enable111planes=true;
        static constexpr bool enable110planes=false;
        
        FCClattice(const MatrixDim& Q);
        static Eigen::Matrix<double,dim,dim> getLatticeBasis();
        static std::vector<LatticePlaneBase> reciprocalPlaneNormals(const Lattice<dim>& lat);
        static std::vector<std::shared_ptr<SlipSystem>> slipSystems(const std::map<std::string,std::shared_ptr<DislocationMobilityBase>>& mobilities,
                                                                    const Lattice<dim>& lat,
                                                                    const PolycrystallineMaterialBase& material,
                                                                    const bool& enablePartials);
        
        
        
    };
    
} // namespace model
#endif

