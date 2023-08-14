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
#include <SingleCrystalBase.h>
#include <DislocationMobilitySelector.h>
//#include <SecondPhase.h>

namespace model
{
    
    template<int dim>
    struct FCClattice
    {
        
    };
    
    template<>
    struct FCClattice<3> : public SingleCrystalBase<3>
    /*                  */,private SingleCrystalBase<3>::PlaneNormalContainerType
    /*                  */,private SingleCrystalBase<3>::SlipSystemContainerType
    /*                  */,private SingleCrystalBase<3>::SecondPhaseContainerType
    {
//        static constexpr auto name="FCC";
        static constexpr int dim=3;
        typedef typename SingleCrystalBase<dim>::MatrixDim MatrixDim;
        typedef typename SingleCrystalBase<dim>::PlaneNormalContainerType PlaneNormalContainerType;
        typedef typename SingleCrystalBase<dim>::SlipSystemContainerType SlipSystemContainerType;
        typedef typename SingleCrystalBase<dim>::SecondPhaseContainerType SecondPhaseContainerType;

        
        static constexpr bool enable111planes=true;
        static constexpr bool enable110planes=false;
        
        FCClattice(const MatrixDim& Q,const PolycrystallineMaterialBase& material,const std::string& polyFile);
        static Eigen::Matrix<double,dim,dim> getLatticeBasis();
        std::vector<std::shared_ptr<LatticePlaneBase>> getPlaneNormals() const;
        std::vector<std::shared_ptr<SlipSystem>> getSlipSystems(const PolycrystallineMaterialBase& material,const std::string& polyFile,const PlaneNormalContainerType& plN) const;
        std::vector<std::shared_ptr<SecondPhase<dim>>> getSecondPhases(const PolycrystallineMaterialBase& material,const SlipSystemContainerType& slipSystems) const ;
        
        const PlaneNormalContainerType& planeNormals() const override;
        const SlipSystemContainerType& slipSystems() const override;
        const SecondPhaseContainerType& secondPhases() const override;
        
    };
    
} // namespace model
#endif

