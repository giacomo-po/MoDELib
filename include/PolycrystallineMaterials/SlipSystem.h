/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_SLIPSYSTEM_H_
#define model_SLIPSYSTEM_H_

#include <memory>
#include <assert.h>
#include <LatticeModule.h>
//#include <LatticePlaneBase.h>
//#include <LatticeVector.h>
//#include <RationalLatticeDirection.h>
#include <DislocationMobilityBase.h>
#include <GammaSurface.h>

namespace model
{
    
    
    struct SlipSystem : public StaticID<SlipSystem>
    {
        
        const LatticePlaneBase& n;
        const RationalLatticeDirection<3>  s;
        const Eigen::Matrix<double,3,1>  unitNormal;
        const Eigen::Matrix<double,3,1>  unitSlip;
        const Eigen::Matrix<double,3,3>  unitTensorSN;
        const Eigen::Matrix<double,3,3>  unitTensorSNsym;
        const Eigen::Matrix<double,3,3>  unitTensorOrthSN;
        const Eigen::Matrix<double,3,3>  unitTensorOrthSNsym;
        const std::shared_ptr<DislocationMobilityBase> mobility;
        const std::shared_ptr<GammaSurface> gammaSurface;
        
        
//        SlipSystem(const LatticeVector<3>& a1,
//                   const LatticeVector<3>& a2,
//                   const LatticeVector<3>& slip_in,
//                   const std::shared_ptr<DislocationMobilityBase>& mobility_in,
//                   const std::shared_ptr<GammaSurface>& gammaSurface_in);
//
//        SlipSystem(const LatticeVector<3>& a1,
//                   const LatticeVector<3>& a2,
//                   const RationalLatticeDirection<3>& slip_in,
//                   const std::shared_ptr<DislocationMobilityBase>& mobility_in,
//                   const std::shared_ptr<GammaSurface>& gammaSurface_in);
        
        SlipSystem(const LatticePlaneBase& n_in,
                   const LatticeVector<3>& slip_in,
                   const std::shared_ptr<DislocationMobilityBase>& mobility_in,
                   const std::shared_ptr<GammaSurface>& gammaSurface_in);
        
        SlipSystem(const LatticePlaneBase& n_in,
                   const RationalLatticeDirection<3>& slip_in,
                   const std::shared_ptr<DislocationMobilityBase>& mobility_in,
                   const std::shared_ptr<GammaSurface>& gammaSurface_in);
        
        bool isPartial() const;
        bool isSameAs(const RationalLatticeDirection<3>& s1,const ReciprocalLatticeDirection<3>& n1);
        double misfitEnergy(const Eigen::Matrix<double,3,1>& b);
        
    };

}
#endif
