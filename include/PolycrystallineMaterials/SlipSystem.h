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
#include <GlidePlaneNoise.h>

namespace model
{


    struct SlipSystem : public StaticID<SlipSystem>
    {
        
        typedef Eigen::Matrix<double,3,1> VectorDim;
        typedef Eigen::Matrix<double,3,3> MatrixDim;
        
        const LatticePlaneBase n;
        const RationalLatticeDirection<3>  s;
        const VectorDim unitNormal;
        const VectorDim unitSlip;
        const VectorDim unitSlipFull;
        const MatrixDim G2Lfull;
        const std::shared_ptr<DislocationMobilityBase> mobility;
        const std::shared_ptr<GammaSurface> gammaSurface;
        const std::shared_ptr<GlidePlaneNoise> planeNoise;
                
        SlipSystem(const LatticePlaneBase& n_in,
                   const LatticeVector<3>& slip_in,
                   const std::shared_ptr<DislocationMobilityBase>& mobility_in,
                   const std::shared_ptr<GammaSurface>& gammaSurface_in,
                   const std::shared_ptr<GlidePlaneNoise>& planeNoise_in);
        
        SlipSystem(const LatticePlaneBase& n_in,
                   const RationalLatticeDirection<3>& slip_in,
                   const std::shared_ptr<DislocationMobilityBase>& mobility_in,
                   const std::shared_ptr<GammaSurface>& gammaSurface_in,
                   const std::shared_ptr<GlidePlaneNoise>& planeNoise_in);
        
        bool isPartial() const;
        bool isSameAs(const RationalLatticeDirection<3>& s1,const ReciprocalLatticeDirection<3>& n1);
        double misfitEnergy(const Eigen::Matrix<double,3,1>& b);
        
        
        Eigen::Matrix<double,2,1> globalToLocal(const VectorDim& x) const;
        VectorDim localToGlobal(const Eigen::Matrix<double,2,1>& x) const;
        std::tuple<Eigen::Matrix<double,3,3>,double,double> gridInterp(const VectorDim& x) const ;
        std::tuple<Eigen::Matrix<double,3,3>,double,double> gridVal(const Eigen::Array<int,2,1>& idx) const;

    };

}
#endif
