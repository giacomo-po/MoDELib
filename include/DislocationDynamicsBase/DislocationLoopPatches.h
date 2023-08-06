/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2018 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2018 by Yinan Cui <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationLoopPatches_H_
#define model_DislocationLoopPatches_H_

#include <map>
#include <vector>
#include <memory>


#include <Eigen/Dense>
#include <PeriodicGlidePlane.h>
#include <PolycrystallineMaterialBase.h>
#include <DislocationFieldBase.h>
#include <SutherlandHodgman.h>

namespace model
{
    
    template <int dim>
    struct DislocationLoopPatches
//: public std::map<std::shared_ptr<PeriodicPlanePatch<dim>>,std::vector<Eigen::Matrix<double,dim-1,1>>>
//    /*                           */,public std::map<std::shared_ptr<PeriodicPlanePatch<dim>>,std::vector<Eigen::Matrix<double,dim,1>>>
    {
        typedef Eigen::Matrix<double,dim-1,1> VectorLowerDim;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef PeriodicGlidePlane<dim> PeriodicGlidePlaneType;
        typedef std::map<std::shared_ptr<PeriodicPlanePatch<dim>>,std::vector<VectorLowerDim>> LocalPatchPositionType;
        typedef std::map<std::shared_ptr<PeriodicPlanePatch<dim>>,std::vector<VectorDim>> GlobalPatchPositionType;
        
        
        
    public:
        
        const std::shared_ptr<PeriodicGlidePlaneType> periodicGlidePlane;

    private:
        
        LocalPatchPositionType _localPatches;
        GlobalPatchPositionType _globalPatches;
    public:
        
        DislocationLoopPatches(const std::shared_ptr<PeriodicGlidePlaneType>& periodicGlidePlane_in);
        void update(const std::vector<VectorDim>& linkShifts, const std::vector<VectorDim>& nodePos);
        const  LocalPatchPositionType&  localPatches() const;
        const GlobalPatchPositionType& globalPatches() const;
        LocalPatchPositionType&  localPatches();
        GlobalPatchPositionType& globalPatches();
        const VectorDim orientedArea() const;
        double solidAngle(const VectorDim& x) const;
        template <typename T> static int sgn(const T& val);

	};	
	
}
#endif

