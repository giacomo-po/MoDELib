/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_SecondPhase_H_
#define model_SecondPhase_H_

#include <memory>
#include <map>
#include <assert.h>

#include <LatticeModule.h>
//#include <LatticePlaneBase.h>
//#include <LatticeVector.h>
//#include <RationalLatticeDirection.h>
#include <SlipSystem.h>

namespace model
{
    
    template<int dim>
    struct SecondPhase
    {
        
        const std::string name;
        const std::map<std::shared_ptr<SlipSystem>,std::shared_ptr<GammaSurface>> gsMap;
                
        SecondPhase(const std::string&,
                    const std::map<std::shared_ptr<SlipSystem>,std::shared_ptr<GammaSurface>>&);
        
        
        double misfitEnergy(const Eigen::Matrix<double,dim,1>& s,const std::shared_ptr<SlipSystem>& ss) const;
        
    };

}
#endif
