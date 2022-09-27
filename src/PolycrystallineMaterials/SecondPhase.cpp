/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_SecondPhase_cpp_
#define model_SecondPhase_cpp_

#include <memory>
//#include <LatticePlaneBase.h>
//#include <LatticeVector.h>
//#include <RationalLatticeDirection.h>
#include <SecondPhase.h>
#include <TerminalColors.h>

namespace model
{
    
    template<int dim>
    SecondPhase<dim>::SecondPhase(const std::string& _name,
                                  const std::map<std::shared_ptr<SlipSystem>,std::shared_ptr<GammaSurface>>& _gsMap) :
    /* init */ name(_name)
    /* init */,gsMap(_gsMap)
    {
        
        std::cout<<greenBoldColor<<"Creating SecondPhase "<<name<<defaultColor<<std::endl;
        
    }

    template<int dim>
    double SecondPhase<dim>::misfitEnergy(const Eigen::Matrix<double,dim,1>& s,const std::shared_ptr<SlipSystem>& ss) const
    {
        
//        std::cout<<"inPlane="<<n<<std::endl;
//
//        std::cout<<"planes"<<std::endl;
//        for(const auto& pair : gsMap)
//        {
//            std::cout<<pair.first<<std::endl;
//        }
        
        const auto gammaIter(gsMap.find(ss));
        if(gammaIter!=gsMap.end())
        {
            return gammaIter->second->operator()(s);
        }
        else
        {
            return 0.0;
        }

    }

    
template struct SecondPhase<3>;
}
#endif
