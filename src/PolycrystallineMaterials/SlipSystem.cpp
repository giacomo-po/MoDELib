/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_SLIPSYSTEM_cpp_
#define model_SLIPSYSTEM_cpp_

#include <memory>
#include <assert.h>
#include <LatticeModule.h>
//#include <LatticePlaneBase.h>
//#include <LatticeVector.h>
//#include <RationalLatticeDirection.h>
#include <DislocationMobilityBase.h>
#include <SlipSystem.h>

namespace model
{

    SlipSystem::SlipSystem(const LatticeVector<3>& a1,
                           const LatticeVector<3>& a2,
                           const LatticeVector<3>& slip_in,
                           const std::shared_ptr<DislocationMobilityBase>& mobility_in,
                           const std::shared_ptr<GammaSurface>& gammaSurface_in):
    /* init */ n(a1,a2)
    /* init */,s(slip_in)
    /* init */,unitNormal(n.cartesian().normalized())
    /* init */,mobility(mobility_in)
    /* init */,gammaSurface(gammaSurface_in)
    {
        
        std::cout<<greenColor<<"Creating SlipSystem "<<this->sID<<defaultColor<<std::endl;
        std::cout<<"  s= "<<std::setprecision(15)<<std::scientific<<s.cartesian().transpose()<<std::endl;
        std::cout<<"  n= "<<std::setprecision(15)<<std::scientific<<n.cartesian().transpose()<<std::endl;
        std::cout<<"  mobility= "<<mobility->name<<std::endl;
        
        
        if(s.dot(n)!=0)
        {
            std::cout<<"PLANE NORMAL AND SLIP DIRECTION ARE NOT ORTHOGONAL. EXITING."<<std::endl;
            exit(EXIT_FAILURE);
        }
        if(!mobility)
        {
            std::cout<<"MOBILITY CANNOT BE A NULLPTR. EXITING."<<std::endl;
            exit(EXIT_FAILURE);
        }
    }

    SlipSystem::SlipSystem(const LatticeVector<3>& a1,
                           const LatticeVector<3>& a2,
                           const RationalLatticeDirection<3>& slip_in,
                           const std::shared_ptr<DislocationMobilityBase>& mobility_in,
                           const std::shared_ptr<GammaSurface>& gammaSurface_in):
    /* init */ n(a1,a2)
    /* init */,s(slip_in)
    /* init */,unitNormal(n.cartesian().normalized())
    /* init */,mobility(mobility_in)
    /* init */,gammaSurface(gammaSurface_in)
    {
        
        std::cout<<greenColor<<"Creating partial SlipSystem "<<this->sID<<defaultColor<<std::endl;
        std::cout<<"  s= "<<std::setprecision(15)<<std::scientific<<s.cartesian().transpose()<<std::endl;
        std::cout<<"  n= "<<std::setprecision(15)<<std::scientific<<n.cartesian().transpose()<<std::endl;
        std::cout<<"  mobility= "<<mobility->name<<std::endl;
        
        
        if(s.dot(n)!=0)
        {
            std::cout<<"PLANE NORMAL AND SLIP DIRECTION ARE NOT ORTHOGONAL. EXITING."<<std::endl;
            exit(EXIT_FAILURE);
        }
        if(!mobility)
        {
            std::cout<<"MOBILITY CANNOT BE A NULLPTR. EXITING."<<std::endl;
            exit(EXIT_FAILURE);
        }
    }

    bool SlipSystem::isPartial() const
    {
        return abs(s.rat.d)!=1;
    }

    bool SlipSystem::isSameAs(const RationalLatticeDirection<3>& s1,const ReciprocalLatticeDirection<3>& n1)
    {
        if(   ((s-s1).squaredNorm()==0 && (n-n1).squaredNorm()==0)
           || ((s+s1).squaredNorm()==0 && (n+n1).squaredNorm()==0)
           )
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    double SlipSystem::misfitEnergy(const Eigen::Matrix<double,3,1>& b)
    {
        return gammaSurface? gammaSurface->operator()(b) : 0.0;
    }

}
#endif
