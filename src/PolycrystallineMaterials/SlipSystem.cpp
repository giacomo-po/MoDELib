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

    SlipSystem::SlipSystem(const LatticePlaneBase& n_in,
                           const LatticeVector<3>& slip_in,
                           const std::shared_ptr<DislocationMobilityBase>& mobility_in,
                           const std::shared_ptr<GammaSurface>& gammaSurface_in,
                           const std::shared_ptr<GlidePlaneNoise>& planeNoise_in):
    /* init */ n(n_in)
    /* init */,s(slip_in)
    /* init */,unitNormal(n.cartesian().normalized())
    /* init */,unitSlip(s.cartesian().normalized())
    /* init */,unitSlipFull(n.primitiveVectors.first.cartesian().normalized())
    /* init */,G2Lfull((MatrixDim()<<unitSlipFull,unitNormal.cross(unitSlipFull),unitNormal).finished().transpose())
    /* init */,mobility(mobility_in)
    /* init */,gammaSurface(gammaSurface_in)
    /* init */,planeNoise(planeNoise_in)
    {
        
        std::cout<<greenColor<<"Creating SlipSystem "<<this->sID<<defaultColor<<std::endl;
        std::cout<<"  s= "<<std::setprecision(15)<<std::scientific<<s.cartesian().transpose()<<std::endl;
        std::cout<<"  n= "<<std::setprecision(15)<<std::scientific<<n.cartesian().transpose()<<std::endl;
        std::cout<<"  mobility= "<<mobility->name<<std::endl;
        
        
        if(s.dot(n)!=0)
        {
            throw std::runtime_error("SlipSystem: PLANE NORMAL AND SLIP DIRECTION ARE NOT ORTHOGONAL.");
        }
        if(!mobility)
        {
            throw std::runtime_error("SlipSystem: MOBILITY CANNOT BE A NULLPTR.");
        }
        
        if((G2Lfull*G2Lfull.transpose()-MatrixDim::Identity()).squaredNorm()>FLT_EPSILON)
        {
            throw std::runtime_error("G2Lfull is not orthogonal.");
        }
    }

    SlipSystem::SlipSystem(const LatticePlaneBase& n_in,
                           const RationalLatticeDirection<3>& slip_in,
                           const std::shared_ptr<DislocationMobilityBase>& mobility_in,
                           const std::shared_ptr<GammaSurface>& gammaSurface_in,
                           const std::shared_ptr<GlidePlaneNoise>& planeNoise_in):
    /* init */ n(n_in)
    /* init */,s(slip_in)
    /* init */,unitNormal(n.cartesian().normalized())
    /* init */,unitSlip(s.cartesian().normalized())
    /* init */,unitSlipFull(n.primitiveVectors.first.cartesian().normalized())
    /* init */,G2Lfull((MatrixDim()<<unitSlipFull,unitNormal.cross(unitSlipFull),unitNormal).finished())
    /* init */,mobility(mobility_in)
    /* init */,gammaSurface(gammaSurface_in)
    /* init */,planeNoise(planeNoise_in)
    {
        
        std::cout<<greenColor<<"Creating partial SlipSystem "<<this->sID<<defaultColor<<std::endl;
        std::cout<<"  s= "<<std::setprecision(15)<<std::scientific<<s.cartesian().transpose()<<std::endl;
        std::cout<<"  n= "<<std::setprecision(15)<<std::scientific<<n.cartesian().transpose()<<std::endl;
        std::cout<<"  mobility= "<<mobility->name<<std::endl;
        
        
        if(s.dot(n)!=0)
        {
            throw std::runtime_error("SlipSystem: PLANE NORMAL AND SLIP DIRECTION ARE NOT ORTHOGONAL.");
        }
        if(!mobility)
        {
            throw std::runtime_error("SlipSystem: MOBILITY CANNOT BE A NULLPTR.");
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

    Eigen::Matrix<double,2,1> SlipSystem::globalToLocal(const VectorDim& x) const
    {
        return (G2Lfull*x).template segment<2>(0);
    }

    Eigen::Matrix<double,3,1> SlipSystem::localToGlobal(const Eigen::Matrix<double,2,1>& x) const
    {
        return (G2Lfull.transpose()).template block<3,2>(0,0)*x;
    }

    std::tuple<Eigen::Matrix<double,3,3>,double,double> SlipSystem::gridInterp(const VectorDim& x) const
    {   // Added by Hyunsoo (hyunsol@g.clemson.edu)
        
        if(planeNoise)
        {
            const std::tuple<double,double,double> gridNoise(planeNoise->gridInterp(globalToLocal(x)));
            
            const Eigen::Matrix<double,3,3> ssStress( std::get<0>(gridNoise)*(G2Lfull.row(2).transpose()*G2Lfull.row(0)+G2Lfull.row(0).transpose()*G2Lfull.row(2))
                                                     +std::get<1>(gridNoise)*(G2Lfull.row(2).transpose()*G2Lfull.row(1)+G2Lfull.row(1).transpose()*G2Lfull.row(2)));
            const double ssRSS((ssStress*unitSlip).dot(unitNormal));
            return std::make_tuple(ssStress,ssRSS,std::get<2>(gridNoise));
        }
        else
        {
            return std::make_tuple(Eigen::Matrix<double,3,3>::Zero(),0.0,0.0);
        }
    }

std::tuple<Eigen::Matrix<double,3,3>,double,double> SlipSystem::gridVal(const Eigen::Array<int,2,1>& idx) const
{   // Added by Hyunsoo (hyunsol@g.clemson.edu)
    
    if(planeNoise)
    {
        const std::tuple<double,double,double> gridNoise(planeNoise->gridVal(idx));
        
        const Eigen::Matrix<double,3,3> ssStress( std::get<0>(gridNoise)*(G2Lfull.row(2).transpose()*G2Lfull.row(0)+G2Lfull.row(0).transpose()*G2Lfull.row(2))
                                                 +std::get<1>(gridNoise)*(G2Lfull.row(2).transpose()*G2Lfull.row(1)+G2Lfull.row(1).transpose()*G2Lfull.row(2)));
        const double ssRSS((ssStress*unitSlip).dot(unitNormal));
        return std::make_tuple(ssStress,ssRSS,std::get<2>(gridNoise));
    }
    else
    {
        return std::make_tuple(Eigen::Matrix<double,3,3>::Zero(),0.0,0.0);
    }
}

}
#endif
