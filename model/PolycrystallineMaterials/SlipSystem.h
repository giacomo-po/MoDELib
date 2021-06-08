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
#include <LatticeMath.h>
//#include <LatticePlaneBase.h>
//#include <LatticeVector.h>
//#include <RationalLatticeDirection.h>
#include <DislocationMobilityBase.h>
#include <PeriodicLatticeInterpolant.h>

namespace model
{
    
    struct GammaSurface : PeriodicLatticeInterpolant<2>
    {
        
        static constexpr int dim=3;
        static constexpr int lowerDim=dim-1;
        
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef Eigen::Matrix<size_t,lowerDim,1> VectorLowerDimI;
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;

        
        typedef Eigen::Matrix<double,lowerDim,lowerDim> MatrixLowerDim;
        
        /**********************************************************************/
        static MatrixDim getG2L(const VectorDim& x,
                                      const VectorDim& z)
        {
            const double xNorm(x.norm());
            const double zNorm(z.norm());
            assert(xNorm>FLT_EPSILON);
            assert(zNorm>FLT_EPSILON);
            assert(fabs(x.dot(z)<FLT_EPSILON*xNorm*zNorm));
            Eigen::Matrix3d temp(Eigen::Matrix3d::Identity());
            temp.col(2)=z/zNorm;
            temp.col(0)=x/xNorm;
            temp.col(1)=temp.col(2).cross(temp.col(0));
            return temp.transpose();
        }
        
        static MatrixLowerDim getLocalBasis(const LatticePlaneBase& n)
        {
            const Eigen::Matrix3d R(getG2L(n.primitiveVectors.first.cartesian(),n.cartesian().normalized()));
            MatrixLowerDim temp(MatrixLowerDim::Zero());
            temp.col(0)=(R*n.primitiveVectors.first.cartesian()).segment<2>(0);
            temp.col(1)=(R*n.primitiveVectors.second.cartesian()).segment<2>(0);
            return temp;
        }

        const LatticePlaneBase latticePlane;
        const Eigen::Matrix3d G2L;
        
//        GammaSurface(const LatticePlaneBase& n,
//                                   const VectorLowerDimI& nums_in,
//                                   const VectorLowerDimI& dens_in,
//                                   const Eigen::Matrix<double,Eigen::Dynamic,lowerDim+1>& f,
//                                   const Eigen::Matrix<double,Eigen::Dynamic,2*lowerDim+1>& df) :
//        /* init */ PeriodicLatticeInterpolant<2>(getLocalBasis(n),nums_in,dens_in,f,df)
//        /* init */,latticePlane(n)
//        /* init */,G2L(getG2L(n.primitiveVectors.first.cartesian(),n.cartesian().normalized()))
//        {
//
//            model::cout<<greenBoldColor<<"Creating GammaSurface on "<<n.cartesian().normalized().transpose()<<" plane"<<std::endl;
//
//
//        }
        
//        GammaSurface(const LatticePlaneBase& n,
//                     const Eigen::Matrix<double,Eigen::Dynamic,lowerDim>& waveVectors_in,
//                     const Eigen::Matrix<double,Eigen::Dynamic,lowerDim+1>& f,
//                     const Eigen::Matrix<double,Eigen::Dynamic,2*lowerDim+1>& df) :
//        /* init */ PeriodicLatticeInterpolant<2>(getLocalBasis(n),waveVectors_in,f,df)
//        /* init */,latticePlane(n)
//        /* init */,G2L(getG2L(n.primitiveVectors.first.cartesian(),n.cartesian().normalized()))
//        {
//
//            model::cout<<greenBoldColor<<"Creating GammaSurface on "<<n.cartesian().normalized().transpose()<<" plane"<<std::endl;
//
//
//        }

        GammaSurface(const LatticePlaneBase& n,
                     const Eigen::Matrix<double,Eigen::Dynamic,lowerDim>& waveVectors,
                     const Eigen::Matrix<double,Eigen::Dynamic,lowerDim+1>& f,
                     const int& rotSymm,
                     const std::vector<Eigen::Matrix<double,lowerDim,1>>& mirSymm) :
        /* init */ PeriodicLatticeInterpolant<2>(getLocalBasis(n),waveVectors,f,rotSymm,mirSymm)
        /* init */,latticePlane(n)
        /* init */,G2L(getG2L(n.primitiveVectors.first.cartesian(),n.cartesian().normalized()))
        {
            
            model::cout<<greenBoldColor<<"Creating GammaSurface on "<<n.cartesian().normalized().transpose()<<" plane"<<std::endl;
            
            
        }
        
        double operator()(const VectorDim& b)
        {
            const VectorDim bL(G2L*b);
            assert(std::fabs(bL(dim-1))<FLT_EPSILON && "SLIP VECTOR NOT ON GAMMA-SURFACE PLANE");
            return PeriodicLatticeInterpolant<2>::operator()(bL.segment<lowerDim>(0));
        }
        
    };
    
    struct SlipSystem : public StaticID<SlipSystem>
    {
        

        const LatticePlaneBase n;
//        const LatticeVector<3>  s;
        const RationalLatticeDirection<3>  s;
        const Eigen::Matrix<double,3,1>  unitNormal;
        const std::shared_ptr<DislocationMobilityBase> mobility;
        const std::shared_ptr<GammaSurface> gammaSurface;
        
//        SlipSystem(const LatticeVector<3>& a1,
//                   const LatticeVector<3>& a2,
//                   const LatticeVector<3>& slip_in,
//                   const std::shared_ptr<DislocationMobilityBase>& mobility_in):
//        /* init */ n(a1,a2)
//        /* init */,s(slip_in)
//        /* init */,unitNormal(n.cartesian().normalized())
//        /* init */,mobility(mobility_in)
//        {
//
//            model::cout<<greenColor<<"Creating SlipSystem "<<this->sID<<defaultColor<<std::endl;
//            model::cout<<"  s= "<<std::setprecision(15)<<std::scientific<<s.cartesian().transpose()<<std::endl;
//            model::cout<<"  n= "<<std::setprecision(15)<<std::scientific<<n.cartesian().transpose()<<std::endl;
//            model::cout<<"  mobility= "<<mobility->name<<std::endl;
//
//
//            if(n.dot(s)!=0)
//            {
//                model::cout<<"PLANE NORMAL AND SLIP DIRECTION ARE NOT ORTHOGONAL. EXITING."<<std::endl;
//                exit(EXIT_FAILURE);
//            }
//            if(!mobility)
//            {
//                model::cout<<"MOBILITY CANNOT BE A NULLPTR. EXITING."<<std::endl;
//                exit(EXIT_FAILURE);
//            }
//        }
        
        SlipSystem(const LatticeVector<3>& a1,
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
            
            model::cout<<greenColor<<"Creating SlipSystem "<<this->sID<<defaultColor<<std::endl;
            model::cout<<"  s= "<<std::setprecision(15)<<std::scientific<<s.cartesian().transpose()<<std::endl;
            model::cout<<"  n= "<<std::setprecision(15)<<std::scientific<<n.cartesian().transpose()<<std::endl;
            model::cout<<"  mobility= "<<mobility->name<<std::endl;
            
            
            if(s.dot(n)!=0)
            {
                model::cout<<"PLANE NORMAL AND SLIP DIRECTION ARE NOT ORTHOGONAL. EXITING."<<std::endl;
                exit(EXIT_FAILURE);
            }
            if(!mobility)
            {
                model::cout<<"MOBILITY CANNOT BE A NULLPTR. EXITING."<<std::endl;
                exit(EXIT_FAILURE);
            }
        }
        
        SlipSystem(const LatticeVector<3>& a1,
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
            
            model::cout<<greenColor<<"Creating partial SlipSystem "<<this->sID<<defaultColor<<std::endl;
            model::cout<<"  s= "<<std::setprecision(15)<<std::scientific<<s.cartesian().transpose()<<std::endl;
            model::cout<<"  n= "<<std::setprecision(15)<<std::scientific<<n.cartesian().transpose()<<std::endl;
            model::cout<<"  mobility= "<<mobility->name<<std::endl;
            
            
            if(s.dot(n)!=0)
            {
                model::cout<<"PLANE NORMAL AND SLIP DIRECTION ARE NOT ORTHOGONAL. EXITING."<<std::endl;
                exit(EXIT_FAILURE);
            }
            if(!mobility)
            {
                model::cout<<"MOBILITY CANNOT BE A NULLPTR. EXITING."<<std::endl;
                exit(EXIT_FAILURE);
            }
        }
        
        bool isPartial() const
        {
            return abs(s.rat.d)!=1;
        }

        
//        bool isSameAs(const LatticeVector<3>& s1,const ReciprocalLatticeDirection<3>& n1)
//        {
//            if(   ((s-s1).squaredNorm()==0 && (n-n1).squaredNorm()==0)
//               || ((s+s1).squaredNorm()==0 && (n+n1).squaredNorm()==0)
//               )
//            {
//                return true;
//            }
//            else
//            {
//                return false;
//            }
//        }
//
//        bool isSameAs(const RationalLatticeDirection<3>& s1,const ReciprocalLatticeDirection<3>& n1)
//        {
//            return isSameAs(s1.dir,n1);
//        }

        bool isSameAs(const RationalLatticeDirection<3>& s1,const ReciprocalLatticeDirection<3>& n1)
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
        
        double misfitEnergy(const Eigen::Matrix<double,3,1>& b)
        {
            return gammaSurface? gammaSurface->operator()(b) : 0.0;
        }
        
    };

}
#endif
