/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_Plane_CPP_
#define model_Plane_CPP_

#include <cfloat>
#include <tuple>
//#include <map>
#include <Eigen/Dense>


#include <Plane.h>


namespace model
{
    
        /**********************************************************************/
        template <int dim>
        Plane<dim>::Plane(const VectorDim& p,const VectorDim& n) :
        /* init */ P(p)
        /* init */,unitNormal(n.normalized())
        /* init */,L2G(getL2G(unitNormal))
        {
            if(n.norm()<FLT_EPSILON)
            {
                std::cout<<"Plane: n="<<n.transpose()<<std::endl;
                throw std::runtime_error("Plane must have non-zero normal");
            }
        }
        
        /**********************************************************************/
        template <int dim>
        bool Plane<dim>::contains(const VectorDim& P0) const
        {
//            const double PP0((P-P0).norm());
//            return PP0<FLT_EPSILON? true : (fabs((P0-P).dot(unitNormal))<FLT_EPSILON*PP0);
            return fabs((P0-P).dot(unitNormal))<FLT_EPSILON;
        }

template <int dim>
bool Plane<dim>::isAbove(const VectorDim& P0) const
{//!\returns true if the plane is above point P0.
    return (P0-P).dot(unitNormal)<0.0;
}

template <int dim>
bool Plane<dim>::isBelow(const VectorDim& P0) const
{//!\returns true if the plane is below point P0.
    return (P0-P).dot(unitNormal)>0.0;
}
        
        /**********************************************************************/
        template <int dim>
        typename Plane<dim>::VectorDim Plane<dim>::snapToPlane(const VectorDim& P0) const
        {
            return P0-(P0-P).dot(unitNormal)*unitNormal;
        }
        
        /**********************************************************************/
        template <int dim>
        double Plane<dim>::distanceTo(const VectorDim& P0) const
        {
            return fabs((P0-P).dot(unitNormal));
        }
                
        /**********************************************************************/
        template <int dim>
        typename Plane<dim>::VectorLowerDim Plane<dim>::localPosition(const VectorDim& point) const
        {
            const VectorDim pointLocal(L2G.transpose()*(point-P));
            if(fabs(pointLocal(2))>FLT_EPSILON)
            {
                std::cout<<"point="<<point.transpose()<<std::endl;
                std::cout<<"P="<<P.transpose()<<std::endl;
                std::cout<<"L2G=\n"<<L2G<<std::endl;
                std::cout<<"pointLocal="<<pointLocal.transpose()<<std::endl;
                throw std::runtime_error("Local point has non-zero z-coordinate");
            }
            return pointLocal.template segment<2>(0);
        }

        /**********************************************************************/
        template <int dim>
        typename Plane<dim>::VectorDim Plane<dim>::globalPosition(const VectorLowerDim& point) const
        {// terurns the position on the plane in global goordinates
            return L2G.template block<dim,dim-1>(0,0)*point+P;
        }

        /**********************************************************************/
        template <int dim>
        typename Plane<dim>::MatrixDim Plane<dim>::getL2G(VectorDim z)
        {
            //            const double xNorm(x.norm());
            const double zNorm(z.norm());
            assert(zNorm>FLT_EPSILON);
            z/=zNorm;
            
            VectorDim x(VectorDim::UnitX().cross(z));
            double xNorm(x.norm());
            if(xNorm>FLT_EPSILON)
            {
                x=x/xNorm;
            }
            else
            {
                x=VectorDim::UnitY().cross(z);
                xNorm=x.norm();
                if(xNorm>FLT_EPSILON)
                {
                    x=x/xNorm;
                }
                else
                {
                    x=VectorDim::UnitZ().cross(z);
                    xNorm=x.norm();
                    if(xNorm>FLT_EPSILON)
                    {
                        x=x/xNorm;
                    }
                    else
                    {
                        throw std::runtime_error("CANNOT FIND VECTOR ORTHOGONAL TO z");
                    }
                }
            }
            
            assert(std::fabs(x.norm()-1.0)<FLT_EPSILON);
            assert(std::fabs(z.norm()-1.0)<FLT_EPSILON);
            assert(fabs(x.dot(z)<FLT_EPSILON));
            MatrixDim temp(Eigen::Matrix3d::Identity());
            temp.col(2)=z;
            temp.col(0)=x;
            temp.col(1)=temp.col(2).cross(temp.col(0));
            return temp;
        }

    template struct Plane<3>;

}
#endif
