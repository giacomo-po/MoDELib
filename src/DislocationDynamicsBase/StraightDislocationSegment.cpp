/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2018 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2018 by Yinan Cui <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_StraightDislocationSegment_cpp_
#define model_StraightDislocationSegment_cpp_

#ifndef _MODEL_NON_SINGULAR_DD_
#define _MODEL_NON_SINGULAR_DD_ 1
#endif

#include <cfloat>
#include <StraightDislocationSegment.h>

namespace model
{
    
        template <int dim>
        typename StraightDislocationSegment<dim>::MatrixDim StraightDislocationSegment<dim>::nonSymmStress_kernel(const VectorDim& r) const
        {
#if _MODEL_NON_SINGULAR_DD_ == 0
            /* Devincre, B., & Condat, M. (1992). Model validation of a 3D
             * simulation of dislocation dynamics: Discretization and line tension
             * effects. Acta Metallurgica Et Materialia, 40(10), 2629–2637.
             */
            
            static_assert(0,"THERE IS NO CONSISTENT REGULARIZATION OF THIS EXPRESSION. USE _MODEL_NON_SINGULAR_DD_=1");
            
//            const Scalar L(r.dot(t));
//            const Scalar R(r.norm()+DislocationFieldBase<dim>::a);
//            const VectorDim rho(r-L*t); // distance vector to the line
//            const VectorDim Y((L+R)*t+rho); // = R*t + r
//            // Y2=R^2+R^2+2R*(r*t)=2R*(R+L)
////            const Scalar Y2(Y.squaredNorm()+DislocationFieldBase<dim>::a2);
////            return (material.C1* b.cross(Y)*t.transpose()
////            /*                 */ -b.cross(t)*Y.transpose()
////            /*                 */ -b.dot(Y.cross(t))*(2.0/Y2*rho*Y.transpose()+0.5*(MatrixDim::Identity()+t*t.transpose()+2.0/Y2*L/R*Y*Y.transpose()))
////                    )*2.0/Y2;
//            
//            
////            const Scalar f1(2.0/Y.squaredNorm());
//            const Scalar f1(1.0/(R*(R+L))); // =2/Y2=2/(2R*(R+L))
//
////            const Scalar f1(2.0/(Y.squaredNorm()+DislocationFieldBase<dim>::a2));
////            const Scalar f1(2.0/Y2);
//            const Scalar f2(material.C1*f1);
//            const Scalar bDotYcrosst(b.dot(Y.cross(t)));
//            const Scalar f3(bDotYcrosst*f1*f1);
//            const Scalar f4(0.5*f1*bDotYcrosst);
//            const Scalar f5(f4*f1*L/R);
//            //            const Scalar f5(0.5*f3*L/R);
//            
//            return  f2*b.cross(Y)*t.transpose()
//            /*  */ -f1*b.cross(t)*Y.transpose()
//            /*  */ -f3*rho*Y.transpose()
//            /*  */ -f4*(MatrixDim::Identity()+t*t.transpose())
//            /*  */ -f5*Y*Y.transpose();

//            const Scalar R(r.norm());
            const Scalar Ra(sqrt(r.squaredNorm()+DislocationFieldBase<dim>::a2));
            const VectorDim Y(r+Ra*t);
            const Scalar Yt(Y.dot(t));
//            const Scalar Yta(Yt+DislocationFieldBase<dim>::a);
            const Scalar Yta(sqrt(Yt*Yt+DislocationFieldBase<dim>::a2));
            const Scalar Y2(Y.squaredNorm()+DislocationFieldBase<dim>::a2);
            const Scalar bYt(b.cross(Y).dot(t));

            
            const Scalar f1(2.0/Y2);
            
            
            return  f1*material.C1*t*(b.cross(Y)).transpose()
            /*   */-f1*Y*bCt.transpose()
            /*   */-f1*bYt/Yta*t*r.transpose()
            /*   */-0.5*f1*bYt*MatrixDim::Identity()
            /*   */-f1*bYt*(R+Yt)/Ra/Y2*r*r.transpose()
            /*   */-0.5*f1*bYt*R/Yta*t*t.transpose();
            
            
#elif _MODEL_NON_SINGULAR_DD_ == 1 /* Cai's non-singular theory */

//            const double Ra2=r.squaredNorm()+DislocationFieldBase<dim>::a2;
//            const double Ra=sqrt(Ra2);
//            const double Ra3=std::pow(Ra,3);
//            const double rdt=r.dot(t);
//            const double rdt2=std::pow(rdt,2);
//            const double A1=-(rdt*(3.0*Ra2-rdt2))/(std::pow(Ra2-rdt2,2)*Ra3);
//            const double A2=1.0/Ra3-rdt*A1;
//            const double A6=-rdt/((Ra2-rdt2)*Ra);
//            const double A3=-rdt/Ra3+A6+rdt2*A1;
//            const double A4=A6+DislocationFieldBase<dim>::a2*A1;
//            const double A5=-material.C1*A6-0.5*DislocationFieldBase<dim>::a2*material.C1*A1;
//            const double A7=material.nu/Ra-rdt*A6-0.5*DislocationFieldBase<dim>::a2*material.C1*A2;
//            
//            const double rbt(r.cross(b).dot(t));
//            
//            return  0.5*rbt*A1*r*r.transpose()
//            /*  */ +    rbt*A2*t*r.transpose()
//            /*  */ +0.5*rbt*A3*t*t.transpose()
//            /*  */ +0.5*rbt*A4*MatrixDim::Identity()
//            /*  */ +A5*r.cross(b)*t.transpose()
//            /*  */ +A6*t.cross(b)*r.transpose()
//            /*  */ +A7*t.cross(b)*t.transpose();

            const Scalar Ra2=r.squaredNorm()+DislocationFieldBase<dim>::a2;
            const Scalar Ra(sqrt(Ra2));
            const VectorDim Ya(r+Ra*t);
            const Scalar Yat(Ya.dot(t));
            const Scalar Ya2a2(Ya.squaredNorm()+DislocationFieldBase<dim>::a2);
            const VectorDim bYa(b.cross(Ya));
            const Scalar bYat(bYa.dot(t));
            
            
            const Scalar f1(2.0/Ya2a2);
            
            return f1*material.C1*(1.0+DislocationFieldBase<dim>::a2/Ya2a2)*t*bYa.transpose()
            /*  */+f1*material.C1*0.5*DislocationFieldBase<dim>::a2/Ra2*t*b.cross(r).transpose()
            /*  */-f1*Ya*bCt.transpose()
            /*  */-f1*bYat/Yat*t*r.transpose()
            /*  */-0.5*f1*bYat*(1.0+2.0*DislocationFieldBase<dim>::a2/Ya2a2+DislocationFieldBase<dim>::a2/Ra2)*MatrixDim::Identity()
            /*  */-f1*bYat*(Ra+Yat)/Ra/Ya2a2*r*r.transpose()
            /*  */-0.5*f1*bYat*Ra/Yat*t*t.transpose();
            
#elif _MODEL_NON_SINGULAR_DD_ == 2 /* Lazar's non-singular theory */
            static_assert(0,"NOT IMPLEMENTED. USE _MODEL_NON_SINGULAR_DD_=1");
#else
#error Unsupported choice of field regularization
#endif
        }
        
        template <int dim>
        typename StraightDislocationSegment<dim>::VectorDim StraightDislocationSegment<dim>::displacement_kernel(const VectorDim& r) const
        {
            const Scalar Ra(sqrt(r.squaredNorm()+DislocationFieldBase<dim>::a2));
            const VectorDim Ya(r+Ra*t);
            const Scalar Yat(Ya.dot(t));
            return -(2.0-0.5/material.C1+(2.0-1.0/material.C1)*log(Yat)-DislocationFieldBase<dim>::a2/Ra/Yat)/8.0/M_PI*bCt
            /*  */ +bCt.dot(r)/Yat/Ra/material.C1/8.0/M_PI*Ya;
        }
        
        template <int dim>
        double StraightDislocationSegment<dim>::elasticInteractionEnergy_kernel(const VectorDim& z,const VectorDim& tA,const VectorDim& bA) const
        {
            const Scalar za(sqrt(z.squaredNorm()+DislocationFieldBase<dim>::a2));
            const VectorDim Ya(z+za*t);
            const Scalar Yat(Ya.dot(t));
            const Scalar logYat(log(Yat));
            const Scalar zaYat(za*Yat);

            const Scalar bAtA(bA.dot(tA));
            const Scalar bt(b.dot(t));
            const Scalar bAt(bA.dot(t));
            const Scalar btA(b.dot(tA));
            const Scalar bAb(bA.dot(b));
            const Scalar tAt(tA.dot(t));
            const Scalar zt(z.dot(t));
            const Scalar zb(z.dot(b));
            const Scalar zbA(z.dot(bA));
            return (0.5*material.C1*bAtA*bt+material.nu*bAt*btA)*(2.0+2.0*logYat-DislocationFieldBase<dim>::a2/zaYat)
            /*  */-bAb*tAt*(1.5+logYat-DislocationFieldBase<dim>::a2/zaYat)
            /*  */-bAt*bt*tAt*(0.5+logYat+zt/Yat)
            /*  */+tAt/Yat*(zbA*zb/za+bAt*zb+bt*zbA);
        }
                
        template <int dim>
        StraightDislocationSegment<dim>::StraightDislocationSegment(const PolycrystallineMaterialBase& mat,
                                   const VectorDim& _P0,
                                   const VectorDim& _P1,
                                   const VectorDim& _b,
                                   const double& _length,
                                   const VectorDim& _t) :
        /* init list */ material(mat),
        /* init list */ P0(_P0),
        /* init list */ P1(_P1),
        /* init list */ b(_b),
        /* init list */ length(_length),
        /* init list */ t(_t),
        /* init list */ bCt(b.cross(t))
        {/*!\param[in] _P0 starting point of the segment
          * \param[in] _P0 ending point of the segment
          * \param[in] _b Burgers vector of the segment
          */
        
        }

        template <int dim>
        void StraightDislocationSegment<dim>::updateGeometry()
        {
            bCt=b.cross(t);
        }

        template <int dim>
        typename StraightDislocationSegment<dim>::MatrixDim StraightDislocationSegment<dim>::nonSymmStress(const VectorDim& x) const
        {
            return length>FLT_EPSILON? nonSymmStress_kernel(P1-x)-nonSymmStress_kernel(P0-x) : MatrixDim::Zero().eval();
        }
        
        template <int dim>
        typename StraightDislocationSegment<dim>::MatrixDim StraightDislocationSegment<dim>::stress(const VectorDim& x) const
        {
            const MatrixDim temp = nonSymmStress(x);
            return material.C2*(temp+temp.transpose());
//            const Eigen::Matrix<float,dim,dim> tempF((material.C2*(temp+temp.transpose())).template cast<float>());
//            return tempF.template cast<double>();
        }
        
        template <int dim>
        typename StraightDislocationSegment<dim>::VectorDim StraightDislocationSegment<dim>::displacement(const VectorDim& x) const
        {/*!\returns the line-integral part of the displacement contribution of this straight segment.
          * Note: the return value  does NOT include the solid-angle contribution to the displacement field.
          */
            return length>FLT_EPSILON? (displacement_kernel(P1-x)-displacement_kernel(P0-x)).eval() : VectorDim::Zero();
        }
        
        template <int dim>
        double StraightDislocationSegment<dim>::elasticInteractionEnergy(const VectorDim& xA,const VectorDim& tA,const VectorDim& bA) const
        {
            
            return -0.5*material.C2*(elasticInteractionEnergy_kernel(P1-xA,tA,bA)-elasticInteractionEnergy_kernel(P0-xA,tA,bA));
        }
        

template class StraightDislocationSegment<3>;
	
}
#endif

