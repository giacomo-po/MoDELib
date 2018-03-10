/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationLinkingNumber_H_
#define model_DislocationLinkingNumber_H_

#include <iostream>
#include <deque>
#include <Eigen/Dense>
#include <cfloat>
#include <math.h>       /* asin */


namespace model
{
    /*! Class template that implements the calculation of pair-wise linking numbers
     *  between dislocation loops. Implementation from ref [1].
     *
     *  [1] Klenin, K., & Langowski, J. (2000). Computation of writhe in modeling of
     *  supercoiled DNA. Biopolymers, 54(5), 307-317.
     */
    template<typename DislocationNetworkType>
    class DislocationLinkingNumber
    {
        typedef typename DislocationNetworkType::LoopType LoopType;
        typedef typename DislocationNetworkType::LoopLinkType LoopLink;
        typedef typename DislocationNetworkType::VectorDim VectorDim;
        
        
        const DislocationNetworkType& DN;
        
        static constexpr double tol=FLT_EPSILON;
        
        /**********************************************************************/
        static double F(const double& t1,
                        const double& t2,
                        const double& cosBeta,
                        const double& sinBeta2,
                        const double& a0)
        {
            
            const double num=t1*t2+a0*a0*cosBeta;
            const double den=a0*sqrt(t1*t1+t2*t2-2.0*t1*t2*cosBeta+a0*a0*sinBeta2);
            //            return -atan2(den+tol,num)/(4.0*M_PI);
            //            return fabs(den)>tol? -atan2(den,num)/(4.0*M_PI) : -atan2(0.0,num)/(4.0*M_PI);
            return fabs(den)>tol? -atan2(den,num)/(4.0*M_PI) : 0.0;
            
        }
        
        /**********************************************************************/
        static VectorDim unitCross(const VectorDim& v1,
                                   const VectorDim& v2)
        {
            VectorDim temp(v1.cross(v2));
            const double tempNorm(temp.norm());
            return (tempNorm<tol)? VectorDim::Zero() : (temp/tempNorm).eval();
        }
        
        /**********************************************************************/
        static double myAsin(const double& x)
        {
            if(x>=1.0)
            {
                return 0.5*M_PI;
            }
            else if(x<=-1.0)
            {
                return -0.5*M_PI;
            }
            else
            {
                return asin(x);
            }
        }
        
        /**********************************************************************/
        static double segmentPairLN_a(const LoopLink& link1,
                                      const LoopLink& link2)
        {
            const VectorDim& R1(link1.source()->get_P());
            const VectorDim& R2(link1.  sink()->get_P());
            const VectorDim& R3(link2.source()->get_P());
            const VectorDim& R4(link2.  sink()->get_P());

            const VectorDim R13(R3-R1);
            const VectorDim R12(R2-R1);
            const VectorDim R34(R4-R3);

            const double R13norm(R13.norm());
            const double R12norm(R12.norm());
            const double R34norm(R34.norm());
            
            if(R13norm>tol && R12norm>tol && R34norm>tol)
            {
                const double prod(R34.cross(R12).dot(R13)/R13norm/R12norm/R34norm);
                if(fabs(prod)>tol)
                {
                    const VectorDim R14(R4-R1);
                    const VectorDim R24(R4-R2);
                    const VectorDim R23(R3-R2);
                    
                    const VectorDim n1(unitCross(R13,R14));
                    const VectorDim n2(unitCross(R14,R24));
                    const VectorDim n3(unitCross(R24,R23));
                    const VectorDim n4(unitCross(R23,R13));
                    
                    
                    const double Os= myAsin(n1.dot(n2))
                    /*            */+myAsin(n2.dot(n3))
                    /*            */+myAsin(n3.dot(n4))
                    /*            */+myAsin(n4.dot(n1));
                    
                    return Os/(4.0*M_PI)*prod/fabs(prod);
                    
                }
                else
                {
                    return 0.0;
                }
            }
            else
            {
                return 0.0;
            }
        }
        
        /**********************************************************************/
        static double segmentPairLN_b(const LoopLink& link1,
                                      const LoopLink& link2)
        {
            const VectorDim R1(link1.source()->get_P());
            const VectorDim R2(link2.source()->get_P());
            
            const VectorDim S1(link1.sink()->get_P()-link1.source()->get_P());
            const VectorDim S2(link2.sink()->get_P()-link2.source()->get_P());
            const double s1(S1.norm());
            const double s2(S2.norm());
            const VectorDim e1(S1.normalized());
            const VectorDim e2(S2.normalized());
            const VectorDim R12(R2-R1);
            
            const double cosBeta(e1.dot(e2));
            const double sinBeta2(1.0-cosBeta*cosBeta);
            const double a1=R12.dot(e2*cosBeta-e1)/(sinBeta2+tol);
            const double a2=R12.dot(e2-e1*cosBeta)/(sinBeta2+tol);
            const double a0=R12.dot(e1.cross(e2))/(sinBeta2+tol);
            
            //                    std::cout<<"R1="<<std::setprecision(15)<<std::scientific<<R1.transpose()<<std::endl;
            //                    std::cout<<"R2="<<std::setprecision(15)<<std::scientific<<R2.transpose()<<std::endl;
            //                    std::cout<<"S1="<<std::setprecision(15)<<std::scientific<<S1.transpose()<<std::endl;
            //                    std::cout<<"S2="<<std::setprecision(15)<<std::scientific<<S2.transpose()<<std::endl;
            std::cout<<"R1="<<std::setprecision(15)<<R1.transpose()<<std::endl;
            std::cout<<"R2="<<std::setprecision(15)<<R2.transpose()<<std::endl;
            std::cout<<"S1="<<std::setprecision(15)<<S1.transpose()<<std::endl;
            std::cout<<"S2="<<std::setprecision(15)<<S2.transpose()<<std::endl;
            const double temp= F(a1+s1,a2+s2,cosBeta,sinBeta2,a0)
            /*             */ -F(a1+s1,a2   ,cosBeta,sinBeta2,a0)
            /*             */ -F(a1   ,a2+s2,cosBeta,sinBeta2,a0)
            /*             */ +F(a1   ,a2   ,cosBeta,sinBeta2,a0);
            std::cout<<"tamp="<<std::setprecision(15)<<std::scientific<<temp<<std::endl;
            
            return F(a1+s1,a2+s2,cosBeta,sinBeta2,a0)
            /* */ -F(a1+s1,a2   ,cosBeta,sinBeta2,a0)
            /* */ -F(a1   ,a2+s2,cosBeta,sinBeta2,a0)
            /* */ +F(a1   ,a2   ,cosBeta,sinBeta2,a0);
        }
        
        
        /**********************************************************************/
        static double linkingNumber(const LoopType& loop1,
                                    const LoopType& loop2)
        {
            double Ln=0.0;
            
            for(const auto& link1 : loop1.links())
            {
                for(const auto& link2 : loop2.links())
                {
                    Ln +=  segmentPairLN_a(*link1.second,*link2.second);
                }
            }
            return Ln;
        }
        
    public:
        /**********************************************************************/
        DislocationLinkingNumber(const DislocationNetworkType& DN_in) :
        /* init */ DN(DN_in)
        {
            
            
        }
        
        
        /**********************************************************************/
        template <class T>
        friend T& operator << (T& os, const DislocationLinkingNumber<DislocationNetworkType>& ln)
        {
            
            for(auto loopIter1=ln.DN.loops().begin();loopIter1!=ln.DN.loops().end();++loopIter1)
            {
                for(auto loopIter2=loopIter1;loopIter2!=ln.DN.loops().end();++loopIter2)
                {
                    if(loopIter1!=loopIter2)
                    {
                        os  <<loopIter1->first<<" "<<loopIter2->first<<" ";
                        os  <<linkingNumber(*loopIter1->second,*loopIter2->second)<<" ";
                        os  <<loopIter1->second->flow().cartesian().transpose()<<" ";
                        os  <<loopIter1->second->glidePlane.n.cartesian().transpose()<<" ";
                        os  <<loopIter2->second->flow().cartesian().transpose()<<" ";
                        os  <<loopIter2->second->glidePlane.n.cartesian().transpose()<<" ";
                        os  <<"\n";
                    }
                }
            }
            
            return os;
        }
        
    };
    
}
#endif

