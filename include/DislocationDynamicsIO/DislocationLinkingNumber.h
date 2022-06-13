/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationLinkingNumber_H_
#define model_DislocationLinkingNumber_H_

#include <numbers>
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <cfloat>
#include <math.h>       /* asin */

#include <SegmentSegmentDistance.h>

namespace model
{
    
    
    template<int dim>
    struct LinkingNumber
    {
        
        static constexpr double tol=FLT_EPSILON;
        
        
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        
        //        /**********************************************************************/
        //        static double F(const double& t1,
        //                        const double& t2,
        //                        const double& cosBeta,
        //                        const double& sinBeta2,
        //                        const double& a0)
        //        {
        //
        //            const double num=t1*t2+a0*a0*cosBeta;
        //            const double den=a0*sqrt(t1*t1+t2*t2-2.0*t1*t2*cosBeta+a0*a0*sinBeta2);
        //            //            return -atan2(den+tol,num)/(4.0*std::numbers::pi);
        //            //            return fabs(den)>tol? -atan2(den,num)/(4.0*std::numbers::pi) : -atan2(0.0,num)/(4.0*std::numbers::pi);
        //            return fabs(den)>tol? -atan2(den,num)/(4.0*std::numbers::pi) : 0.0;
        //
        //        }
        
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
                return 0.5*std::numbers::pi;
            }
            else if(x<=-1.0)
            {
                return -0.5*std::numbers::pi;
            }
            else
            {
                return asin(x);
            }
        }
        
        /**********************************************************************/
        static double segmentPairLN(const VectorDim& R1,
                                    const VectorDim& R2,
                                    const VectorDim& R3,
                                    const VectorDim& R4)
        {
            //            const VectorDim& R1(link1.source->get_P());
            //            const VectorDim& R2(link1.  sink->get_P());
            //            const VectorDim& R3(link2.source->get_P());
            //            const VectorDim& R4(link2.  sink->get_P());
            
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
                    
                    return Os/(4.0*std::numbers::pi)*prod/fabs(prod);
                    
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
        static double loopPairLN(const std::vector<VectorDim>& loop1,
                                 const std::vector<VectorDim>& loop2)
        {
            
            double temp(0.0);
            
            for(size_t i=0;i<loop1.size();++i)
            {
                
                const size_t i1(i==loop1.size()-1? 0 : i+1);
                
                for(size_t j=0;j<loop2.size();++j)
                {
                    
                    const size_t j1(j==loop2.size()-1? 0 : j+1);
                    
                    temp+=segmentPairLN(loop1[i],
                                        loop1[i1],
                                        loop2[j],
                                        loop2[j1]);
                }
            }
            return temp;
        }
        
        /**********************************************************************/
        static double loopPairHelicity(const std::vector<VectorDim>& loop1,
                                       const VectorDim& b1,
                                       const std::vector<VectorDim>& loop2,
                                       const VectorDim& b2)
        {
            
            
            return loopPairLN(loop1,loop2)*b1.dot(b2);
            
        }
        
    };
    
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
        static constexpr int dim=TypeTraits<DislocationNetworkType>::dim;
        
        
        
        /**********************************************************************/
        static double segmentPairLN(const LoopLink& link1,
                                    const LoopLink& link2)
        {
            
            const VectorDim& R1(link1.source->get_P());
            const VectorDim& R2(link1.  sink->get_P());
            const VectorDim& R3(link2.source->get_P());
            const VectorDim& R4(link2.  sink->get_P());
            
            if(   (link1.source->sID==link2.source->sID && link1.sink->sID==link2.sink->sID)
               || (link1.source->sID==link2.sink->sID && link1.sink->sID==link2.source->sID))
            {
                return 0.0;
            }
            else
            {
                return LinkingNumber<dim>::segmentPairLN(R1,R2,R3,R4);
            }
        }
        
        
        
        
        
        //        /**********************************************************************/
        //        static double segmentPairLN(const VectorDim& R1,
        //                                      const VectorDim& R2,
        //                                      const VectorDim& R3,
        //                                      const VectorDim& R4)
        //        {
        ////            const VectorDim& R1(link1.source->get_P());
        ////            const VectorDim& R2(link1.  sink->get_P());
        ////            const VectorDim& R3(link2.source->get_P());
        ////            const VectorDim& R4(link2.  sink->get_P());
        //
        //            const VectorDim R13(R3-R1);
        //            const VectorDim R12(R2-R1);
        //            const VectorDim R34(R4-R3);
        //
        //            const double R13norm(R13.norm());
        //            const double R12norm(R12.norm());
        //            const double R34norm(R34.norm());
        //
        //            if(R13norm>tol && R12norm>tol && R34norm>tol)
        //            {
        //                const double prod(R34.cross(R12).dot(R13)/R13norm/R12norm/R34norm);
        //                if(fabs(prod)>tol)
        //                {
        //                    const VectorDim R14(R4-R1);
        //                    const VectorDim R24(R4-R2);
        //                    const VectorDim R23(R3-R2);
        //
        //                    const VectorDim n1(unitCross(R13,R14));
        //                    const VectorDim n2(unitCross(R14,R24));
        //                    const VectorDim n3(unitCross(R24,R23));
        //                    const VectorDim n4(unitCross(R23,R13));
        //
        //
        //                    const double Os= myAsin(n1.dot(n2))
        //                    /*            */+myAsin(n2.dot(n3))
        //                    /*            */+myAsin(n3.dot(n4))
        //                    /*            */+myAsin(n4.dot(n1));
        //
        //                    return Os/(4.0*std::numbers::pi)*prod/fabs(prod);
        //
        //                }
        //                else
        //                {
        //                    return 0.0;
        //                }
        //            }
        //            else
        //            {
        //                return 0.0;
        //            }
        //        }
        
        //        /**********************************************************************/
        //        static double segmentPairLN_b(const LoopLink& link1,
        //                                      const LoopLink& link2)
        //        {
        //            const VectorDim R1(link1.source->get_P());
        //            const VectorDim R2(link2.source->get_P());
        //
        //            const VectorDim S1(link1.sink->get_P()-link1.source->get_P());
        //            const VectorDim S2(link2.sink->get_P()-link2.source->get_P());
        //            const double s1(S1.norm());
        //            const double s2(S2.norm());
        //            const VectorDim e1(S1.normalized());
        //            const VectorDim e2(S2.normalized());
        //            const VectorDim R12(R2-R1);
        //
        //            const double cosBeta(e1.dot(e2));
        //            const double sinBeta2(1.0-cosBeta*cosBeta);
        //            const double a1=R12.dot(e2*cosBeta-e1)/(sinBeta2+tol);
        //            const double a2=R12.dot(e2-e1*cosBeta)/(sinBeta2+tol);
        //            const double a0=R12.dot(e1.cross(e2))/(sinBeta2+tol);
        //
        ////            //                    std::cout<<"R1="<<std::setprecision(15)<<std::scientific<<R1.transpose()<<std::endl;
        ////            //                    std::cout<<"R2="<<std::setprecision(15)<<std::scientific<<R2.transpose()<<std::endl;
        ////            //                    std::cout<<"S1="<<std::setprecision(15)<<std::scientific<<S1.transpose()<<std::endl;
        ////            //                    std::cout<<"S2="<<std::setprecision(15)<<std::scientific<<S2.transpose()<<std::endl;
        ////            std::cout<<"R1="<<std::setprecision(15)<<R1.transpose()<<std::endl;
        ////            std::cout<<"R2="<<std::setprecision(15)<<R2.transpose()<<std::endl;
        ////            std::cout<<"S1="<<std::setprecision(15)<<S1.transpose()<<std::endl;
        ////            std::cout<<"S2="<<std::setprecision(15)<<S2.transpose()<<std::endl;
        ////            const double temp= F(a1+s1,a2+s2,cosBeta,sinBeta2,a0)
        ////            /*             */ -F(a1+s1,a2   ,cosBeta,sinBeta2,a0)
        ////            /*             */ -F(a1   ,a2+s2,cosBeta,sinBeta2,a0)
        ////            /*             */ +F(a1   ,a2   ,cosBeta,sinBeta2,a0);
        ////            std::cout<<"tamp="<<std::setprecision(15)<<std::scientific<<temp<<std::endl;
        //
        //            return F(a1+s1,a2+s2,cosBeta,sinBeta2,a0)
        //            /* */ -F(a1+s1,a2   ,cosBeta,sinBeta2,a0)
        //            /* */ -F(a1   ,a2+s2,cosBeta,sinBeta2,a0)
        //            /* */ +F(a1   ,a2   ,cosBeta,sinBeta2,a0);
        //        }
        
        
        /**********************************************************************/
        static std::pair<double,bool> linkingNumber(const LoopType& loop1,
                                                    const LoopType& loop2)
        {
            double Ln=0.0;
            
            bool loopsTouch=false;
            
            for(const auto& link1 : loop1.loopLinks())
            {
                for(const auto& link2 : loop2.loopLinks())
                {
                    if(SegmentSegmentDistance<dim>(link1->source->get_P(),
                                                   link1->sink->get_P(),
                                                   link2->source->get_P(),
                                                   link2->sink->get_P()).dMin<FLT_EPSILON)
                    {
                        loopsTouch=true;
                    }
                    
                    Ln +=  segmentPairLN(*link1,*link2);
                    
                }
            }
            
            return std::make_pair(Ln,loopsTouch);
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
                        const std::pair<double,bool> Ln(linkingNumber(*loopIter1->second.lock(),*loopIter2->second.lock()));
                        
                        os  <<loopIter1->first<<" "<<loopIter2->first;
                        os  <<" "<<Ln.second<<" "<<Ln.first;
                        os  <<" "<<loopIter1->second.lock()->burgers().dot(loopIter2->second.lock()->burgers());
                        //                        os  <<" "<<loopIter1->second->burgers().transpose();
                        //                        os  <<" "<<loopIter2->second->burgers().transpose();
                        os  <<"\n";
                    }
                }
            }
            
            return os;
        }
        
    };
    
}
#endif

