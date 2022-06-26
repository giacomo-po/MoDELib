/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po       <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_CrossSlipModels_cpp
#define model_CrossSlipModels_cpp

#include <DislocationNetwork.h>
#include <CrossSlipModels.h>

namespace model
{

//    template <typename DislocationNetworkType>
//    struct BaseCrossSlipModel
//    {
//
//        typedef TypeTraits<DislocationNetworkType> TraitsType;
//        static constexpr int dim=TraitsType::dim;
//        typedef typename TraitsType::NetworkLinkType NetworkLinkType;
//        typedef typename TraitsType::NetworkNodeType NetworkNodeType;
//        typedef std::tuple<std::shared_ptr<NetworkNodeType>,std::shared_ptr<NetworkNodeType>,size_t,size_t> CrossSlipTupleType;
//        typedef std::deque<CrossSlipTupleType> CrossSlipContainerType;
//
//
//        virtual void addToCrossSlip(const NetworkLinkType& link,
//                                    CrossSlipContainerType& crossSlipDeq,
//                                    const double& r) =0;
//
//    };

    
        template <typename DislocationNetworkType>
        AthermalCrossSlipModel<DislocationNetworkType>::AthermalCrossSlipModel()
        {
            std::cout<<greenBoldColor<<"Creating AthermalCrossSlipModel"<<defaultColor<<std::endl;
        }
        
        template <typename DislocationNetworkType>
        void AthermalCrossSlipModel<DislocationNetworkType>::addToCrossSlip(const NetworkLinkType& link,
                            CrossSlipContainerType& crossSlipDeq,
                            const double& r)
        {
            if(link.grains().size()==1)
            {
                // Find slip system on which screw segment has highest glide force
                const auto& grain(**link.grains().begin());
                const auto pki=link.pkIntegral();
                std::map<double,int> ssMap;
                for(size_t s=0;s<grain.slipSystems().size();++s)
                {
                    const auto& slipSystem(grain.slipSystems()[s]);
                    //const VectorDim ncs=slipSystem->unitNormal.normalized();
                    if((slipSystem->s.cartesian()-link.burgers()).norm()<FLT_EPSILON)
                    {
                        const double pkGlide=(pki-pki.dot(slipSystem->unitNormal)*slipSystem->unitNormal).norm();
                        ssMap.emplace(pkGlide,s);
                    }
                }
                
                if(ssMap.size())
                {
                    const auto& crosSlipSystem(grain.slipSystems()[ssMap.rbegin()->second]); // last element in map has highest pkGlide
                    if(crosSlipSystem->unitNormal.cross(link.glidePlaneNormal()).squaredNorm()>FLT_EPSILON)
                    {// cross slip system is different than original system
                        crossSlipDeq.emplace_back(link.source,
                                                  link.sink,
                                                  grain.grainID,
                                                  ssMap.rbegin()->second);
                    }
                }
            }
        }
        
    


        template <typename DislocationNetworkType>
        EscaigCrossSlipModelHEX<DislocationNetworkType>::EscaigCrossSlipModelHEX(const PolycrystallineMaterialBase& material) :
        /* init */ enableBasalSlipSystems(TextFileParser(material.materialFile).readScalar<int>("enableBasalSlipSystems",false))
        /* init */,enablePrismaticSlipSystems(TextFileParser(material.materialFile).readScalar<int>("enablePrismaticSlipSystems",false))
        /* init */,enablePyramidalSlipSystems(TextFileParser(material.materialFile).readScalar<int>("enablePyramidalSlipSystems",false))
        /* init */,csw(TextFileParser(material.materialFile).readScalar<double>("csw_SI",true)*material.b_SI/material.cs_SI)
        /* init */,csLb(enableBasalSlipSystems? TextFileParser(material.materialFile).readScalar<double>("csLb_SI",true)/material.b_SI : 0.0)
        /* init */,csLp(enablePrismaticSlipSystems? TextFileParser(material.materialFile).readScalar<double>("csLp_SI",true)/material.b_SI : 0.0)
        /* init */,csb2pE(enableBasalSlipSystems && enablePrismaticSlipSystems? TextFileParser(material.materialFile).readScalar<double>("csb2pE_eV",true)*material.eV2J/(material.mu_SI*std::pow(material.b_SI,3)) : 0.0)
        /* init */,csb2pVeg(enableBasalSlipSystems && enablePrismaticSlipSystems? TextFileParser(material.materialFile).readScalar<double>("csb2pVeg_SI",true)/std::pow(material.b_SI,3) : 0.0)
        /* init */,csb2pVsc(enableBasalSlipSystems && enablePrismaticSlipSystems? TextFileParser(material.materialFile).readScalar<double>("csb2pVsc_SI",true)/std::pow(material.b_SI,3) : 0.0)
        /* init */,csb2pVec(enableBasalSlipSystems && enablePrismaticSlipSystems? TextFileParser(material.materialFile).readScalar<double>("csb2pVec_SI",true)/std::pow(material.b_SI,3) : 0.0)
        /* init */,csp2bE(enableBasalSlipSystems && enablePrismaticSlipSystems? TextFileParser(material.materialFile).readScalar<double>("csp2bE_eV",true)*material.eV2J/(material.mu_SI*std::pow(material.b_SI,3)) : 0.0)
        /* init */,csp2bVeg(enableBasalSlipSystems && enablePrismaticSlipSystems? TextFileParser(material.materialFile).readScalar<double>("csp2bVeg_SI",true)/std::pow(material.b_SI,3) : 0.0)
        /* init */,csp2bVsc(enableBasalSlipSystems && enablePrismaticSlipSystems? TextFileParser(material.materialFile).readScalar<double>("csp2bVsc_SI",true)/std::pow(material.b_SI,3) : 0.0)
        /* init */,csp2bVec(enableBasalSlipSystems && enablePrismaticSlipSystems? TextFileParser(material.materialFile).readScalar<double>("csp2bVec_SI",true)/std::pow(material.b_SI,3) : 0.0)
        {
            std::cout<<greenBoldColor<<"Creating EscaigCrossSlipModelHEX"<<defaultColor<<std::endl;
        }
        
        template <typename DislocationNetworkType>
        double EscaigCrossSlipModelHEX<DislocationNetworkType>::crossSlipRateKernel(const int& k,const NetworkLinkType& link,const SlipSystem& css,
                                   const double& w,
                                   const double& L,
                                   const double& dE,
                                   const double& Veg,
                                   const double& Vsc,
                                   const double& Vec) const
        { /*!@param[in] k the current quadrature point
           *
           */
            const double eg(link.slipSystem()->unitNormal.transpose()*link.quadraturePoint(k).stress*link.slipSystem()->unitNormal.cross(link.slipSystem()->s.cartesian().normalized()));
            const double ss(css.unitNormal.transpose()*link.quadraturePoint(k).stress*css.s.cartesian().normalized());
            const double ec(css.unitNormal.transpose()*link.quadraturePoint(k).stress*css.unitNormal.cross(css.s.cartesian().normalized()));
            const double dG(dE-std::fabs(ss)*Vsc+eg*Veg-ec*Vec);
            return w*std::exp(-dG/link.network().poly.kB/link.network().poly.T)*link.quadraturePoint(k).j/L;
        }
        
        template <typename DislocationNetworkType>
        void EscaigCrossSlipModelHEX<DislocationNetworkType>::addToCrossSlip(const NetworkLinkType& link,
                            CrossSlipContainerType& crossSlipDeq,
                            const double& r)
        {
            if(link.grains().size()==1)
            {
                const auto& grain(**link.grains().begin());
                std::map<double,int> ssMap;
                for(size_t s=0;s<grain.slipSystems().size();++s)
                {
                    const auto& slipSystem(grain.slipSystems()[s]);
                    if(   (slipSystem->s.cartesian()-link.burgers()).norm()<FLT_EPSILON // slip system contains Burgers
                       && slipSystem->unitNormal.cross(link.glidePlaneNormal()).squaredNorm()>FLT_EPSILON) // not current slip system
                    {
                        
                        if(link.slipSystem()->mobility && slipSystem->mobility)
                        {
                            
                            double csRate(0.0);
                            if(   link.slipSystem()->mobility->name.find("basal") != std::string::npos
                               && slipSystem->mobility->name.find("prismatic") != std::string::npos)
                            {// basal to prismatic crossSlip
                                NetworkLinkType::QuadratureDynamicType::integrate(link.quadraturePoints().size(),this,csRate,&EscaigCrossSlipModelHEX<DislocationNetworkType>::crossSlipRateKernel,link,*slipSystem,
                                                                                  csw,csLb,csb2pE,csb2pVeg,csb2pVsc,csb2pVec);
                            }
                            else if(   link.slipSystem()->mobility->name.find("prismatic") != std::string::npos
                                    && slipSystem->mobility->name.find("basal") != std::string::npos)
                            {// prismatic to basal  crossSlip
                                NetworkLinkType::QuadratureDynamicType::integrate(link.quadraturePoints().size(),this,csRate,&EscaigCrossSlipModelHEX<DislocationNetworkType>::crossSlipRateKernel,link,*slipSystem,
                                                                                  csw,csLp,csp2bE,csp2bVeg,csp2bVsc,csp2bVec);
                            }
                            const double csProb(1.0-std::exp(-csRate*link.network().simulationParameters.dt));
                            std::cout<<"Link "<<link.tag()<<" with SS "<<link.slipSystem()->sID<<" and Length "<<(link.sink->get_P()-link.source->get_P()).norm()<<" has CS "<<s<<" and csRate "<<csRate<<std::endl;
                            ssMap.emplace(csProb,s);
                        }
                    }
                }
                
                if(ssMap.size())
                {
                    std::cout<<"Link "<<link.tag()<<" with SS "<<link.slipSystem()->sID<<" has CS "<<ssMap.rbegin()->second<<" and csProb "<<ssMap.rbegin()->first<<std::endl;
                    std::cout<<"The random number is "<<r<<std::endl;
                    if(ssMap.rbegin()->first > r)
                    {
                        const auto& crosSlipSystem(grain.slipSystems()[ssMap.rbegin()->second]); // last element in map has highest cross-slip probability
                        std::cout<<"crossSlipDeq size "<<crossSlipDeq.size()<<std::endl;
                        crossSlipDeq.emplace_back(link.source,
                                                  link.sink,
                                                  grain.grainID,
                                                  ssMap.rbegin()->second);
                    }
                    
                }
            }
        }
        
    template struct AthermalCrossSlipModel<DislocationNetwork<3,0>>;
    template struct EscaigCrossSlipModelHEX<DislocationNetwork<3,0>>;


}
#endif
