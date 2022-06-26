/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po       <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_CrossSlipModels_H
#define model_CrossSlipModels_H

#include <deque>
#include <map>
#include <tuple>

#include <Eigen/Dense>

#include <TypeTraits.h>
#include <LatticeModule.h>
#include <BCClattice.h>
#include <FCClattice.h>
#include <Grain.h>
#include <TerminalColors.h>

namespace model
{

    template <typename DislocationNetworkType>
    struct BaseCrossSlipModel
    {
        
        typedef TypeTraits<DislocationNetworkType> TraitsType;
        static constexpr int dim=TraitsType::dim;
        typedef typename TraitsType::NetworkLinkType NetworkLinkType;
        typedef typename TraitsType::NetworkNodeType NetworkNodeType;
        typedef std::tuple<std::shared_ptr<NetworkNodeType>,std::shared_ptr<NetworkNodeType>,size_t,size_t> CrossSlipTupleType;
        typedef std::deque<CrossSlipTupleType> CrossSlipContainerType;
        
        
        virtual void addToCrossSlip(const NetworkLinkType& link,
                                    CrossSlipContainerType& crossSlipDeq,
                                    const double& r) =0;
        
    };

    template <typename DislocationNetworkType>
    struct AthermalCrossSlipModel : public BaseCrossSlipModel<DislocationNetworkType>
    {
        
        typedef TypeTraits<DislocationNetworkType> TraitsType;
        static constexpr int dim=TraitsType::dim;
        typedef typename TraitsType::NetworkLinkType NetworkLinkType;
        typedef typename TraitsType::NetworkNodeType NetworkNodeType;
        typedef std::tuple<std::shared_ptr<NetworkNodeType>,std::shared_ptr<NetworkNodeType>,size_t,size_t> CrossSlipTupleType;
        typedef std::deque<CrossSlipTupleType> CrossSlipContainerType;
        
        
        AthermalCrossSlipModel();
        
        void addToCrossSlip(const NetworkLinkType& link,
                            CrossSlipContainerType& crossSlipDeq,
                            const double& r) override;
    };

    template <typename DislocationNetworkType>
    struct EscaigCrossSlipModelHEX : public BaseCrossSlipModel<DislocationNetworkType>
    {
        
        typedef TypeTraits<DislocationNetworkType> TraitsType;
        static constexpr int dim=TraitsType::dim;
        typedef typename TraitsType::NetworkLinkType NetworkLinkType;
        typedef typename TraitsType::NetworkNodeType NetworkNodeType;
        typedef std::tuple<std::shared_ptr<NetworkNodeType>,std::shared_ptr<NetworkNodeType>,size_t,size_t> CrossSlipTupleType;
        typedef std::deque<CrossSlipTupleType> CrossSlipContainerType;
        
        const bool enableBasalSlipSystems;
        const bool enablePrismaticSlipSystems;
        const bool enablePyramidalSlipSystems;
        
        
        const double csw;
        const double csLb;
        const double csLp;
        
        const double csb2pE;
        const double csb2pVeg;
        const double csb2pVsc;
        const double csb2pVec;
        
        const double csp2bE;
        const double csp2bVeg;
        const double csp2bVsc;
        const double csp2bVec;
        
        EscaigCrossSlipModelHEX(const PolycrystallineMaterialBase& material);
        
        
        double crossSlipRateKernel(const int& k,const NetworkLinkType& link,const SlipSystem& css,
                                   const double& w,
                                   const double& L,
                                   const double& dE,
                                   const double& Veg,
                                   const double& Vsc,
                                   const double& Vec) const;
        
        void addToCrossSlip(const NetworkLinkType& link,
                            CrossSlipContainerType& crossSlipDeq,
                            const double& r) override;
        
    };



}
#endif
