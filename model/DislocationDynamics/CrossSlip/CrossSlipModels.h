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
#include <model/LatticeMath/LatticeMath.h>
#include <model/DislocationDynamics/Materials/BCClattice.h>
#include <model/DislocationDynamics/Materials/FCClattice.h>
#include <model/DislocationDynamics/Polycrystals/Grain.h>

namespace model
{
    
    template<typename CrystalType>
    struct CrossSlipModels
    {
        
    };
    

    template<>
    struct CrossSlipModels<FCClattice<3>>
    {
        typedef std::tuple<size_t,size_t,size_t,size_t> CrossSlipTupleType;
        typedef std::deque<CrossSlipTupleType> CrossSlipContainerType;
        
        /**********************************************************************/
        template<typename LinkType>
        static void addToCrossSlip(const LinkType& link,
                                   CrossSlipContainerType& crossSlipDeq,
                                   const int& modelID)
        {
            
            
            switch (modelID)
            {
                case 0:
                {// no cross-slip
                    
                    break;
                }
                    
                case 1:
                {// deterministic cross-slip model based on max glide PK force
                    
                    assert(link.grains().size()==1 && "THERE SHOULD BE ONLY ONE GRAIN FOR AN INTERNAL SEGMENT");
                    
                    // Find slip system on which screw segment has highest glide force
                    const auto& grain(**link.grains().begin());
                    const auto pki=link.pkIntegral();
                    std::map<double,int> ssMap;
                    for(size_t s=0;s<grain.slipSystems().size();++s)
                    {
                        const auto& slipSystem(grain.slipSystems()[s]);
                        //const VectorDim ncs=slipSystem.unitNormal.normalized();
                        if((slipSystem.s.cartesian()-link.burgers()).norm()<FLT_EPSILON)
                        {
                            const double pkGlide=(pki-pki.dot(slipSystem.unitNormal)*slipSystem.unitNormal).norm();
                            ssMap.emplace(pkGlide,s);
                        }
                    }
                    
                    if(ssMap.size())
                    {
                        const auto& crosSlipSystem(grain.slipSystems()[ssMap.rbegin()->second]); // last element in map has highest pkGlide
                        if(crosSlipSystem.unitNormal.cross(link.glidePlaneNormal()).squaredNorm()>FLT_EPSILON)
                        {// cross slip system is different than original system
                            crossSlipDeq.emplace_back(link.source->sID,
                                                      link.sink->sID,
                                                      grain.grainID,
                                                      ssMap.rbegin()->second);
                        }
                    }
                    break;
                }
                    
                    
                default:
                {
                    assert(0 && "CROSS-SLIP MODEL NOT IMPLEMENTED");
                    break;
                }
            }
            
            
        }
    };
    
    template<>
    struct CrossSlipModels<BCClattice<3>>
    {
        
        typedef std::tuple<size_t,size_t,size_t,size_t> CrossSlipTupleType;
        typedef std::deque<CrossSlipTupleType> CrossSlipContainerType;
        
        /**********************************************************************/
        template<typename LinkType>
        static void addToCrossSlip(const LinkType& link,
                                   CrossSlipContainerType& crossSlipDeq,
                                   const int& modelID)
        {
            
            
            switch (modelID)
            {
                case 0:
                {// no cross-slip
                    
                    break;
                }
                    
                case 1:
                {// deterministic cross-slip model based on max glide PK force
                    
                    assert(link.grains().size()==1 && "THERE SHOULD BE ONLY ONE GRAIN FOR AN INTERNAL SEGMENT");
                    
                    // Find slip system on which screw segment has highest glide force
                    const auto& grain(**link.grains().begin());
                    const auto pki=link.pkIntegral();
                    std::map<double,int> ssMap;
                    for(size_t s=0;s<grain.slipSystems().size();++s)
                    {
                        const auto& slipSystem(grain.slipSystems()[s]);
                        //const VectorDim ncs=slipSystem.unitNormal.normalized();
                        if((slipSystem.s.cartesian()-link.burgers()).norm()<FLT_EPSILON)
                        {
                            const double pkGlide=(pki-pki.dot(slipSystem.unitNormal)*slipSystem.unitNormal).norm();
                            ssMap.emplace(pkGlide,s);
                        }
                    }
                    
                    if(ssMap.size())
                    {
                        const auto& crosSlipSystem(grain.slipSystems()[ssMap.rbegin()->second]); // last element in map has highest pkGlide
                        if(crosSlipSystem.unitNormal.cross(link.glidePlaneNormal()).squaredNorm()>FLT_EPSILON)
                        {// cross slip system is different than original system
                            crossSlipDeq.emplace_back(link.source->sID,
                                                      link.sink->sID,
                                                      grain.grainID,
                                                      ssMap.rbegin()->second);
                        }
                    }
                    break;
                }

                    
                default:
                {
                    assert(0 && "CROSS-SLIP MODEL NOT IMPLEMENTED");
                    break;
                }
            }
            

        }
        
    };
    
//    struct CrossSlipModels
//    {
//        
//        /**********************************************************************/
//        template <typename DislocationSegmentType>
//        static const LatticePlaneBase& maxRSS(const DislocationSegmentType& ds)
//        {/*!@param[in] ds const reference to a DislocationSegment
//          *\returns the plane on which the segment moves according to the cross-slip model.
//          */
//            
//            typedef Eigen::Matrix<double,DislocationSegmentType::dim,1> VectorDim;
//            
//            const VectorDim pkForce=ds.integralPK();
//            
//            std::deque<const LatticePlaneBase*> allNormals(ds.conjugatePlaneNormals);
//            // add normalPrimary at the beginning of allNormals.
//            // This way, in case of duplicate keys, normalPrimary is inserted in argMap
//            allNormals.insert(allNormals.begin(),&ds.glidePlane.n);
//            std::map<double,int> argMap; // map automatically sorts keys
//            const VectorDim nc=allNormals[0]->cartesian().normalized();
//            const double trss0((pkForce-pkForce.dot(nc)*nc).norm());
//            argMap.insert(std::make_pair(trss0,0)); // normalPrimary
//            
//            
//            for (unsigned int i=1; i< allNormals.size(); i++)
//            {
//                const VectorDim nc=allNormals[i]->cartesian().normalized();
//                const double trss((pkForce-pkForce.dot(nc)*nc).norm());
//                argMap.insert(std::make_pair(trss,i)); // normalPrimary
//            }
//            
//            return *allNormals[argMap.rbegin()->second];
//        }
//        
//    };
    
} // namespace model
#endif



///**********************************************************************/
//VectorDim getConjugateNormal(const DislocationSegmentType& ds, const double& sinThetaCrossSlipCr,const double& crossSlipLength) const
//{
//    VectorDim temp(normalPrimary);
//    if ( !sourceOnMeshBoundary && !sinkOnMeshBoundary
//        && chord.normalized().cross(Burgers.normalized()).norm()<=sinThetaCrossSlipCr
//        && !isSessile
//        && chord.norm()>1.1*crossSlipLength
//        && Material<Isotropic>::kT > 0.0 )
//    {
//        std::set<double> probabilities;
//        double ptotal(0.0);
//        VectorVectorDim allNormals(ds.conjugatePlaneNormal());
//        // add normalPrimary at the beginning of allNormals.
//        // This way, in case of duplicate keys, normalPrimary is inserted in argMap
//        allNormals.insert(allNormals.begin(),normalPrimary);
//        //                allNormals.push_back(normalPrimary);
//        //                std::cout<<std::endl;
//        std::map<double,int> argMap; // map automatically sorts keys
//
//        for (unsigned int i=0; i< allNormals.size(); i++)
//        {
//            const double trss((pkForce-pkForce.dot(allNormals[i])*allNormals[i]).norm());
//            const double arg(-Material<Isotropic>::vAct*(Material<Isotropic>::tauIII-trss)/( Material<Isotropic>::kT ));
//            const double ptemp( exp(arg));
//            if(!std::isfinite(ptemp)) // arg makes exp(arg) blow up, so store arg itself
//            {
//                argMap.insert(std::make_pair(arg,i)); // normalPrimary
//            }
//
//            ptotal+=ptemp;
//            //                    std::cout<<" ptemp = "<<ptemp<<" conjugate normal = "<<allNormals[i]<<std::endl;
//            probabilities.insert(ptotal);
//        }
//
//        if (argMap.size()>1) // at least one probability is inf
//        {
//            // Pick the highest arg
//            temp= allNormals[argMap.rbegin()->second];
//
//        }
//        else // none of the  probabilities are inf
//        {
//            //                    double random_number(roll_die());
//            //                    double r(0.1*random_number*ptotal);
//            //                    double r(0.1*random_number*ptotal);
//
//            double r(static_cast<double>(std::rand()) / RAND_MAX * ptotal);
//
//            std::set<double>::iterator it(probabilities.lower_bound(r));
//            int n(std::distance(probabilities.begin(),it));
//
//            temp= allNormals[n];
//        }
//        //                std::cout<<"r = "<<r<<std::endl;
//        //                std::cout<<" random number = "<<random_number<<std::endl;
//        //                std::cout<<"ptotal = "<<ptotal<<std::endl;
//        //                std::cout<<"primary normal " <<normalPrimary.transpose()<<std::endl;
//        //                std::cout<<" conjugate normal final = "<<temp.transpose()<<std::endl;
//    }
//
//    return temp;
//}


