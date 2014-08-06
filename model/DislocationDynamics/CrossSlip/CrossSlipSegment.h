/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Tamer Crosby     <tcrosby@ucla.edu>.
 * Copyright (C) 2011 by Can Erel         <canerel55@gmail.com>,
 * Copyright (C) 2011 by Giacomo Po       <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_CROSSSLIPSEGMENT_H
#define model_CROSSSLIPSEGMENT_H

#include <math.h> // isfinite
#include <float.h>
#include <list>
#include <stdlib.h> // rand()
#include <ctime> 
#include <time.h>
#include <iterator>
#include <vector>
#include <set>
#include <map>
//#include <boost/random/mersenne_twister.hpp>
//#include <boost/random/uniform_int_distribution.hpp>
#include <model/DislocationDynamics/Materials/Material.h>

namespace model {
    
//    boost::random::mt19937 gen(time(0));
    
    template <typename DislocationSegmentType>
	class CrossSlipSegment
    {
		
		typedef typename DislocationSegmentType::VectorDim VectorDim;
        typedef std::vector<Eigen::Matrix<double,3,1>> VectorVectorDim;
        typedef std::pair<std::vector<double>,std::vector<double> > VectorPairType;
        typedef std::vector<VectorPairType> vector_VectorPairType;
        
        /**********************************************************************/
        VectorDim getConjugateNormalDeterministic(const DislocationSegmentType& ds,
                                                  const double& sinThetaCrossSlipCr,
                                                  const double& crossSlipLength) const
        {/*!@param[in] ds const reference to a DislocationSegment
          * @param[in] sinThetaCrossSlipCr sine (in rads) of the tolerance in the cross-slip angle
          * @param[in] crossSlipLength minimum length of the segment required for cross-slip
          *\returns the unit normal of the cross-slip plane. If this vector is the current
          * plane normal of ds, cross-slip does not take place.
          */
            VectorDim temp(normalPrimary);
            
            if ( !sourceOnMeshBoundary && !sinkOnMeshBoundary 
                && chord.normalized().cross(Burgers.normalized()).norm()<=sinThetaCrossSlipCr
                && !isSessile
                && chord.norm()>1.1*crossSlipLength)
            {
                VectorVectorDim allNormals(ds.conjugatePlaneNormal());
                // add normalPrimary at the beginning of allNormals.
                // This way, in case of duplicate keys, normalPrimary is inserted in argMap
                allNormals.insert(allNormals.begin(),normalPrimary);
                std::map<double,int> argMap; // map automatically sorts keys
                
                for (unsigned int i=0; i< allNormals.size(); i++)
                {
                    const double trss((pkForce-pkForce.dot(allNormals[i])*allNormals[i]).norm());
                    argMap.insert(std::make_pair(trss,i)); // normalPrimary
                }
                
                temp= allNormals[argMap.rbegin()->second];
            }
            return temp;
        }
        
    public:
        
		/*const*/  size_t sourceID;
		/*const*/  size_t   sinkID;
		/*const*/  VectorDim chord;
		/*const*/  VectorDim midPoint;
		/*const*/  VectorDim Burgers;
		/*const*/  bool sourceOnMeshBoundary;
		/*const*/  bool sinkOnMeshBoundary;
		/*const*/  VectorDim pkForce;
		/*const*/  VectorDim normalPrimary;
        /*const*/  bool isSessile;
		/*const*/  VectorDim normalConjugate;
		/*const*/  bool isCrossSlipSegment;
        
		/* Constructor *******************************************************/		
		CrossSlipSegment(const DislocationSegmentType& ds, const double& sinThetaCrossSlipCr,const double& crossSlipLength) : 
		/*init list   */ sourceID(ds.source->sID),
		/*init list   */ sinkID(ds.sink->sID),
		/*init list   */ chord(ds.chord()),
		/*init list   */ midPoint(ds.get_r(0.5)),
		/*init list   */ Burgers(ds.Burgers),
		/*init list   */ sourceOnMeshBoundary(ds.source->meshLocation() == onMeshBoundary),
		/*init list   */ sinkOnMeshBoundary  (ds.sink  ->meshLocation() == onMeshBoundary),
		/*init list   */ pkForce(ds.integralPK()),
		/*init list   */ normalPrimary(ds.glidePlaneNormal),
        /*init list   */ isSessile(std::fabs(Burgers.dot(normalPrimary))>FLT_EPSILON),
        /*init list   */ normalConjugate(getConjugateNormalDeterministic(ds,sinThetaCrossSlipCr,crossSlipLength)),
		/*init list   */ isCrossSlipSegment(normalConjugate != normalPrimary) 
        {
            
                        
        }
        
    };
    
    //////////////////////////////////////////////////////////////s
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


