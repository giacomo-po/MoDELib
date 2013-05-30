/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Can Erel         <canerel55@gmail.com>,
 * Copyright (C) 2011 by Giacomo Po       <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_CROSSSLIPSEGMENT_H
#define model_CROSSSLIPSEGMENT_H

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <model/DislocationDynamics/Materials/Material.h>
#include <math.h>
#include <float.h>
#include <list>
#include <stdlib.h>
#include <ctime> 
#include <time.h>
#include <iterator>

namespace model {
    
    boost::random::mt19937 gen(time(0));
    
    template <typename DislocationSegmentType>
	class CrossSlipSegment {
		
		typedef typename DislocationSegmentType::VectorDim VectorDim;
        typedef std::vector<Eigen::Matrix<double,3,1>> vector_VectorDim;
        typedef std::pair<std::vector<double>,std::vector<double> > vector_pair;
        typedef std::vector<vector_pair> vector_vector_pair;
        
        
        VectorDim getConjugateNormal(const DislocationSegmentType& ds, const double& sinThetaCrossSlipCr,const double& crossSlipLength) const
        {
            VectorDim temp(normalPrimary);
            if ( !sourceOnMeshBoundary && !sinkOnMeshBoundary 
                && chord.normalized().cross(Burgers.normalized()).norm()<=sinThetaCrossSlipCr
                && !isSessile
                && chord.norm()>1.1*crossSlipLength
                && Material<Isotropic>::kT > 0.0 ) {
                
                std::set<double> probabilities;
                double ptotal(0.0);
                
                vector_VectorDim AllNormals(ds.conjugatePlaneNormal());
                AllNormals.push_back(normalPrimary);
                std::cout<<std::endl;

                for (int i=0; i< AllNormals.size(); i++) {
                    const double trss((pkForce-pkForce.dot(AllNormals[i])*AllNormals[i]).norm()); 
                    const double ptemp( exp( -Material<Isotropic>::vAct*(Material<Isotropic>::tauIII-trss)/( Material<Isotropic>::kT ) ));  
                    ptotal+=ptemp;
//                    std::cout<<" ptemp = "<<ptemp<<" conjugate normal = "<<AllNormals[i]<<std::endl;
                    probabilities.insert(ptotal);
                }
                
                double random_number(roll_die());
                double r(0.1*random_number*ptotal);
                
                std::set<double>::iterator it(probabilities.lower_bound(r));
                int n(std::distance(probabilities.begin(),it));
        
                temp= AllNormals[n];
                
//                std::cout<<"r = "<<r<<std::endl;
//                std::cout<<" random number = "<<random_number<<std::endl; 
//                std::cout<<"ptotal = "<<ptotal<<std::endl;
//                std::cout<<"primary normal " <<normalPrimary.transpose()<<std::endl;
//                std::cout<<" conjugate normal final = "<<temp.transpose()<<std::endl; 
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
		/*init list   */ sourceOnMeshBoundary(ds.source->nodeMeshLocation == onMeshBoundary),
		/*init list   */ sinkOnMeshBoundary(ds.sink  ->nodeMeshLocation == onMeshBoundary),
		/*init list   */ pkForce(ds.integralPK()),
		/*init list   */ normalPrimary(ds.glidePlaneNormal),
        /*init list   */ isSessile(std::fabs(Burgers.dot(normalPrimary))>FLT_EPSILON),
        /*init list   */ normalConjugate(getConjugateNormal(ds,sinThetaCrossSlipCr,crossSlipLength)),
		/*init list   */ isCrossSlipSegment(normalConjugate != normalPrimary) 

        {
            
                        
        }
        
        /*-----------------------------------------------------------------------*/
        int roll_die() const {
            boost::random::uniform_int_distribution<> dist(1,10);
            return dist(gen);
        }
        
    };
    
    //////////////////////////////////////////////////////////////s
} // namespace model
#endif

