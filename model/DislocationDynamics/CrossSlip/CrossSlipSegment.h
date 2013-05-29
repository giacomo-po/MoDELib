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
        
        
        VectorDim getConjugateNormal(const DislocationSegmentType& ds, const double& sinThetaCrossSlipCr,const double& crossSlipLength)
        {
            VectorDim temp(normalPrimary);
            if ( !sourceOnMeshBoundary && !sinkOnMeshBoundary 
                && chord.normalized().cross(Burgers.normalized()).norm()<=sinThetaCrossSlipCr
                && !isSessile
                && chord.norm()>1.1*crossSlipLength) {
                
                std::set<double> probabilities;
                double ptotal(0.0);
                double random_number(roll_die());
                double mu(160.0e3);
                double b(0.2722e-9);
                double K(1.38e-23/mu/std::pow(b,3));
                double T(6000);
                double tau_o(32.0/mu);
                double Vact(300);
                
                vector_VectorDim AllNormals(ds.conjugatePlaneNormal());
                AllNormals.push_back(normalPrimary);
                
                for (int i=0; i< AllNormals.size(); i++) {
                    double trss((pkForce-pkForce.dot(AllNormals[i])*AllNormals[i]).norm()); 
                    double ptemp( exp( -Vact*(tau_o-trss)/( K*T ) ));  
                    ptotal+=ptemp;
                    probabilities.insert(ptotal);
                }
                
                double r(0.1*random_number*ptotal);
                
                std::set<double>::iterator it(probabilities.lower_bound(r));
                int n(std::distance(probabilities.begin(),it));
        
                temp= AllNormals[n];
                
                std::cout<<"primary normal " <<normalPrimary.transpose()<<std::endl;
                std::cout<<" conjugate normal = "<<temp.transpose()<<std::endl; 
                std::cout<<" random number = "<<random_number<<std::endl; 
                std::cout<<"r = "<<r<<std::endl;
            }
            return temp;
        }
        
    public:
        
		/*const*/ /*const*/ size_t sourceID;
		/*const*/ /*const*/ size_t   sinkID;
		/*const*/ /*const*/ VectorDim chord;
		/*const*/ /*const*/ VectorDim midPoint;
		/*const*/ /*const*/ VectorDim Burgers;
		/*const*/ /*const*/ bool sourceOnMeshBoundary;
		/*const*/ /*const*/ bool sinkOnMeshBoundary;
		/*const*/ /*const*/ VectorDim pkForce;
		/*const*/ /*const*/ VectorDim normalPrimary;
        /*const*/ /*const*/ bool isSessile;
		/*const*/ /*const*/ VectorDim normalConjugate;
		/*const*/ /*const*/ bool isCrossSlipSegment;
        
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
        int roll_die() {
            boost::random::uniform_int_distribution<> dist(1,10);
            return dist(gen);
        }
        
    };
    
    //////////////////////////////////////////////////////////////s
} // namespace model
#endif

