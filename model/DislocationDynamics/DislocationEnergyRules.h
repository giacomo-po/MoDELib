/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DISLOCATIONENERGYRULES_H_
#define model_DISLOCATIONENERGYRULES_H_

//#include <boost/tuple/tuple.hpp>
#include <tuple>
#include <Eigen/Dense>
#include <Eigen/StdVector>

#include <model/DislocationDynamics/EdgeConfigs.h>

namespace model
{
	
	template<short unsigned int dim>
	class DislocationEnergyRules
    {
		
		typedef Eigen::Matrix<double,dim,1> VectorDimD;
		
	public:
		
		
		static double interactionEnergy(const VectorDimD& b1, const VectorDimD& b2, const int& s1, const int& s2, const VectorDimD& t, const double& nu){
			const VectorDimD t1(static_cast<double>(s1)*t);
			const VectorDimD t2(static_cast<double>(s2)*t);
			return b1.dot(t1)*b2.dot(t2)+1.0/(1.0-nu)*(b1.dot(b2)*t1.dot(t2)-b1.dot(t2)*b2.dot(t1));
		}

		
		/* findEdgeConfiguration *********************************************/
		template <typename DislocationNodeType>
		static void findEdgeConfiguration(DislocationNodeType& dN)
        {
//            static void findEdgeConfiguration(DislocationNodeType& dN, const double& nu){
			
            	Eigen::VectorXi edgeConfiguration;

            
			if (dN.is_isolated()){
//				dN.edgeConfiguration.setZero(0);
                edgeConfiguration.setZero(0);
			}
			else{
				// 1- Collect all the Burgers from the neighbors considering them as out-neighbors.
                typedef typename DislocationNodeType::FlowType FlowType;
				std::vector<FlowType,Eigen::aligned_allocator<FlowType> >  BsV{}; // the vector of Burgers as out-neighbors
				for (typename DislocationNodeType::constNeighborIteratorType neighborIter=dN.neighborhood().begin();neighborIter!=dN.neighborhood().end();++neighborIter){
					const int dir=std::get<2>(neighborIter->second);
					switch ( dir ) {
						case   1:
							BsV.push_back((std::get<1>(neighborIter->second)->flow)*(+1));
							break;
						case  -1:
							BsV.push_back((std::get<1>(neighborIter->second)->flow)*(-1));
							break;
						default:	// self
							break;
					}
				}
				
                
                
                const size_t nN(BsV.size());
                
                // 5- reset edgeConfiguration
//				dN.edgeConfiguration.setZero(nN);
				edgeConfiguration.setZero(nN);

                if (nN<=EdgeDynamicConfigs::maxEdge)
                {
                
                    // 2- Get all the possible edge configurations
                    Eigen::MatrixXi Ci(EdgeDynamicConfigs::getCi(BsV.size()));
                    //				assert(Ci.cols()==BsV.size() && "SOMETHING WENT WRONG");
                    
                    
                    // 3- Change the sign of Ci according to the direction of the first neighbor
                    int sigCi=1;
                    for (typename DislocationNodeType::constNeighborIteratorType neighborIter=dN.neighborhood().begin();neighborIter!=dN.neighborhood().end();++neighborIter)
                    {
                        const int dir(std::get<2>(neighborIter->second));
                        if (dir!=0)
                        {
                            sigCi=dir;
                            break;
                        }
                    }
                    Ci*=sigCi;
                    
                    // 4- Store total chord parametric length
                    double gT(0.0);
                    for (typename DislocationNodeType::constNeighborIteratorType neighborIter=dN.neighborhood().begin();neighborIter!=dN.neighborhood().end();++neighborIter)
                    {
                        if (std::get<2>(neighborIter->second)!=0)
                        {
                            gT+=std::get<1>(neighborIter->second)->chordParametricLength();
                        }
                    }
                    
                    int chi=1+!(dN.neighborhood().size()>2);
                    gT*=chi;
                    
                    // 5- Compute and sort energy levels
                    std::multimap<double,Eigen::VectorXi> ELev{};	// MAP DOES NOT ACCEPT EQUAL VALUES SO multimap IS USED TO STORE  DEGENERATE STATES
                    for (int k=0;k<Ci.rows();++k)
                    {
                        
                        // 5a- Compute the Catmull-Rom tangent using the current edge configuration
                        VectorDimD tTemp(VectorDimD::Zero());
                        int j(0);
                        for (typename DislocationNodeType::constNeighborIteratorType neighborIter=dN.neighborhood().begin();neighborIter!=dN.neighborhood().end();++neighborIter)
                        {
                            if (std::get<2>(neighborIter->second)!=0)
                            {
                                double gkj(std::get<1>(neighborIter->second)->chordParametricLength());
                                tTemp+=Ci(k,j)*(std::get<0>(neighborIter->second)->get_P()-dN.get_P())/gkj*(gT-gkj)/gT;
                                j++;
                            }
                        }
                        
                        //					const double tTempNorm(tTemp.norm());
                        //					if (tTempNorm>FLT_EPSILON){
                        //					tTemp/=tTempNorm; //! TANGENTS MUST BE ALL UNIT FOR THE COMPARISON OF THE ENRGY LEVELS TO MAKE SENSE!!!!
                        //					}
                        //					else{
                        //					tTemp.setZero();
                        //					}
                        //
                        //					// 5b- Compute the energy for the current configuration and current tangent
                        //					double Ek(0.0);
                        //					for(unsigned int a=0; a<BsV.size();++a){
                        //						for(unsigned int b=0; b<BsV.size();++b){
                        //							Ek+=interactionEnergy(BsV[a],BsV[b],Ci(k,a),Ci(k,b),tTemp,nu);
                        //						}
                        //					}
                        //					ELev.insert(std::make_pair(Ek,Ci.row(k)));
                        
                        
                        bool isTrivialConfig(true);
                        for (int c=0;c<Ci.cols();++c)
                        {
                            isTrivialConfig*=(Ci(k,c)==Ci(k,0));
                        }
                        
                        ELev.insert(std::make_pair(tTemp.squaredNorm()*(!isTrivialConfig),Ci.row(k))); // the norm of the tangent is the measure 
                    }
                    
                    
                    //				// 5- reset edgeConfiguration
                    //				dN.edgeConfiguration.setZero(Ci.cols());
                    
                    if (dN.is_balanced())
                    {
                        assert(ELev.size()>1 && "MORE THAN ONE ENERGY LEVELS MUST BE FOUND FOR A BALANCED NODE");	
                        //std::multimap<double,Eigen::VectorXi>::const_iterator firstNonZero(ELev.lower_bound(FLT_EPSILON)); // the first element that compares >=FLT_EPSILON
                        //assert(firstNonZero!=ELev.end() && "AT LEAST ONE POSITIVE ENERGY LEVEL MUST EXIST");					
                    //    dN.edgeConfiguration=ELev.rbegin()->second;
                        edgeConfiguration=ELev.rbegin()->second;
                    }
                
                
                }
                

												
				
				// 7- store tangent coefficients
				unsigned int kk(0);
				typename DislocationNodeType::NeighborContainerType tempNeigh=dN.neighborhood(); // DESIGN FLAW: neighborhood should have a non-const version
				for (typename DislocationNodeType::NeighborIteratorType neighborIter=tempNeigh.begin();neighborIter!=tempNeigh.end();++neighborIter){
					int dir=std::get<2>(neighborIter->second);
					switch ( dir ) {
						case   1:	// out
//							std::get<1>(neighborIter->second)->sourceTfactor=dN.edgeConfiguration(kk);
                            std::get<1>(neighborIter->second)->sourceTfactor=edgeConfiguration(kk);
							kk++;
							break;
						case  -1:	// in
//							std::get<1>(neighborIter->second)->  sinkTfactor=dN.edgeConfiguration(kk);
                            std::get<1>(neighborIter->second)->  sinkTfactor=edgeConfiguration(kk);
							kk++;
							break;
						default:	// self
							break;
					} // end switch
				} // end for
			} // end else (not isolated)
		}
		
	};
	/*********************************************************************/
	/*********************************************************************/
} // end namespace
#endif


