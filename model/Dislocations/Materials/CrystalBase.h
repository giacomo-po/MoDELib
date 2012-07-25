/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_CRYSTALBASE_H_
#define model_CRYSTALBASE_H_

#include <set>
#include <math.h>
#include <float.h>
#include <assert.h>

#include <Eigen/Core>
#include <model/Dislocations/Materials/SlipSystem.h>


namespace model {
	

	
	
	
	///////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////
	// GENERAL CASE
	template<short unsigned int dim, short unsigned int Nslips>
	class CrystalBase /*: public std::set<model::SlipSystem<dim,Nslips> >*/{
		
		typedef Eigen::Matrix<double,dim,1> VectorDim;
		typedef Eigen::Matrix<double,dim,Nslips> MatrixDimNslips;
		
		typedef SlipSystem<dim,Nslips> SlipSystemType;
        typedef std::set<SlipSystemType> SlipSystemContainerType;
        static SlipSystemContainerType slipSystemContainer;
        
	public:
		
        /* conjugatePlaneNormal ************************************************************************************/
		static void find_slipSystem(const VectorDim  & chord, const VectorDim & Burgers,
							 std::set<SlipSystem<dim,Nslips> > & allowedSlipSystems, const double& tol = FLT_EPSILON){

			assert(chord.norm()>tol && "CHORD HAS ZERO NORM");
			assert(Burgers.norm()>tol && "BURGERS HAS ZERO NORM");
			
			
			VectorDim normalizedChord(chord.normalized());
	
			allowedSlipSystems.clear();
				
            // Try to find a plane which has normal orthogonal to both the chord and the Burgers
			for (typename std::set<SlipSystem<dim,Nslips> >::iterator iter=slipSystemContainer.begin();iter!=slipSystemContainer.end();++iter){
				if(	  std::fabs( iter->normal.dot(normalizedChord))<tol && std::fabs( iter->normal.dot(Burgers.normalized()))<tol){
					allowedSlipSystems.insert( *iter );
				}
			}
			
			unsigned int N(allowedSlipSystems.size());
			switch (N) {
				case 0: // CHECK FOR SESSILE
					for (typename SlipSystemContainerType::const_iterator iter=slipSystemContainer.begin();iter!=slipSystemContainer.end();++iter){
//						std::cout<<"c*n="<< std::fabs( iter->normal.dot(normalizedChord)) << " tol is "<<tol<<std::endl;
						if(	std::fabs( iter->normal.dot(normalizedChord))<tol ){
							allowedSlipSystems.insert( *iter );
						}
					}
//					if (allowedSlipSystems.size()==0){
                    if (allowedSlipSystems.size()<2){
						std::cout<<" chord is"<<chord.transpose()<<std::endl;
						for (typename SlipSystemContainerType::const_iterator iter=slipSystemContainer.begin();iter!=slipSystemContainer.end();++iter){
							std::cout<<"n="<<iter->normal.transpose()<<" |c*n|="<< std::fabs( iter->normal.dot(normalizedChord)) << " tol is "<<tol<<std::endl;
						}
						assert(allowedSlipSystems.size()>=2 && "SESSILE SEGMENTS MUST FORM ON THE INTERSECTION OF TWO CRYSTALLOGRAPHIC PLANES.");					
					}
				//	std::cout<<"CRYSTAL.h: FOUND SESSILE SEGMENT"<<std::endl;
					break;
				case 1: // OK
					break;
				case 2: // CROSS-SLIP SEGMENT
					break;
				default:
					assert(0 && "NUMBER OF SLIP SYSTEMS>2");
					break;
			}

			
		}
		
        
        /* conjugatePlaneNormal ************************************************************************************/
        static VectorDim conjugatePlaneNormal(const VectorDim& B, const VectorDim& N, const double& tol=FLT_EPSILON) {
            
			int count(0);
			VectorDim temp(VectorDim::Zero());
			for (typename SlipSystemContainerType::const_iterator iter=slipSystemContainer.begin();iter!=slipSystemContainer.end();++iter){
				
				for (int k=0;k<Nslips;++k){
					if(	 B.normalized().cross(iter->slip.col(k)).norm()<tol 
					   && N.normalized().cross(iter->normal).norm()>tol){
						temp=iter->normal;
						++count;
					}
				}
			}
			assert(count==1 && "FOUND MORE THAN ONE CONJUGATE PLANES"); // IN BCC THERE IS MORE THAN ONE PLANE!
			return temp;
		}
        
        
        
        /* insert *************************************************************/
        static std::pair<typename SlipSystemContainerType::iterator,bool> insert(const SlipSystemType& x){
            return slipSystemContainer.insert(x);
        }
        
        /* clear *************************************************************/ // CHANGE THIS AND ONLY ALLOW "ROTATE"
        static void clear(){
            slipSystemContainer.clear();
        }        
		
		
		
	};
	
		// declare static data members
    template<short unsigned int dim, short unsigned int Nslips>
    std::set<model::SlipSystem<dim,Nslips> > CrystalBase<dim,Nslips>::slipSystemContainer;
	
	
	
	//////////////////////////////////////////////////////////////
} // namespace model 
#endif




