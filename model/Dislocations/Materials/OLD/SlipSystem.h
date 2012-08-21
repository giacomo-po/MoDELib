/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_SLIPSYSTEM_H_
#define model_SLIPSYSTEM_H_

#include <assert.h>
#include <Eigen/Dense>
#include <model/Utilities/StaticID.h>

namespace model {
	
	//////////////////////////////////////////////////////////
	// SlipSystem
	template<short unsigned int dim, short unsigned int Nslips>
	class SlipSystem : public model::StaticID<SlipSystem<dim,Nslips> >{
		
		typedef Eigen::Matrix<double,dim,1>			VectorDim;
		typedef Eigen::Matrix<double,dim,Nslips>	MatrixDimNslips;
		
		//////////////////////////////////////////////////////////////////////////
		static MatrixDimNslips normalizeByColumns(const MatrixDimNslips& slip_in){
			MatrixDimNslips temp;
			for(int k=0; k<Nslips;++k){
				assert(slip_in.col(k).squaredNorm()>0.0);
				temp.col(k)=slip_in.col(k).normalized();
			}
			return temp;
		}

		
	public:
		
		const VectorDim normal;
		const MatrixDimNslips slip;
//		const MatrixDimNslips slip;

		
		//////////////////////////////////
		//Constructor
		SlipSystem(const VectorDim & normal_in, const MatrixDimNslips & slip_in) : normal(normal_in.normalized()), slip(normalizeByColumns(slip_in)){
			
			double tol=1.0e-14;
			
			for(int k=0; k<Nslips;++k){
//				std::cout<<"Creating slip system ("<<normal.transpose()<<")["<<slip.col(k).transpose()<<"]"<<std::endl;
				assert(std::fabs(normal.dot(slip.col(k)))<tol);
			}
			
			
		}
		
//		//////////////////////////////////
//		//Destructor
//		~SlipSystem(){
//			
//			
//			for(int k=0; k<Nslips;++k){
//				std::cout<<"Deleting slip system ("<<normal.transpose()<<")["<<slip.col(k).transpose()<<"]"<<std::endl;
//			}
//			
//			
//		}
		
		bool operator<(const SlipSystem & Other) const {
			return this->sID<Other.sID;
		}
		
		
	};
	//////////////////////////////////////////////////////////////
} // namespace model 
#endif
