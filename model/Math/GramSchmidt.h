/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_GRAMSCHMIDT_H_
#define model_GRAMSCHMIDT_H_


#include <vector>
#include <Eigen/Dense>
#include <Eigen/StdVector>


namespace model {
	
	
	template <short unsigned int dim>
	class GramSchmidt : public std::vector<Eigen::Matrix<double,dim,1>,Eigen::aligned_allocator<Eigen::Matrix<double,dim,1> > >{
		
		/*! \brief A class template that performs Gram-Schmidt orthonormalization 
		 *  of a set of vectors in arbitrary dimension.
		 */
		
		
		typedef Eigen::Matrix<double,dim,1> VectorDim;		
		typedef std::vector<VectorDim,Eigen::aligned_allocator<VectorDim> > VectorOfNormalsType;
//		VectorOfNormalsType planeNormals;

		


		
	public:
		
//		GramSchmidt(const std::vector<Eigen::Matrix<double,dim,1> > & NV, std::vector<Eigen::Matrix<double,dim,1> > & GS) : VectorOfNormalsType(NV) {

		/********************************************************************/
		GramSchmidt() {}

		/********************************************************************/		
		GramSchmidt(const VectorOfNormalsType& NV) {
			orthoNormalize(NV);
		}
		
		
		/********************************************************************/		
		void orthoNormalize(const VectorOfNormalsType& NV){
			VectorOfNormalsType::operator=(NV);
			
			for (size_t i=0;i<NV.size();++i){
				
				for (size_t j=0;j<i;++j){
					this->operator[](i)-= NV[i].dot(this->operator[](j))*this->operator[](j);
				}
				if (this->operator[](i).squaredNorm()>FLT_EPSILON){
					this->operator[](i).normalize();
				}
				else{
					this->operator[](i).setZero();
				}
			}
			
			//			std::cout<<"Node "<<this->sID<< " size GS = "<<GS.size()<<std::endl;
			
			
			for (typename VectorOfNormalsType::iterator iter=this->begin();iter!=this->end();){
				if (iter->squaredNorm()==0.0){
					//					std::cout<<"I'm here 5"<<std::endl;
					
					iter=this->erase(iter); // TEST THIS !!!
					//					std::cout<<"I'm here 6"<<std::endl;
					
				}
				else{
					++iter;
				}
			}
			
			assert(this->size()<=dim);
		}

		
	};
	
}
#endif

