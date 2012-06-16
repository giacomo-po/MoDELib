/* This file is part of mmdl, the Mechanics of Material Defects Library.
 *
 * Copyright (C) 2011 by Giacomo Po <giacomopo@gmail.com>.
 *
 * mmdl is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef mmdl_GRAMSCHMIDT_H_
#define mmdl_GRAMSCHMIDT_H_


#include <vector>
#include <Eigen/Dense>


namespace mmdl {
	
	
	template <short unsigned int dim>
	class GramSchmidt {
		
		/*! \brief A class template that performs Gram-Schmidt orthonormalization 
		 *  of a set of vectors in arbitrary dimension.
		 */
		
		
		typedef Eigen::Matrix<double,dim,1> VectorDim;


		
	public:
		
		GramSchmidt(const std::vector<Eigen::Matrix<double,dim,1> > & NV, std::vector<Eigen::Matrix<double,dim,1> > & GS)  {
			
			GS=NV;
			
			for (size_t i=0;i<NV.size();++i){
				
				for (size_t j=0;j<i;++j){
					GS[i]-= NV[i].dot(GS[j])*GS[j];
				}
				if (GS[i].squaredNorm()>FLT_EPSILON){
					GS[i].normalize();
				}
				else{
					GS[i].setZero();
				}
			}
			
			//			std::cout<<"Node "<<this->sID<< " size GS = "<<GS.size()<<std::endl;
			
			
			for (typename std::vector<VectorDim>::iterator iter=GS.begin();iter!=GS.end();){
				if (iter->squaredNorm()==0.0){
					//					std::cout<<"I'm here 5"<<std::endl;
					
					iter=GS.erase(iter); // TEST THIS !!!
					//					std::cout<<"I'm here 6"<<std::endl;
					
				}
				else{
					++iter;
				}
			}
			
			assert(GS.size()<=dim);
		}

		
	};
	
}
#endif

