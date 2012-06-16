/* This file is part of mmdl, the Mechanics of Material Defects Library.
 *
 * Copyright (C) 2011 by Giacomo Po <giacomopo@gmail.com>.
 *
 * mmdl is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef  mmdl_GIVENSHT_H_
#define  mmdl_GIVENSHT_H_



#include <Eigen/Dense>
#include <mmdl/Math/GeneralizedEigenSolver_Giacomo/GivensQR.h>

namespace mmdl {
	
	template<int N, int ImJ, int NmJ> 
	struct BaseGivensHTstep{
		enum {J=N-NmJ}; // J=N-(N-J)
		enum {I=ImJ+J}; // I=ImJ+J=ImJ+N-NmJ
		enum {NmI=N-I}; 
		static	void step(Eigen::Matrix<double,N,N>& H, Eigen::Matrix<double,N,N>& T){
			Eigen::Matrix<double,2,2> G(Givens::Q<double>(H(I,J),H(I+1,J)));
			H.template block<2,NmJ>(I,J)=(G*H.template block<2,NmJ>(I,J)).eval();	// use .eval() to avoid aliasing
			T.template block<2,NmI>(I,I)=(G*T.template block<2,NmI>(I,I)).eval();	// use .eval() to avoid aliasing
			G=Givens::Q<double>(-T(I+1,I+1),T(I+1,I)).transpose();
			T.template block<I+2,2>(0,I)=(T.template block<I+2,2>(0,I)*G).eval();	// use .eval() to avoid aliasing
			H.template block<N,2>(0,I)=(H.template block<N,2>(0,I)*G).eval();	    // use .eval() to avoid aliasing
		}
	};
	
	template<int N, int ImJ, int NmJ>
	struct GivensHTstep{
		/* Define ImJ=I-J 
		 *    and NmJ=N-J
		 * therefore: J=N-NmJ 
		 *            I=ImJ+J=ImJ+N-NmJ
		 
		 * N-I=N-ImJ-N+NmJ=NmJ-ImJ
		 * when the current element is trasformed call recursivly decreasing ImJ by one, that is
		 * GivensQRstep<N,ImJ-1,NmJ>
		 */
		static	void step(Eigen::Matrix<double,N,N>& H, Eigen::Matrix<double,N,N>& T){
			BaseGivensHTstep<N,ImJ,NmJ>::step(H,T);			
			GivensHTstep<N,ImJ-1,NmJ>::step(H,T);
		}
	};
	
	template<int N, int NmJ>
	struct GivensHTstep<N,1,NmJ>{
		enum {ImJ=1};
		static void step(Eigen::Matrix<double,N,N>& H, Eigen::Matrix<double,N,N>& T){			
			BaseGivensHTstep<N,ImJ,NmJ>::step(H,T);
			GivensHTstep<N,NmJ-1-2,NmJ-1>::step(H,T);	// restart with next column
		}
	};
	
	template<int N>
	struct GivensHTstep<N,1,3>{ 
		enum {NmJ=3};
		enum {ImJ=1};
		static void step(Eigen::Matrix<double,N,N>& H, Eigen::Matrix<double,N,N>& T){
			BaseGivensHTstep<N,ImJ,NmJ>::step(H,T);
		}
	};
	
	template<int N>
	class GivensHT : public GivensQR<N>{
		/*! Given fixed-size and square matrices A and B, GivensHT finds the matrices H and T such that:
		 *	H=Q*A*Z is upper Hessenberg
		 *  T=Q*B*Z is upper triangular
		 *  Q and Z are orthogonal
		 */
		
		Eigen::Matrix<double,N,N> Z;
		Eigen::Matrix<double,N,N> T;
		Eigen::Matrix<double,N,N> H;

	public:
		GivensHT(const Eigen::Matrix<double,N,N>& A, 
		/*    */ const Eigen::Matrix<double,N,N>& B) : GivensQR<N>(B), // this will initialize Q and R
		/* init list                              */   Z(Eigen::Matrix<double,N,N>::Identity()),
		/* init list                              */   T(this->matrixR()),
		/* init list                              */   H(this->matrixQ().transpose()*A){
			GivensHTstep<N,N-2,N>::step(H,T);
		}
		
		const Eigen::Matrix<double,N,N>& matrixH(){return H;}
	    const Eigen::Matrix<double,N,N>& matrixT(){return T;}
	    const Eigen::Matrix<double,N,N>& matrixZ(){return Z;}
	};
	
}
#endif

//		const Eigen::Matrix<double,N,N>& matrixQ(){return Q;}


//			T.template block<N,2>(0,I)=(T.template block<N,2>(0,I)*G).eval();	// use .eval() to avoid aliasing


//			std::cout<<"(I,J)=("<<I<<","<<J<<")"<<std::endl;
//			std::cout<<"(-T(I+1,I+1),T(I+1,I))=("<<-T(I+1,I+1)<<","<<T(I+1,I)<<")"<<std::endl;
//			std::cout<<"T="<<std::endl<<T<<std::endl;
//			
//			std::cout<<"Tblock="<<std::endl<<T.template block<I,2>(0,I)<<std::endl;

