/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef  model_GIVENSQR_H_
#define  model_GIVENSQR_H_

#include <Eigen/Dense>
#include <model/Math/GeneralizedEigenSolver_Giacomo/Givens.h>

namespace model {
	
	template<int N, int ImJ, int NmJ>
	struct GivensQRkernel{

		enum {J=N-NmJ};	// J=N-(N-J)
		enum {I=ImJ+J}; // I=(I-J)+J
		enum {NmI=N-I}; // 
		
		static	void step(Eigen::Matrix<double,N,N>& Q, Eigen::Matrix<double,N,N>& R){
			Eigen::Matrix<double,2,2> G(Givens::Q<double>(R(I,J),R(I+1,J)));
			R.template block<2,NmJ>(I,J)=(G*R.template block<2,NmJ>(I,J)).eval();			// use .eval() to avoid aliasing
			Q.template block<N,2>(0,I)=(Q.template block<N,2>(0,I)*G.transpose()).eval();   // use .eval() to avoid aliasing
		}
	};
	
	template<int N, int ImJ, int NmJ>
	struct GivensQRstep{
		/* Define ImJ=I-J 
		 *    and NmJ=N-J
		 * therefore: J=N-NmJ 
		 *            I=ImJ+J=ImJ+N-NmJ
		 * when the current element is trasformed call recursivly decreasing ImJ by one, that is
		 * GivensQRstep<N,ImJ-1,NmJ>
		 */
		static	void step(Eigen::Matrix<double,N,N>& Q, Eigen::Matrix<double,N,N>& R){
			Eigen::Matrix<double,2,2> G(Givens::Q<double>(R(ImJ+N-NmJ,N-NmJ),R(ImJ+N-NmJ+1,N-NmJ)));
			GivensQRkernel<N,ImJ,NmJ>::step(Q,R);
			GivensQRstep<N,ImJ-1,NmJ>::step(Q,R);
		}
	};
	
	template<int N, int NmJ>
	struct GivensQRstep<N,0,NmJ>{
		enum  {ImJ=0};
		static void step(Eigen::Matrix<double,N,N>& Q,Eigen::Matrix<double,N,N>& R){
			GivensQRkernel<N,ImJ,NmJ>::step(Q,R);
			GivensQRstep<N,NmJ-3,NmJ-1>::step(Q,R);	// restart with next column, that is ImJ=I-J=NmJ-1-2 
		}
	};
	
	
	template<int N>
	struct GivensQRstep<N,0,2>{
		enum  {ImJ=0};
		enum  {NmJ=2};
		static void step(Eigen::Matrix<double,N,N>& Q, Eigen::Matrix<double,N,N>& R){
			GivensQRkernel<N,ImJ,NmJ>::step(Q,R);
		}		
	};
	
	
	/*! \brief Class template that performs and stores the QR decomposition of a NxN matrix based on Givens rotation.
	 *
	 * \f[
	 *    a=2
	 *  \f]
	 */
	template<int N>
	class GivensQR{

	protected:
		Eigen::Matrix<double,N,N> Q;
		Eigen::Matrix<double,N,N> R;
		
	public:
		GivensQR(const Eigen::Matrix<double,N,N>& A) : Q(Eigen::Matrix<double,N,N>::Identity()), 
		/* init list                              */   R(A){
			GivensQRstep<N,N-2,N>::step(Q,R);
		}
		
		const Eigen::Matrix<double,N,N>& matrixQ() const {return Q;}
	  const Eigen::Matrix<double,N,N>& matrixR() const {return R;}
	};
	
} // end namespace
#endif
