/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_CombinationWithoutRepetition_H_
#define model_CombinationWithoutRepetition_H_

#include <assert.h>
#include <Eigen/Core>
#include <model/CompileTimeMath/Binomial.h>

namespace model {
	
	
	// THIS SHOULD BE REWRITTEN USING FIXED-SIZW ARRAYS ONLY, NO EIGEN SO THAT 
	
	
	template<typename T, size_t N, size_t k>
	class CombinationWithoutRepetitionBase{
	
	
	public:
		CombinationWithoutRepetitionBase(){
			assert(N>=k);
		}
	};

	


	
	
	template<typename T, size_t N, size_t k>
	class CombinationWithoutRepetition : public CombinationWithoutRepetitionBase<T,N,k> {
	
		Eigen::Matrix<T,Binomial<N,k>::value,k> P;
		
		
	public:
		/////////////////////////////////////////////////////
		// constructor
		CombinationWithoutRepetition(Eigen::Matrix<T,1,N> V){
			// From Pascal Rule:
			// bin(N,k) = bin(n-1,k-1) + bin(n-1,k) 
			P <<	(Eigen::Matrix<T,k,Binomial<N-1,k-1>::value>() <<	Eigen::Matrix<T,1,Binomial<N-1,k-1>::value>().Constant(V(0)), 
																				CombinationWithoutRepetition<T,N-1,k-1>(Eigen::Matrix<T,N-1,1>(V.template segment<N-1>(1)))().transpose()).finished().transpose(),
					CombinationWithoutRepetition<T,N-1,k>(Eigen::Matrix<T,1,N-1>(V.template segment<N-1>(1)))();

		
		}
		
		/////////////////////////////////////////////////////
		// operator()
		Eigen::Matrix<T,Binomial<N,k>::value,k> operator()(){
			return P;
		}
	
	
	};
	
	/////////////////////////////////////////////////////
	/////////////////////////////////////////////////////
	// SPECIALIZATION FOR THE TERM BIN(N-1,k). The stop condition is reached at BIN(k,k)
	template<typename T, size_t k>
	class CombinationWithoutRepetition<T,k,k> : public CombinationWithoutRepetitionBase<T,k,k> {
		
		enum {N=k};
		const Eigen::Matrix<T,Binomial<N,k>::value,k> P;
		
		
	public:
		/////////////////////////////////////////////////////
		// constructor
		CombinationWithoutRepetition(Eigen::Matrix<T,1,N> V) : P(V){}
		
		/////////////////////////////////////////////////////
		// operator()
		Eigen::Matrix<T,Binomial<N,k>::value,k> operator()(){
			return P;
		}
		
	};
	
	
	/////////////////////////////////////////////////////
	/////////////////////////////////////////////////////
	// SPECIALIZATION FOR THE TERM BIN(N-1,k-1). The stop condition is reached at BIN(N,1)
	template<typename T, size_t N>
	class CombinationWithoutRepetition<T,N,1> : public CombinationWithoutRepetitionBase<T,N,1> {
		
		enum {k=1};
		const Eigen::Matrix<T,Binomial<N,k>::value,k> P;
		
		
	public:
		/////////////////////////////////////////////////////
		// constructor
		CombinationWithoutRepetition(Eigen::Matrix<T,1,N> V) : P(V.transpose()){}
		
		/////////////////////////////////////////////////////
		// operator()
		Eigen::Matrix<T,Binomial<N,k>::value,k> operator()(){
			return P;
		}
	};
	
	//////////////////////////////////////////////////////////////////////////
} // end namespace model

#endif


//////////////////////////////////////
//template <typename T, int N>
//Eigen::Matrix<T,N,1> cycle(Eigen::Matrix<T,N,1> const & V, int k = 1){
//	return	(k>0) ? cycle((Eigen::Matrix<T,N,1>() << Eigen::Matrix<T,N-1,1>(V.template segment<N-1>(1)), V(0)).finished(),k-1) : V;		
//}
