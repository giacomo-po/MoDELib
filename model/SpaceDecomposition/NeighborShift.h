/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_NEIGHBORSHIFT_H_
#define model_NEIGHBORSHIFT_H_

#include <Eigen/Dense>
#include <model/Math/CompileTimeMath/Pow.h>

namespace model {
	
	
	template <short unsigned int dim>
	struct NeighborShift{
		//enum {poolSize=3};
		enum {Nneighbors= Pow<3,dim>::value};
		static const Eigen::Matrix<int,dim, Pow<3,dim>::value> shifts;
	};
	

	
	template <>
	struct NeighborShift<3>{
		enum {dim=3};
		//enum {poolSize=3};
		enum {Nneighbors=  Pow<3,dim>::value};
		// First column changes every 3^2=9 values
		// Second column changes every 3^1 values
		// Third column changes every 3^0 values
		
		// In general the m-th column changes every N^(k-m) values m=0...k
		
		/*
		 for (int i=0;i<Pow<N,K>::value;++i){
		 
		 for (j=0;j<K;++j){
		 
		 shifts(i,j)= pool();
		 
		 }
		 }
		 
		 */
		
		
		static Eigen::Matrix<int,dim, Nneighbors> neighborIDs(const Eigen::Matrix<int,dim,1>& cellID){
			Eigen::Matrix<int,dim, Nneighbors> temp(Eigen::Matrix<int,dim, Nneighbors>::Zero());
			for (int n=0;n<Nneighbors;++n){
				temp.col(n)=cellID+shifts.col(n);
			}
			return temp;
		}
		
		
		static Eigen::Matrix<int,dim, Nneighbors> getShifts(){
			return (Eigen::Matrix<int, Nneighbors,dim>()<< 0, 0, 0,
					/*                                                  */ 0, 0, 1,
					/*                                                  */ 0, 0,-1,
					/*                                                  */ 0, 1, 0,
					/*                                                  */ 0, 1, 1,
					/*                                                  */ 0, 1,-1,
					/*                                                  */ 0,-1, 0,
					/*                                                  */ 0,-1, 1,
					/*                                                  */ 0,-1,-1,
					
					/*                                                  */ 1, 0, 0,
					/*                                                  */ 1, 0, 1,
					/*                                                  */ 1, 0,-1,
					/*                                                  */ 1, 1, 0,
					/*                                                  */ 1, 1, 1,
					/*                                                  */ 1, 1,-1,
					/*                                                  */ 1,-1, 0,
					/*                                                  */ 1,-1, 1,
					/*                                                  */ 1,-1,-1,
					
					/*                                                  */-1, 0, 0,
					/*                                                  */-1, 0, 1,
					/*                                                  */-1, 0,-1,
					/*                                                  */-1, 1, 0,
					/*                                                  */-1, 1, 1,
					/*                                                  */-1, 1,-1,
					/*                                                  */-1,-1, 0,
					/*                                                  */-1,-1, 1,
					/*                                                  */-1,-1,-1).finished().transpose();
		}
		
		
		
		//	enum {dim=3};
		//enum {poolSize=3};
		static const Eigen::Matrix<int,dim, Pow<3,dim>::value> shifts;
		
	};
	

	//	declare static data members	
	const Eigen::Matrix<int,3, Pow<3,3>::value> NeighborShift<3>::shifts(NeighborShift<3>::getShifts());
	
	////////////////////////////////////////////////////////////////////////////////
}	// close namespace model
#endif

