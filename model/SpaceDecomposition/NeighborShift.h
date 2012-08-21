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
	
    template <int dim, int neighborOrder>
    struct NeighborShiftTraits{
        enum{Nneighbors=Pow<2*neighborOrder+1,dim>::value};
        typedef Eigen::Matrix<int,dim,Nneighbors> MatrixType; // WORK WITH TRANSPOSE
    };
    
    /**************************************************************************/
    /**************************************************************************/
    template <int dim, int neighborOrder>
    struct NeighborShift{
        static_assert(dim>=1, "dim MUST BE >=1.");
        static_assert(neighborOrder>=0, "neighborOrder MUST BE >=1.");
        
        enum{Nneighbors=NeighborShiftTraits<dim,neighborOrder>::Nneighbors};
        static const typename NeighborShiftTraits<dim,neighborOrder>::MatrixType shifts;
        
        /**************************************************************************/
        static typename NeighborShiftTraits<dim,neighborOrder>::MatrixType getShifts(){
            const typename NeighborShiftTraits<dim-1,neighborOrder>::MatrixType shiftsLower(NeighborShift<dim-1,neighborOrder>::getShifts());
            const typename NeighborShiftTraits<1,neighborOrder>::MatrixType shifts1(NeighborShift<1,neighborOrder>::getShifts());
            typename NeighborShiftTraits<dim,neighborOrder>::MatrixType temp;
            for (int k=0;k<NeighborShift<1,neighborOrder>::Nneighbors;++k){
                temp.template block<dim-1,NeighborShiftTraits<dim-1,neighborOrder>::Nneighbors>(0,k*NeighborShiftTraits<dim-1,neighborOrder>::Nneighbors)=shiftsLower;
                temp.template block<1,NeighborShiftTraits<dim-1,neighborOrder>::Nneighbors>(dim-1,k*NeighborShiftTraits<dim-1,neighborOrder>::Nneighbors).setConstant(shifts1(k));
            }
            return temp;
        }
        
        /**************************************************************************/
        static typename NeighborShiftTraits<dim,neighborOrder>::MatrixType neighborIDs(const Eigen::Matrix<int,dim,1>& cellID){
            typename NeighborShiftTraits<dim,neighborOrder>::MatrixType temp(NeighborShiftTraits<dim,neighborOrder>::MatrixType::Zero());
            for (int n=0;n<Nneighbors;++n){
                temp.col(n)=cellID+shifts.col(n);
            }
            return temp;
        }
        
    };
    
    /**************************************************************************/
    /**************************************************************************/
    template <int neighborOrder>
    struct NeighborShift<1,neighborOrder>{
        enum{dim=1};
        enum{Nneighbors=NeighborShiftTraits<dim,neighborOrder>::Nneighbors};
        static const typename NeighborShiftTraits<dim,neighborOrder>::MatrixType shifts;
        
        
        /**************************************************************************/
        static typename NeighborShiftTraits<dim,neighborOrder>::MatrixType getShifts(){
            return Eigen::Matrix<int,NeighborShiftTraits<dim,neighborOrder>::Nneighbors,dim>::LinSpaced(-neighborOrder,neighborOrder);
        }
        
        /**************************************************************************/
        static typename NeighborShiftTraits<dim,neighborOrder>::MatrixType neighborIDs(const Eigen::Matrix<int,dim,1>& cellID){
            typename NeighborShiftTraits<dim,neighborOrder>::MatrixType temp(NeighborShiftTraits<dim,neighborOrder>::MatrixType::Zero());
            for (int n=0;n<Nneighbors;++n){
                temp.col(n)=cellID+shifts.col(n);
            }
            return temp;
        }
        
    };
    
    // static data members
    template <int dim, int neighborOrder>
    const typename NeighborShiftTraits<dim,neighborOrder>::MatrixType NeighborShift<dim,neighborOrder>::shifts=NeighborShift<dim,neighborOrder>::getShifts();
    

	
	////////////////////////////////////////////////////////////////////////////////
}	// close namespace model
#endif

