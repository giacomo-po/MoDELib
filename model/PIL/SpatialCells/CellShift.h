/* This file is part of PIL, the Particle Interaction Library.
 *
 * Copyright (C) 2012 by Giacomo Po <gpo@ucla.edu>
 * Copyright (C) 2012 by Shao-Ching Huang <sch@ucla.edu>
 * Copyright (C) 2012 by Tajendra Singh <tvsingh@ucla.edu>
 * Copyright (C) 2012 by Tamer Crosby <tcrosby@ucla.edu>
 *
 * PIL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef pil_CellShift_H_
#define pil_CellShift_H_

#include <Eigen/Dense>
#include <model/Math/CompileTimeMath/Pow.h>

namespace pil {
	
    template <int dim, int neighborOrder>
    struct CellShiftTraits{
        /*! A trait class defining how many neighbors of order neighborOrder a 
         *  SpatialCell has in dim dimensions.
         */
        
        //! The number of "neighbors" that a spatial cell has in dim dimensions
        enum{Nneighbors=model::Pow<2*neighborOrder+1,dim>::value};
        
        //! The type of matrix that contains the cellID(s) of all the neighbor cells
        typedef Eigen::Matrix<int,dim,Nneighbors> MatrixType; // WORK WITH TRANSPOSE
    };
    
    /**************************************************************************/
    /**************************************************************************/
    template <int dim, int neighborOrder>
    struct CellShift{
        /*! A class to compute the neighbors (or shifts) of order neighborOrder 
         *  for a SpatialCell in dim dimensions.
         */
        
        static_assert(dim>=1, "dim MUST BE >=1.");
        static_assert(neighborOrder>=0, "neighborOrder MUST BE >=0.");
        
        enum{Nneighbors=CellShiftTraits<dim,neighborOrder>::Nneighbors};
        static const typename CellShiftTraits<dim,neighborOrder>::MatrixType shifts;
        
        /**************************************************************************/
        static typename CellShiftTraits<dim,neighborOrder>::MatrixType getShifts(){
            /*! 
             */
            const typename CellShiftTraits<dim-1,neighborOrder>::MatrixType shiftsLower(CellShift<dim-1,neighborOrder>::getShifts());
            const typename CellShiftTraits<1,neighborOrder>::MatrixType shifts1(CellShift<1,neighborOrder>::getShifts());
            typename CellShiftTraits<dim,neighborOrder>::MatrixType temp;
            for (int k=0;k<CellShift<1,neighborOrder>::Nneighbors;++k){
                temp.template block<dim-1,CellShiftTraits<dim-1,neighborOrder>::Nneighbors>(0,k*CellShiftTraits<dim-1,neighborOrder>::Nneighbors)=shiftsLower;
                temp.template block<1,CellShiftTraits<dim-1,neighborOrder>::Nneighbors>(dim-1,k*CellShiftTraits<dim-1,neighborOrder>::Nneighbors).setConstant(shifts1(k));
            }
            return temp;
        }
        
        /**************************************************************************/
        static typename CellShiftTraits<dim,neighborOrder>::MatrixType neighborIDs(const Eigen::Matrix<int,dim,1>& cellID){
            typename CellShiftTraits<dim,neighborOrder>::MatrixType temp(CellShiftTraits<dim,neighborOrder>::MatrixType::Zero());
            for (int n=0;n<Nneighbors;++n){
                temp.col(n)=cellID+shifts.col(n);
            }
            return temp;
        }
        
    };
    
    /**************************************************************************/
    /**************************************************************************/
    template <int neighborOrder>
    struct CellShift<1,neighborOrder>{
        enum{dim=1};
        enum{Nneighbors=CellShiftTraits<dim,neighborOrder>::Nneighbors};
        static const typename CellShiftTraits<dim,neighborOrder>::MatrixType shifts;
        
        
        /**************************************************************************/
        static typename CellShiftTraits<dim,neighborOrder>::MatrixType getShifts(){
            return Eigen::Matrix<int,CellShiftTraits<dim,neighborOrder>::Nneighbors,dim>::LinSpaced(-neighborOrder,neighborOrder);
        }
        
        /**************************************************************************/
        static typename CellShiftTraits<dim,neighborOrder>::MatrixType neighborIDs(const Eigen::Matrix<int,dim,1>& cellID){
            typename CellShiftTraits<dim,neighborOrder>::MatrixType temp(CellShiftTraits<dim,neighborOrder>::MatrixType::Zero());
            for (int n=0;n<Nneighbors;++n){
                temp.col(n)=cellID+shifts.col(n);
            }
            return temp;
        }
        
    };
    
    // static data members
    template <int dim, int neighborOrder>
    const typename CellShiftTraits<dim,neighborOrder>::MatrixType CellShift<dim,neighborOrder>::shifts=CellShift<dim,neighborOrder>::getShifts();
    

	
	////////////////////////////////////////////////////////////////////////////////
}	// close namespace pil
#endif

