/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DISLOCATIONCELL_H_
#define model_DISLOCATIONCELL_H_

#include <Eigen/Dense>
#include <model/SpaceDecomposition/SpaceCell.h>
#include <model/Dislocations/Materials/Material.h>
#include <model/Dislocations/NearestNeighbor/MultipoleExpansion.h>
#include <model/Utilities/TypeTraits.h>

namespace model {	
	
    double cellSize= 1000.0;
    
    /**************************************************************************/
	/* TRAITS *****************************************************************/
	template<short unsigned int dim, double & cellSize>
	struct DislocationQuadratureParticle;

	template<short unsigned int dim, double & cellSize>
	struct DislocationCell;
    
    template<>
	struct TypeTraits<DislocationQuadratureParticle<3,cellSize> > {
        //enum {dim=3;}
        typedef DislocationQuadratureParticle<3,cellSize> ParticleType;
        typedef               DislocationCell<3,cellSize> CellType;
    };

    template<>
	struct TypeTraits<DislocationCell<3,cellSize> > {
        typedef DislocationQuadratureParticle<3,cellSize> ParticleType;
        typedef               DislocationCell<3,cellSize> CellType;
    };
    
    
	/********************************************************************************************/
	/********************************************************************************************/
	/*! \brief A dim-dimensional cell occupying the spatial region cellID<= x/cellSize < (cellID+1). 
	 *  DislocationCell is aware off all ParticleType objects present inside it. 
	 */
	template<short unsigned int dim, double & cellSize>
	struct DislocationCell : public SpaceCell<DislocationCell<dim,cellSize>,dim,cellSize>{
		
        typedef typename TypeTraits<DislocationCell<dim,cellSize> >::ParticleType ParticleType;
		typedef SpaceCell<DislocationCell<dim,cellSize>,dim,cellSize> Base;
        typedef typename Base::VectorDimI VectorDimI;
        typedef typename Base::VectorDimD VectorDimD;
        
        typedef Eigen::Matrix<double,dim,dim>  MatrixDimD;	// remove this with Dislocation Stuff
        MatrixDimD alpha;	// the dislocation density tensor !! remove this with Dislocation Stuff
        
	public:
        
			
		/* Constructor *******************************************/
		DislocationCell(const VectorDimI& cellID_in) : Base::SpaceCell(cellID_in), alpha(MatrixDimD::Zero()){
        }
		
		
        /* multipoleStress ****************************************/
        MatrixDimD multipoleStress(const VectorDimD& Rfield) const {
            return MultipoleExpansion<dim>::multipoleStress(Rfield,alpha,this->cellID,cellSize);
        }

	};
    
  	////////////////////////////////////////////////////////////////////////////////
}	// close namespace model
#endif

