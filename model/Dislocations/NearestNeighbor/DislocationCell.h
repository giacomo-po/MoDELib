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


namespace model {	
	
	/********************************************************************************************/
	/********************************************************************************************/
	/*! \brief A dim-dimensional cell occupying the spatial region cellID<= x/cellSize < (cellID+1). 
	 *  DislocationCell is aware off all ParticleType objects present inside it. 
	 */
	template<typename ParticleType, short unsigned int dim, double & cellSize>
	struct DislocationCell : public SpaceCell<ParticleType,dim,cellSize>{
		
		typedef SpaceCell<ParticleType,dim,cellSize> SpaceCellType;

        typedef Eigen::Matrix<double,dim,dim>  MatrixDimD;	// remove this with Dislocation Stuff
        MatrixDimD alpha;	// the dislocation density tensor !! remove this with Dislocation Stuff
        
	public:
        
			
		/* Constructor *******************************************/
		DislocationCell(const VectorDimI& cellID_in) : SpaceCellType::SpaceCell(const VectorDimI& cellID_in), alpha(MatrixDimD::Zero()){}
		

		
		/* addParticle *******************************************/
		void addParticle(const ParticleType* const pP){
            SpaceCellType::addParticle(pP);
            alpha+=pP->B * pP->T.transpose() * pP->quadWeight;
		}
		

		
//		/* removeParticle ****************************************/
//		void removeParticle(const ParticleType* const pP){
//            SpaceCellType::removeParticle(pP);
//            alpha-=pP->B * pP->T.transpose() * pP->quadWeight;
//
//		}
		
        /* multipoleStress ****************************************/
        MatrixDimD multipoleStress(const VectorDimD& Rfield) const {
            return MultipoleExpansion<dim>::multipoleStress(Rfield,alpha);
        }

        
       		
	};
    
  	////////////////////////////////////////////////////////////////////////////////
}	// close namespace model
#endif

