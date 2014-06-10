/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _MODEL_DislocationParticle_H_
#define _MODEL_DislocationParticle_H_

#include <Eigen/Dense>

#include <model/SpaceDecomposition/SpatialCellParticle.h>
#include <model/DislocationDynamics/Materials/Material.h>
#include <model/DislocationDynamics/NearestNeighbor/DislocationStress.h>
#include <model/DislocationDynamics/NearestNeighbor/DislocationEnergy.h>
#include <model/ParticleInteraction/PointSource.h>
#include <model/ParticleInteraction/FieldPoint.h>
#include <model/Mesh/SimplicialMesh.h>



namespace model {
	
	/********************************************************************************************/
	/********************************************************************************************/
	template<short unsigned int _dim>
	struct DislocationParticle :
    /* inheritance           */ public PointSource<DislocationParticle<_dim>,
    /*                                    */ _dim,
    /*                                    */ DislocationStress<_dim>,
    /*                                    */ DislocationEnergy<_dim> >,
    /* inheritance           */ public FieldPoint<DislocationParticle<_dim>,
    /*                                    */ _dim,
    /*                                    */ DislocationStress<_dim>,
    /*                                    */ DislocationEnergy<_dim> >
    {
        
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        
#ifdef UserStressFile
#include UserStressFile
#endif
        
        enum{dim=_dim};
		
		typedef DislocationParticle<_dim> DislocationParticleType;
		typedef DislocationStress<_dim>   StressField;
		typedef DislocationEnergy<_dim>   ElasticEnergy;

        typedef PointSource<DislocationParticleType,_dim,StressField,ElasticEnergy> PointSourceType;
        typedef FieldPoint <DislocationParticleType,_dim,StressField,ElasticEnergy> FieldPointType;
        typedef typename PointSourceType::VectorDimD VectorDimD;

        
        
		//! The tangent vector of the parent DislocationSegment at this point
		const VectorDimD  T;
		
		//! A const reference to the Burgers vector of the parent DislocationSegment
		const VectorDimD& B;
        
        //! A const reference to Quadrature abscissa corresponding to this particle
        const double& quadAbscissa;
        
        //! A const reference to Quadrature weight corresponding to this particle
		const double& quadWeight;
        
//        enum{FULL=0,CELL_PARTICLE=1,CELL_CELL=2};
//        
//        static int nearCellStressApproximation;
//        static int  farCellStressApproximation;
//
//        MatrixDim _stress;
		
		/********************************************************/
		DislocationParticle(const VectorDimD& Pin, const VectorDimD& Tin, const VectorDimD& Bin,
                            const double& qA,const double& qW) :
		/* base init     */ PointSourceType(Pin),
//		/* base init     */  FieldPointType(Pin),
//		/* base init     */  FieldPointType(PointSourceType::P),
		/* init list     */ T(Tin),
		/* init list     */ B(Bin),
		/* init list     */ quadAbscissa(qA),
		/* init list     */ quadWeight(qW)
        {/*! Constructor updates the alpha-tensor of the cell containing this
          */
            
            //this->pCell->alpha += B * T.transpose() * quadWeight;
//            this->pCell->alpha += B * T.transpose() * quadWeight;
		}
        
        

        /********************************************************/
		typename StressField::MatrixType stress() const
        {/*!\returns the total stress field on this DislocationQuadratureParticle
          */

            
            
//            switch (nearCellStressApproximation)
//            {
//                case FULL: // quadrature-quadrature
//                    // finish here
//                    assert(0 && "FINISH HERE");
//                    break;
//                case CELL_PARTICLE: // cell-quadrature
//                    for (typename CellMapType::const_iterator nearCellIter=this->nearCellsBegin();nearCellIter!=this->nearCellsEnd();++nearCellIter){
//                        temp+= nearCellIter->second->multipoleStress(P);
//                    }
//                    break;
//                case CELL_CELL: // cell-cell
//                    //temp+= this->pCell->nearStress;
//                    temp+= this->pCell->nearStress;
//                    break;
//                    
//                default: // no computation
//                    break;
//            }
            
            
//            switch (farCellStressApproximation)
//            {
//                case FULL: // quadrature-quadrature
//                    // finish here
//                    assert(0 && "FINISH HERE");
//                    break;
//                case CELL_PARTICLE: // cell-quadrature
//                    for (typename CellMapType::const_iterator farCellIter=this->farCellsBegin();farCellIter!=this->farCellsEnd();++farCellIter){
//                        temp+= farCellIter->second->multipoleStress(P);
//                    }
//                    break;
//                case CELL_CELL: // cell-cell
//                    //temp+= this->pCell->farStress;
//                    temp+= this->pCell->farStress;
//                    break;
//                    
//                default: // no computation
//                    break;
//            }
            
            
            typename StressField::MatrixType temp(this->template field<StressField>());
            
            // stress is in fraction of mu. Keep only resolution of 8 digits
            const double dig(1.0e+08);
            temp = ((temp*dig).template cast<long int>().template cast<double>()/dig);
 
#ifdef UserStressFile
            //			temp+=userStress(k);
#endif
            
            return Material<Isotropic>::C2*(temp+temp.transpose());

		}


	};
    
    /**************************************************************************/
    /**************************************************************************/
}	// close namespace
#endif
