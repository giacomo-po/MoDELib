/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _MODEL_DislocationParticle_H_
#define _MODEL_DislocationParticle_H_

#include <tuple>
#include <Eigen/Dense>

#include <model/SpaceDecomposition/CellPropertiesBase.h>
#include <model/DislocationDynamics/Materials/Material.h>
#include <model/DislocationDynamics/NearestNeighbor/DislocationStress.h>
#include <model/DislocationDynamics/NearestNeighbor/DislocationEnergy.h>
#include <model/DislocationDynamics/NearestNeighbor/DislocationDisplacement.h>
#include <model/ParticleInteraction/PointSource.h>
#include <model/ParticleInteraction/FieldPoint.h>
#include <model/Mesh/SimplicialMesh.h>
#include <model/Utilities/TypeTraits.h>


namespace model
{
    
    /**************************************************************************/
	/**************************************************************************/
    template<short unsigned int _dim>
	struct DislocationParticle;
    
    /*!\brief Class template specialization which defines the properties (data
     * members) of SpatialCell(s) containing DislocationParticle(s).
     *
     * The properties are:
     * - the alpha tensor (MatrixDim)
     */
    template<short unsigned int _dim>
	struct TypeTraits<DislocationParticle<_dim>>
    {
    
        typedef Eigen::Matrix<double,_dim,_dim> MatrixDim;
        typedef typename CellPropertiesBase<MatrixDim>::SpatialCellProperties SpatialCellProperties;
        
        /**********************************************************************/
        static SpatialCellProperties init()
        {
            return std::make_tuple(MatrixDim::Zero());
        }
        
    };
	
	/**************************************************************************/
	/**************************************************************************/
	template<short unsigned int _dim>
	struct DislocationParticle :
    /* inheritance           */ public PointSource<DislocationParticle<_dim>,
    /*                                    */ _dim,
    /*                                    */ DislocationStress<_dim>,
    /*                                    */ DislocationDisplacement<_dim>,
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
		
		typedef DislocationParticle<_dim>       DislocationParticleType;
		typedef DislocationStress<_dim>         StressField;
		typedef DislocationDisplacement<_dim>   DisplacementField;
		typedef DislocationEnergy<_dim>         ElasticEnergy;
        
        
        typedef Eigen::Matrix<double,_dim,_dim> MatrixDim;
        typedef CellPropertiesBase<MatrixDim> CellPropertiesBaseType;
        typedef typename CellPropertiesBaseType::SpatialCellProperties SpatialCellProperties;
        
        typedef PointSource<DislocationParticleType,_dim,StressField,DisplacementField,ElasticEnergy> PointSourceType;
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
        {/*! Constructor updates the alpha-tensor of the cell containing *this
          */
            
            cellAlphaTensor() += B * T.transpose() * quadWeight; // add to alpha tensor of the cell
		}
        
        
        
        /**********************************************************************/
		typename StressField::MatrixType stress() const
        {/*!\returns the total stress field on this DislocationQuadratureParticle
          */
            typename StressField::MatrixType temp(this->template field<StressField>());
            
            // stress is in fraction of mu. Keep only resolution of 8 digits
            const double dig(1.0e+08);
            temp = ((temp*dig).template cast<long int>().template cast<double>()/dig);
            
#ifdef UserStressFile
            			temp+=userStress(*this);
#endif
            
//            return Material<Isotropic>::C2*(temp+temp.transpose());
            return temp;
            
		}
        
        /**********************************************************************/
        Eigen::Matrix<double,_dim,_dim>& cellAlphaTensor()
        {
            return std::get<0>(*(this->pCell));
        }
        
        /**********************************************************************/
        const Eigen::Matrix<double,_dim,_dim>& cellAlphaTensor() const
        {
            return std::get<0>(*(this->pCell));
        }
        
        
	};
    
    /**************************************************************************/
    /**************************************************************************/
}	// close namespace
#endif



//        enum{FULL=0,CELL_PARTICLE=1,CELL_CELL=2};
//
//        static int nearCellStressApproximation;
//        static int  farCellStressApproximation;
//
//        MatrixDim _stress;




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
