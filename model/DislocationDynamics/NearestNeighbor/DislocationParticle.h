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

#include <model/ParticleInteraction/PointSource.h>
#include <model/ParticleInteraction/FieldPoint.h>


namespace model {
	
	/********************************************************************************************/
	/********************************************************************************************/
	template<short unsigned int _dim>
	struct DislocationParticle :
    /* inheritance           */ public PointSource<DislocationParticle<_dim>,
    /*                                    */ _dim,
    /*                                    */ DislocationStress<_dim> >,
    /* inheritance           */ public FieldPoint<DislocationParticle<_dim>,
    /*                                    */ _dim,
    /*                                    */ DislocationStress<_dim> >
    {
        
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        
#ifdef UserStressFile
#include UserStressFile
#endif
        
        enum{dim=_dim};
		
		typedef DislocationParticle<_dim> DislocationParticleType;
		typedef DislocationStress<_dim>   StressField;

        typedef PointSource<DislocationParticleType,_dim,StressField> PointSourceType;
        typedef FieldPoint <DislocationParticleType,_dim,StressField> FieldPointType;
        typedef typename PointSourceType::VectorDimD VectorDimD;

        
		//! A const reference to Quadrature weight corresponding to this particle
        
		//const int k;			// *this is the k-th quadrature point on the segment
		const VectorDimD  P;
		const VectorDimD  T;
		
		//! A const reference to the Burgers vector of the parent DislocationSegment
		const VectorDimD B;
        
        const double quadAbscissa;
		const double quadWeight;
        
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
		/* base init     */  FieldPointType(Pin),
		/* init list     */ P(Pin),
		/* init list     */ T(Tin),
		/* init list     */ B(Bin),
		/* init list     */ quadAbscissa(qA),
		/* init list     */ quadWeight(qW)
//        /* init list     */ _stress(MatrixDim::Zero())
        {/*! Constructor updates the alpha-tensor of the cell containing this
          */
            
            //this->pCell->alpha += B * T.transpose() * quadWeight;
//            this->pCell->alpha += B * T.transpose() * quadWeight;
		}
        
        

        /********************************************************/
		typename StressField::MatrixType stress() const
        {/*! The total stress field on this DislocationQuadratureParticle
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
            
            
            typename StressField::MatrixType temp(this->template getField<StressField>());
            
            
            const double dig(1.0e+08);

            

            // stress is in fraction of mu. Keep only resolution of 6 digits
            temp = ((temp*dig).template cast<long int>().template cast<double>()/dig);
             
//            std::cout<<std::setprecision(15)<<std::scientific<<temp<<std::endl;


            
#ifdef UserStressFile
            //			temp+=userStress(k);
#endif
            
//			return _stress;
            return Material<Isotropic>::C2*(temp+temp.transpose());
//            return temp;

		}
//
        /********************************************************/
		typename StressField::MatrixType stress(const typename StressField::MatrixType& sourceBvpStress, const typename StressField::MatrixType& sinkBvpStress) const {
			/*! The total stress field on this DislocationQuadratureParticle
			 */
			return stress()+sourceBvpStress*(1.0-quadAbscissa)+sinkBvpStress*quadAbscissa;
		}
	
	};
    
    /**************************************************************************/
    /**************************************************************************/
}	// close namespace
#endif


//        void computeStress()
//        {
//
//             //_stress=MatrixDim::Zero();
//            for (typename CellMapType::const_iterator neighCellIter=this->neighborCellsBegin();neighCellIter!=this->neighborCellsEnd();++neighCellIter){
//                for (typename ParticleContainerType::const_iterator partIter=neighCellIter->second->particleBegin();partIter!=neighCellIter->second->particleEnd();++partIter){
//                    _stress+=(*partIter)->stress_at(P);
//                }
//            }
//        }
//




//
//		/********************************************************/
//        //		template<short unsigned int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule = GaussLegendre>
//		MatrixDim stress_at(const VectorDimD& Rfield) const
//        {/*! Returns the asymmetric (and dimensionless) part of the stress integrand generated by the current quadrature point.
//          * @param[in] k			the current quadrature point
//          * @param[in] Rfield	the vector connecting source point (corresponding to the current quadrature point) to field point
//          *
//          * The return value is calculated according to:
//          * Cai, W., Arsenlis, A., Weinberger, C., & Bulatov, V. (2006). A non-singular continuum theory of dislocations. Journal Of The Mechanics And Physics Of Solids, 54(3), 561–587.
//          *	\f[
//          *		d\mathbf{s} = (1-\nu) \left(1+\frac{3}{2}\frac{a^2}{R_a^2}\right)\frac{\partial \mathbf{r}}{\partial u}\otimes \left(\mathbf{b}\times \mathbf{R}\right)+
//          *		\mathbf{R}\otimes\left(\frac{\partial \mathbf{r}}{\partial u}\times\mathbf{b}\right)+
//          *		\frac{1}{2} \left[ \left(\mathbf{R}\times\mathbf{b}\right)\cdot \frac{\partial \mathbf{r}}{\partial u} \right]\left[\mathbf{I}\left(1+\frac{3a^2}{R_a^2}\right)+\frac{3}{R_a^2}\mathbf{R}\otimes\mathbf{R}\right]
//          *	\f]
//          *  where \f$R_a^2=|\mathbf{R}|^2+a^2\f$ is the modified squared norm of \f$\mathbf{R}\f$.
//          *
//          *	The return value is asymmetric and dimensionless in the sense that the actual stress integrand is:
//          *	\f[
//          *		d\mathbf{\sigma}=\frac{\mu}{4\pi (1-\nu)}\left(d\mathbf{s}+d\mathbf{s}^T\right)
//          *	\f]
//          */
//			VectorDimD R(Rfield-P);
//			double RaSquared (R.squaredNorm() + a2);
//			MatrixDim temp(   (Material<Isotropic>::C1*(1.0+1.5*a2/RaSquared)*T*(B.cross(R)).transpose()
//                               + 	R*(T.cross(B)).transpose()
//                               +   0.5* R.cross(B).dot(T) * (I*(1.0+3.0*a2/RaSquared) + 3.0/RaSquared*R*R.transpose())
//                               )/std::pow(sqrt(RaSquared),3)*quadWeight);
//
////			return Material<Isotropic>::C2*(temp+temp.transpose());
//			return temp;
//
//		}
////
////
////
////
////		/********************************************************/
////		MatrixDim stress() const
////        {/*! The total stress field on this DislocationParticle
////			 */
////            //            MatrixDim temp(MatrixDim::Zero());
////            //            for (typename CellMapType::const_iterator nearCellIter=this->nearCellsBegin();nearCellIter!=this->nearCellsEnd();++nearCellIter){
////            //                for (typename ParticleContainerType::const_iterator partIter=nearCellIter->second->particleBegin();partIter!=nearCellIter->second->particleEnd();++partIter){
////            //                    temp+=(*partIter)->stress_at(P);
////            //                }
////            //            }
////            //            if(useMultipoleStress){
////            //                for (typename CellMapType::const_iterator farCellIter=this->farCellsBegin();farCellIter!=this->farCellsEnd();++farCellIter){
////            //                    temp+= farCellIter->second->multipoleStress(P);
////            //                }
////            //            }
////
////            MatrixDim temp(MatrixDim::Zero());
////            for (typename CellMapType::const_iterator neighCellIter=this->neighborCellsBegin();neighCellIter!=this->neighborCellsEnd();++neighCellIter){
////                for (typename ParticleContainerType::const_iterator partIter=neighCellIter->second->particleBegin();partIter!=neighCellIter->second->particleEnd();++partIter){
////                    temp+=(*partIter)->stress_at(P);
////                }
////            }
////
////
////            switch (nearCellStressApproximation)
////            {
////                case FULL: // quadrature-quadrature
////                    // finish here
////                    assert(0 && "FINISH HERE");
////                    break;
////                case CELL_PARTICLE: // cell-quadrature
////                    for (typename CellMapType::const_iterator nearCellIter=this->nearCellsBegin();nearCellIter!=this->nearCellsEnd();++nearCellIter){
////                        temp+= nearCellIter->second->multipoleStress(P);
////                    }
////                    break;
////                case CELL_CELL: // cell-cell
////                    //temp+= this->pCell->nearStress;
////                    temp+= this->pCell->nearStress;
////                    break;
////
////                default: // no computation
////                    break;
////            }
////
////
////            switch (farCellStressApproximation)
////            {
////                case FULL: // quadrature-quadrature
////                    // finish here
////                    assert(0 && "FINISH HERE");
////                    break;
////                case CELL_PARTICLE: // cell-quadrature
////                    for (typename CellMapType::const_iterator farCellIter=this->farCellsBegin();farCellIter!=this->farCellsEnd();++farCellIter){
////                        temp+= farCellIter->second->multipoleStress(P);
////                    }
////                    break;
////                case CELL_CELL: // cell-cell
////                    //temp+= this->pCell->farStress;
////                    temp+= this->pCell->farStress;
////                    break;
////
////                default: // no computation
////                    break;
////            }
////
////
////#ifdef UserStressFile
////            //			temp+=userStress(k);
////#endif
////
////			return temp;
////		}
////
////		/********************************************************/
////		MatrixDim stress(const MatrixDim& sourceBvpStress, const MatrixDim& sinkBvpStress) const {
////			/*! The total stress field on this DislocationParticle
////			 */
////            MatrixDim temp(stress());
////			temp+=sourceBvpStress*(1.0-quadAbscissa)+sinkBvpStress*quadAbscissa;
////			return temp;
////		}
//
