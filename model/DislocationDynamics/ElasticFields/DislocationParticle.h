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
#include <model/DislocationDynamics/ElasticFields/DislocationStress.h>
#include <model/DislocationDynamics/ElasticFields/DislocationEnergy.h>
#include <model/DislocationDynamics/ElasticFields/DislocationDisplacement.h>
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
        
#ifdef _MODEL_ENABLE_CELL_VERTEX_ALPHA_TENSORS_
        typedef std::array<MatrixDim,Pow<2,_dim>::value> CellAlphaType;
//        typedef typename CellPropertiesBase<VertexAlphaTensorType,VertexAlphaTensorType,VertexAlphaTensorType>::SpatialCellProperties SpatialCellProperties;
#else
        typedef MatrixDim CellAlphaType;
#endif
        typedef typename CellPropertiesBase<CellAlphaType,CellAlphaType,CellAlphaType>::SpatialCellProperties SpatialCellProperties;
        
//        static const CellAlphaType emptyAlpha;
        
        /**********************************************************************/
        static SpatialCellProperties init()
        {
#ifdef _MODEL_ENABLE_CELL_VERTEX_ALPHA_TENSORS_
            CellAlphaType temp;
            temp.fill(MatrixDim::Zero());
            return std::make_tuple(temp,
                                   temp,
                                   temp
                                   );
#else
            return std::make_tuple(MatrixDim::Zero(),
                                   MatrixDim::Zero(),
                                   MatrixDim::Zero()
                                   );
#endif
        }
        
    };
    
//    template<short unsigned int _dim>
//    struct TypeTraits<DislocationParticle<_dim>>
    
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
        
#ifdef userStressFile
#include userStressFile
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
        
        typedef typename TypeTraits<DislocationParticleType>::CellAlphaType CellAlphaType;
        //! The tangent vector of the parent DislocationSegment at this point
        
        const size_t&  sourceID;
        const size_t&  sinkID;
        const size_t   quadID;
        
        const VectorDimD  T;
        
        //! A const reference to the Burgers vector of the parent DislocationSegment
        const VectorDimD& B;
        
        //! A const reference to Quadrature abscissa corresponding to this particle
        const double quadAbscissa;
        //
        //! A const reference to Quadrature weight corresponding to this particle
        const double& quadWeight;
        
        /********************************************************/
        template <typename...EnabledFields>
        DislocationParticle(const VectorDimD& Pin,
                            const size_t& sourceID_in, const size_t& sinkID_in, const size_t& quadID_in,
                            const VectorDimD& Tin, const VectorDimD& Bin,
                            const double& qA,const double& qW,
                            //                            const EnabledFields&...enabled
                            const bool& enbStressSource,const bool& enbStressField,
                            const bool& enbDispSource,  const bool&,
                            const bool& enbEnrgSource,  const bool& enbEnrgField) :
        /* base init     */ PointSourceType(Pin,enbStressSource,enbDispSource,enbEnrgSource),
        /* base init     */ FieldPointType(enbStressField,enbEnrgField),
        //		/* base init     */  FieldPointType(Pin),
        //		/* base init     */  FieldPointType(PointSourceType::P),
        /* init list     */ sourceID(sourceID_in),
        /* init list     */ sinkID(sinkID_in),
        /* init list     */ quadID(quadID_in),
        /* init list     */ T(Tin),
        /* init list     */ B(Bin),
        /* init list     */ quadAbscissa(qA),
        /* init list     */ quadWeight(qW)
        {/*! Constructor updates the alpha-tensor of the cell containing *this
          */
            
#ifdef _MODEL_ENABLE_CELL_VERTEX_ALPHA_TENSORS_
            Eigen::Matrix<double,1,Pow<2,dim>::value> cellWeights=this->vertexWeigths();
            assert(std::fabs(cellWeights.sum()-1.0)<FLT_EPSILON);
            for(size_t v=0;v<Pow<2,dim>::value;++v)
            {
                if(static_cast<SingleSourcePoint<DislocationParticleType,StressField>* const>(this)->enabled)
                {
                    cellAlphaTensorStress()[v] += B * T.transpose() * quadWeight * cellWeights[v]; // add to alpha tensor of the cell
                }
                if(static_cast<SingleSourcePoint<DislocationParticleType,DisplacementField>* const>(this)->enabled)
                {
                    cellAlphaTensorDisp()[v] += B * T.transpose() * quadWeight * cellWeights[v]; // add to alpha tensor of the cell
                }
                if(static_cast<SingleSourcePoint<DislocationParticleType,ElasticEnergy>* const>(this)->enabled)
                {
                    cellAlphaTensorEnergy()[v] += B * T.transpose() * quadWeight * cellWeights[v]; // add to alpha tensor of the cell
                }
            }
            
#else
            if(static_cast<SingleSourcePoint<DislocationParticleType,StressField>* const>(this)->enabled)
            {
                cellAlphaTensorStress() += B * T.transpose() * quadWeight; // add to alpha tensor of the cell
            }
            if(static_cast<SingleSourcePoint<DislocationParticleType,DisplacementField>* const>(this)->enabled)
            {
                cellAlphaTensorDisp() += B * T.transpose() * quadWeight; // add to alpha tensor of the cell
            }
            if(static_cast<SingleSourcePoint<DislocationParticleType,ElasticEnergy>* const>(this)->enabled)
            {
                cellAlphaTensorEnergy() += B * T.transpose() * quadWeight; // add to alpha tensor of the cell
            }
#endif
            
        }
        
        
        
        /**********************************************************************/
        typename StressField::MatrixType stress() const
        {/*!\returns the total stress field on this DislocationQuadratureParticle
          */
            typename StressField::MatrixType temp(this->template field<StressField>());
            
            
#ifdef userStressFile
            temp+=userStress(*this);
#endif
            
            // stress is in fraction of mu. Keep only resolution of 8 digits
            //            const double dig(1.0e+08);
            //            temp = ((temp*dig).template cast<long int>().template cast<double>()/dig);
            
            //            return Material<Isotropic>::C2*(temp+temp.transpose());
            return temp;
            
        }
        
        /**********************************************************************/
        CellAlphaType& cellAlphaTensorStress()
        {
            return std::get<0>(*(this->pCell));
        }
        
        /**********************************************************************/
        const CellAlphaType& cellAlphaTensorStress() const
        {
            return std::get<0>(*(this->pCell));
        }
        
        /**********************************************************************/
        CellAlphaType& cellAlphaTensorDisp()
        {
            return std::get<1>(*(this->pCell));
        }
        
        /**********************************************************************/
        const CellAlphaType& cellAlphaTensorDisp() const
        {
            return std::get<1>(*(this->pCell));
        }
        
        /**********************************************************************/
        CellAlphaType& cellAlphaTensorEnergy()
        {
            return std::get<2>(*(this->pCell));
        }
        
        /**********************************************************************/
        const CellAlphaType& cellAlphaTensorEnergy() const
        {
            return std::get<2>(*(this->pCell));
        }
        
        /**********************************************************************/
        template <class T>
        friend T& operator << (T& os, const DislocationParticleType& p)
        {
//            typename StressField::MatrixType temp(p.template field<StressField>());
            //os  << p.source->sID<<"\t"<< ds.sink->sID<<"\t"
            os  << p.sourceID<< "\t"<< p.sinkID<<"\t"<<p.quadID<<"\t"
            /**/<< std::setprecision(15)<<std::scientific<<p.P.transpose()<<"\t"
            /**/<< std::setprecision(15)<<std::scientific<<p.B.transpose()<<"\t"
            /**/<< std::setprecision(15)<<std::scientific<<p.quadAbscissa<<"\t"
            /**/<< std::setprecision(15)<<std::scientific<<p.quadWeight<<"\t"
            /**/<< static_cast<const SingleSourcePoint<DislocationParticleType,StressField>* const>(&p)->enabled<<"\t"
            /**/<< static_cast<const FieldPointBase<DislocationParticleType,StressField>* const>(&p)->enabled<<"\t"
            /**/<< static_cast<const SingleSourcePoint<DislocationParticleType,DisplacementField>* const>(&p)->enabled<<"\t"
            //            /**/<< static_cast<const FieldPointBase<DislocationParticleType,DisplacementField>* const>(&p)->enabled<<"\t"
            /**/<< static_cast<const SingleSourcePoint<DislocationParticleType,ElasticEnergy>* const>(&p)->enabled<<"\t"
            /**/<< static_cast<const FieldPointBase<DislocationParticleType,ElasticEnergy>* const>(&p)->enabled;
//            /**/<<temp.row(0)<<" "<<temp.row(1)<<" "<<temp.row(2);
            return os;
        }
        
    };
    
    /**************************************************************************/
    /**************************************************************************/
}	// close namespace
#endif

