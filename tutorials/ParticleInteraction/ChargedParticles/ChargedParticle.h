/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2012 by Giacomo Po <gpo@ucla.edu>
 * Copyright (C) 2012 by Shao-Ching Huang <sch@ucla.edu>
 * Copyright (C) 2012 by Tajendra Singh <tvsingh@ucla.edu>
 * Copyright (C) 2012 by Tamer Crosby <tcrosby@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _ChargedParticle_h
#define _ChargedParticle_h

#include <model/SpaceDecomposition/CellPropertiesBase.h>
#include <model/ParticleInteraction/PointSource.h>
#include <model/ParticleInteraction/FieldPoint.h>

#include <tutorials/ParticleInteraction/ChargedParticles/ElectricField.h>
#include <tutorials/ParticleInteraction/ChargedParticles/MagneticField.h>

namespace model
{
    
    template <short unsigned int _dim>
    class ChargedParticle;
    
    template<short unsigned int _dim>
    struct TypeTraits<ChargedParticle<_dim>>
    {
        
        typedef Eigen::Matrix<double,_dim,_dim> MatrixDim;
        typedef typename CellPropertiesBase<MatrixDim>::SpatialCellProperties SpatialCellProperties;
        
        /**********************************************************************/
        static SpatialCellProperties init()
        {
            return std::make_tuple(MatrixDim::Zero());
        }
        
    };
    
    template <short unsigned int _dim>
    class ChargedParticle :
    /* inheritance     */ public PointSource<ChargedParticle<_dim>,
    /*                                    */ _dim,
    /*                                    */ ElectricField<_dim>,
    /*                                    */ MagneticField<_dim> >,
    /* inheritance     */ public  FieldPoint<ChargedParticle<_dim>,
    /*                                    */ _dim,
    /*                                    */ ElectricField<_dim>,
    /*                                    */ MagneticField<_dim> >
    {
    public:
        
        enum{dim=_dim};
        
        typedef ElectricField<_dim> Efield;
        typedef MagneticField<_dim> Bfield;

        
        typedef PointSource<ChargedParticle<_dim>,_dim,Efield,Bfield> PointSourceType;
        typedef FieldPoint <ChargedParticle<_dim>,_dim,Efield,Bfield> FieldPointType;

        typedef typename PointSourceType::VectorDimD VectorDimD;

    private:

        VectorDimD _p; // THIS SHOULD BE STORED IN SpatialCellParticle
        VectorDimD _v; // THIS SHOULD BE STORED IN SpatialCellParticle

    public:

        const double q; // the electric charge of the particle
        
        /*****************************************/
        ChargedParticle(const VectorDimD& pIN, const VectorDimD& vIN, const double& qIN) :
        /* init list */ PointSourceType(pIN), //  SpatialCellParticle must be constructed with initial position
//        /* init list */  FieldPointType(pIN), //  SpatialCellParticle must be constructed with initial position
//        /* init list */ _p(pIN),
        /* init list */ _v(vIN),
        /* init list */ q(qIN)//,
        {/*!@param[in] pIN position of this ChargedParticle
          * @param[in] qIN charge of this ChargedParticle
          * 
          * Constructor with input position and charge
          */
        }
                
//        /*****************************************/
//        const VectorDimD& P() const
//        {/*! The position of this ChargedParticle
//          */
//            return _p;
//        }
        
        /*****************************************/
        const VectorDimD& V() const
        {/*! The velocity of this ChargedParticle
          */
            return _v;
        }
        
        /*****************************************/
        VectorDimD force() const
        {/*! The charge of this ChargedParticle
          */
            return q*this->template field<Efield>();
        }
        
        /*****************************************/
        template <class T>
        friend T& operator << (T& os, const ChargedParticle& cP)
        {/*! operator << is used to output ChargedParticle info
          *  Example:
          *  ChargedParticle p;
          *  std::cout<<p;
          */
            os<<cP.sID<<"\t"
            <<cP.P.transpose()<<"\t"
            <<cP.q<<"\t"
            <<cP.force().transpose()<<"\t"
#ifdef _MODEL_MPI_
            <<cP.mpiID()
#else
            <<cP.sID      
#endif
            <<std::endl;
            return os;
        }
            
    };
            
}
#endif
