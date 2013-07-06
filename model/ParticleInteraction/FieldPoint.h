/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2012 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _model_FieldPoint_h
#define _model_FieldPoint_h

//#include <model/Utilities/CRTP.h>
#include <model/SpaceDecomposition/SpatialCellObserver.h>
#include <model/ParticleInteraction/FieldPointBase.h>


namespace model {
    
    /*!\brief A class template that contains the underlying FieldType data of
     * the Derived particle.
     *
     * If MPI is not used, the data is stored locally in the based class
     * Eigen::Matrix<scalar,rows,cols>.
     *
     * In case that MPI is used, the underlying data is not stored with the 
     * object but is deallocated to a static vector that can efficiently be 
     * syncronized using MPI_Allgatherv, and conveniently used as an Eigen 
     * object via the Eigen::Map class.
     */
    template<typename Derived,
    /*    */ short unsigned int _dim,
    /*    */ typename ...FieldTypes>
    struct FieldPoint;
    
    /**************************************************************************/
    /**************************************************************************/
    // template specialization: one or more fields
    template<typename Derived,
    /*    */ short unsigned int _dim,
    /*    */ typename FieldType,
    /*    */ typename ...MoreFieldTypes>
    struct FieldPoint<Derived,_dim,FieldType,MoreFieldTypes...> :
    /* inheritance  */ public FieldPoint<Derived,_dim,MoreFieldTypes...>,
    /* inheritance  */ public FieldPointBase<Derived,FieldType>
    {
        
        typedef typename FieldPoint<Derived,_dim,MoreFieldTypes...>::VectorDimD VectorDimD;
        
//        /**********************************************************************/
//        FieldPoint(const VectorDimD& Pin) :
//        /* base init */ FieldPoint<Derived,_dim,MoreFieldTypes...>(Pin)
//        {/*! @param[in] Pin the position of this FieldPoint
//          * Constructor initializes base class
//          */
//        }
        
        template <typename OtherFieldType>
        FieldPointBase<Derived,OtherFieldType>& field()
        {
            return *static_cast<FieldPointBase<Derived,OtherFieldType>* const>(this);
        }
        
        template <typename OtherFieldType>
        const FieldPointBase<Derived,OtherFieldType>& field() const
        {
            return *static_cast<const FieldPointBase<Derived,OtherFieldType>* const>(this);
        }
        
#ifdef _MODEL_MPI_
        
        typedef typename FieldPointBase<Derived,FieldType>::EigenMapType EigenMapType;
        
        void set_mpiID(const size_t& k)
        {
            FieldPoint<Derived,_dim,MoreFieldTypes...>::set_mpiID(k);
            FieldPointBase<Derived,FieldType>::set_mpiID(k);
        }
#endif
        
    };
    
    /**************************************************************************/
    /**************************************************************************/
    // template specialization: zero fields
    template<typename Derived,
    /*    */ short unsigned int _dim>
    struct FieldPoint<Derived,_dim>
    {
        
        typedef SpatialCellObserver<Derived,_dim> SpatialCellObserverType;
        typedef typename SpatialCellObserverType::VectorDimD VectorDimD;

#ifdef _MODEL_MPI_
    private:
        size_t _mpiID;
#endif

    public:
        
        //! A const reference to the Derived position
//        const VectorDimD& fieldP;
        
        
#ifdef _MODEL_MPI_
        
        /**********************************************************************/
        const size_t& mpiID() const
        {
            return _mpiID;
        }
        
        void set_mpiID(const size_t& k)
        {
            _mpiID=k;
        }
        
//        /**********************************************************************/
//        FieldPoint(const VectorDimD& p) :
//        /* init list */ fieldP(p),
//        /* init list */ _mpiID(0)
//        {/*
//          */
//        }
#else
//        const size_t mpiID;
        
//        /**********************************************************************/
//        FieldPoint(const VectorDimD& p) :
//        /* init list */ fieldP(p)
//        {/*
//          */
//        }
#endif
        
        
//        template <typename SourcePointType>
//        SpatialCellObserver<SourcePointType,_dim> asObserver() const
//        {
////            return SpatialCellObserver<SourcePointType,_dim>(fieldP);
//            return SpatialCellObserver<SourcePointType,_dim>(static_cast<const Derived*>(this)->P);
//        }

        template <typename SourcePointType>
        typename SpatialCellObserver<SourcePointType,_dim>::CellMapType neighborCells() const
        {
            //            return SpatialCellObserver<SourcePointType,_dim>(fieldP);
            return SpatialCellObserver<SourcePointType,_dim>::neighborCells(static_cast<const Derived*>(this)->P);
        }
        
        
    };
        
} // end namespace
#endif
