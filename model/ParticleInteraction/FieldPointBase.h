/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2012 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _model_FieldPointBase_h
#define _model_FieldPointBase_h

namespace model
{
    
    /*!\brief A class template that contains the underlying FieldType data of
     * the Derived particle.
     *
     * If MPI is not used, the data is stored locally in the base class
     * Eigen::Matrix<scalar,rows,cols>.
     *
     * In case that MPI is used, the underlying data is not stored with the 
     * object but is delocalized to a static vector that can efficiently be
     * syncronized using MPI_Allgatherv. However, the base data is conveniently 
     * used as regular Eigen object via the Eigen::Map class.
     */    
    template<typename Derived, typename FieldType>
    struct FieldPointBase :
#ifdef _MODEL_MPI_
    /* inheritance  */ public Eigen::Map<Eigen::Matrix<typename FieldType::Scalar,FieldType::rows,FieldType::cols> >
#else
    /* inheritance  */ public Eigen::Matrix<typename FieldType::Scalar,FieldType::rows,FieldType::cols>
#endif
    {        
        enum{rows=FieldType::rows};
        enum{cols=FieldType::cols};
        enum{DataPerParticle=rows*cols};
        typedef typename FieldType::Scalar Scalar;
        typedef FieldPointBase<Derived,FieldType> FieldPointBaseType;
        
        const bool enabled;

        
#ifdef _MODEL_MPI_
        static std::vector<Scalar> resultVector;
        typedef Eigen::Map<Eigen::Matrix<typename FieldType::Scalar,FieldType::rows,FieldType::cols> > EigenMapType;
        typedef EigenMapType BaseEigenType;
        
        
        /**********************************************************************/
        static void resize(const unsigned int&  k, const Scalar& val = Scalar())
        {
            resultVector.resize(k*DataPerParticle,val);
        }
        
        /**********************************************************************/
        void set_mpiID(const size_t& k)
        {
            new (static_cast<EigenMapType*>(this)) EigenMapType(&resultVector[k*DataPerParticle]);
        }
        
        /**********************************************************************/
        FieldPointBase(const bool& enb) :
        /* init list */ Eigen::Map<Eigen::Matrix<Scalar,rows,cols> >(NULL),
        /* init list */ enabled(enb)
        {/*
          */
        }
        
#else
        typedef Eigen::Matrix<typename FieldType::Scalar,FieldType::rows,FieldType::cols> BaseEigenType;
        
        /**********************************************************************/
        FieldPointBase(const bool& enb) :
        /* init list */ Eigen::Matrix<Scalar,rows,cols>(Eigen::Matrix<Scalar,rows,cols>::Zero()),
        /* init list */ enabled(enb)
        {/*
          */
        }
#endif
        
        template<typename OtherDerived>
        FieldPointBaseType & operator= (const Eigen::MatrixBase<OtherDerived>& other)
        {
            BaseEigenType::operator=(other);
            return *this;
        }
        
    };
    
    // declare statica data
#ifdef _MODEL_MPI_
    template<typename Derived,
    /*    */ typename FieldType>
    std::vector<typename FieldType::Scalar> FieldPointBase<Derived,FieldType>::resultVector;
#endif
    
} // end namespace
#endif
