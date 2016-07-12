/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_Grain_H_
#define model_Grain_H_

#include <cfloat> // FLT_EPSILON
#include <Eigen/Core>
#include <model/Math/RoundEigen.h>
#include <model/Mesh/SimplicialMesh.h>
#include <model/Mesh/MeshRegionObserver.h>
#include <model/LatticeMath/LatticeMath.h>
//#include <model/LatticeMath/LatticeVector.h>
//#include <model/LatticeMath/ReciprocalLatticeVector.h>


namespace model
{
    
    
    
    template <int dim>
    class Grain
    {
        
        //typedef Simplex<dim,dim> SimplexType;
        typedef MeshRegion<Simplex<dim,dim> > MeshRegionType;
        typedef MeshRegionObserver<MeshRegionType> MeshRegionObserverType;
        
        
        typedef Eigen::Matrix<long int,dim,1> VectorDimI;

        typedef Eigen::Matrix<  double,dim,1> VectorDimD;
        typedef Eigen::Matrix<double,dim,dim> MatrixDimD;
        typedef LatticeVector<dim> LatticeVectorType;
        typedef ReciprocalLatticeVector<dim> ReciprocalLatticeVectorType;
        
        
        //! The static column matrix of lattice vectors
        MatrixDimD    _covBasis;
//        MatrixDimD    AT;
//        MatrixDimD invA;
        MatrixDimD _contraBasis;
        
    public:
        static constexpr double roundTol=FLT_EPSILON;

        const MeshRegionType& region;
        
        /**********************************************************************/
        Grain(const MeshRegionType& region_in) :
        /* init */ _covBasis(MatrixDimD::Identity()),
//        /* init */ AT(A.transpose()),
//        /* init */ invA(A.inverse()),
        /* init */ _contraBasis(MatrixDimD::Identity()),
        /* init */ region(region_in)
        {
            model::cout<<"Creating Grain "<<region.regionID<<std::endl;
            
//            for(const auto& simplex_p : region)
//            {
//            if(simplex_p)
//            {
//            
//            }
//            }
            
        }
        
        /**********************************************************************/
        void setLatticeBasis(const Eigen::Matrix<double,dim,dim,1>& A)
        {
            _covBasis=A;
//            AT=A.transpose();
//            invA=A.inverse();
            _contraBasis=A.inverse().transpose();
            //            cofA=invA*A.determinant();
            
            std::cout<<"Lattice basis (in columns) =\n"<<_covBasis<<std::endl;
            std::cout<<"Lattice reciprocal basis (in columns) =\n"<<_contraBasis<<std::endl;
            
        }
        
//        /**********************************************************************/
//        VectorDimI d2contra(const VectorDimD& d) const
//        {
//            const VectorDimD nd(invA*d);
//            const VectorDimD rd(RoundEigen<double,dim>::round(nd));
//            if((nd-rd).norm()>roundTol)
//            {
//                std::cout<<"d2contra, nd="<<nd.transpose()<<std::endl;
//                std::cout<<"d2contra, rd="<<rd.transpose()<<std::endl;
//                assert(0 && "Input vector is not a lattice vector");
//            }
//            return rd.template cast<long int>();
//        }
        
        /**********************************************************************/
        VectorDimI latticeDirection(const VectorDimD& d) const
        {
            bool found=false;
            VectorDimD rdk(VectorDimD::Zero());
            
            const VectorDimD nd(_contraBasis.transpose()*d);
            
            
            for(int k=0;k<dim;++k)
            {
                const VectorDimD ndk(nd/nd(k));
                rdk=RoundEigen<double,dim>::round(ndk);
                if((ndk-rdk).norm()<roundTol)
                {
                    found=true;
                    break;
                }
                
            }
            assert(found && "Input vector is not on a lattice direction");
            return rdk.template cast<long int>();
        }
        
        /**********************************************************************/
        VectorDimI snapToLattice(const VectorDimD& d) const
        {
            VectorDimD nd(_contraBasis.transpose()*d);
            return RoundEigen<double,dim>::round(nd).template cast<long int>();
        }
        
//        /**********************************************************************/
//        VectorDimI d2cov(const VectorDimD& d) const
//        {
//            const VectorDimD nd(AT*d);
//            const VectorDimD rd(RoundEigen<double,dim>::round(nd));
//            if((nd-rd).norm()>roundTol)
//            {
//                std::cout<<"d2cov, nd="<<nd.transpose()<<std::endl;
//                std::cout<<"d2cov, rd="<<rd.transpose()<<std::endl;
//                assert(0 && "Input vector is not a reciprocal lattice vector");
//            }
//            //            assert((nd-rd).norm()<roundTol && "Input vector is not a lattice vector");
//            return rd.template cast<long int>();
//        }
        
        /**********************************************************************/
        VectorDimI reciprocalLatticeDirection(const VectorDimD& d) const
        {
            bool found=false;
            VectorDimD rdk(VectorDimD::Zero());
            
            const VectorDimD nd(_covBasis.tranpose()*d);
            
            
            for(int k=0;k<dim;++k)
            {
                const VectorDimD ndk(nd/nd(k));
                rdk=RoundEigen<double,dim>::round(ndk);
                if((ndk-rdk).norm()<roundTol)
                {
                    found=true;
                    break;
                }
                
            }
            assert(found && "Input vector is not on a lattice direction");
            return rdk.template cast<long int>();
        }
        
        /**********************************************************************/
        const MatrixDimD& covBasis() const
        {
            return _covBasis;
        }
        
//        /**********************************************************************/
//        const MatrixDimD& invCovBasis() const
//        {
//            return invA;
//        }
        
        /**********************************************************************/
        const MatrixDimD& contraBasis() const
        {
            return _contraBasis;
        }
        
//        /**********************************************************************/
//        const MatrixDimD& invContraBasis() const
//        {
//            return AT;
//        }
        
        /**********************************************************************/
        LatticeVectorType latticeVectorFromPosition(const VectorDimD& p) const
        {
            std::cout<<"Grain "<<region.regionID<<" creating LatticeVector"<<std::endl;
            return LatticeVectorType(p,_covBasis,_contraBasis);
        }
        
        /**********************************************************************/
        ReciprocalLatticeVectorType reciprocalLatticeVectorFromPosition(const VectorDimD& p) const
        {
            std::cout<<"Grain "<<region.regionID<<" creating ReciprocalLatticeVector"<<std::endl;
            return ReciprocalLatticeVectorType(p,_covBasis,_contraBasis);
        }
        
    };
    
    
} // end namespace
#endif

