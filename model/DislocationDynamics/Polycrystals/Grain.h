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
#include <model/DislocationDynamics/Materials/PeriodicElement.h>
#include <model/DislocationDynamics/Materials/FCCcrystal.h>
#include <model/DislocationDynamics/Materials/BCCcrystal.h>

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
        
        typedef std::vector<LatticePlaneBase> PlaneNormalContainerType;

        static constexpr PeriodicElement<13,Isotropic> Al=PeriodicElement<13,Isotropic>();
        static constexpr PeriodicElement<28,Isotropic> Ni=PeriodicElement<28,Isotropic>();
        static constexpr PeriodicElement<29,Isotropic> Cu=PeriodicElement<29,Isotropic>();
        static constexpr PeriodicElement<74,Isotropic>  W=PeriodicElement<74,Isotropic>();

        /**********************************************************************/
        void setLatticeBasis()
        {
            model::cout<<"Grain"<<grainID<<", material="<<materialZ<<std::endl;

            Eigen::Matrix<double,dim,dim,1> A(Eigen::Matrix<double,dim,dim,1>::Identity());
            
            switch (materialZ)
            {
                case Al.Z:
                    A=PeriodicElement<Al.Z,Isotropic>::CrystalStructure::template getLatticeBasis<dim>();
                    break;
                case Ni.Z:
                    A=PeriodicElement<Ni.Z,Isotropic>::CrystalStructure::template getLatticeBasis<dim>();
                    break;
                case Cu.Z:
                    A=PeriodicElement<Cu.Z,Isotropic>::CrystalStructure::template getLatticeBasis<dim>();
                    break;
                case W.Z:
                    A=PeriodicElement<W.Z,Isotropic>::CrystalStructure::template getLatticeBasis<dim>();
                    break;
                    //                case Fe:
                    //                    CrystalOrientation<dim>::template rotate<PeriodicElement<Fe,Isotropic>::CrystalStructure>(C2G);
                    //                    break;
                default:
                    assert(0 && "Material not implemented.");
                    break;
            }
            
            _covBasis=C2G*A;
            _contraBasis=_covBasis.inverse().transpose();
            
            std::cout<<"Lattice basis (in columns) =\n"<<_covBasis<<std::endl;
            std::cout<<"Lattice reciprocal basis (in columns) =\n"<<_contraBasis<<std::endl;
        }
        
//        /**********************************************************************/
//        void setLatticeBasis(const Eigen::Matrix<double,dim,dim,1>& A)
//        {
//            _covBasis=A;
//            _contraBasis=A.inverse().transpose();
//            
//            std::cout<<"Lattice basis (in columns) =\n"<<_covBasis<<std::endl;
//            std::cout<<"Lattice reciprocal basis (in columns) =\n"<<_contraBasis<<std::endl;
//        }
        
//        /**********************************************************************/
//        template<size_t Z>
//        void select()
//        {/*! Z is atomic number
//          */
//            setLatticeBasis(C2G*PeriodicElement<Z,Isotropic>::CrystalStructure::template getLatticeBasis<dim>());
//            planeNormalContainer=CrystalStructure::template reciprocalPlaneNormals<dim>();
//            //slipSystemContainer=CrystalStructure::slipSystems();
//        }
        
        PlaneNormalContainerType planeNormalContainer;
        
        //! The static column matrix of lattice vectors
        MatrixDimD    _covBasis;
        MatrixDimD _contraBasis;
        
        Eigen::Matrix<double,dim,dim> C2G;
        
        int materialZ;
        
    public:
        static constexpr double roundTol=FLT_EPSILON;

        const MeshRegionType& region;
        const int& grainID;
        
        /**********************************************************************/
        Grain(const MeshRegionType& region_in) :
        /* init */ _covBasis(MatrixDimD::Identity()),
        /* init */ _contraBasis(MatrixDimD::Identity()),
        /* init */ C2G(Eigen::Matrix<double,dim,dim>::Identity()),
        /* init */ materialZ(29),
        /* init */ region(region_in),
        /* init */ grainID(region.regionID)
        {
            model::cout<<"Creating Grain "<<grainID<<std::endl;
            selectMaterial(materialZ);
        }
        
        /**********************************************************************/
        void selectMaterial(const int& Z)
        {
            model::cout<<greenColor<<"Grain "<<grainID<<", selecting material"<<defaultColor<<std::endl;

            materialZ=Z;
            setLatticeBasis();
            
            switch (materialZ)
            {
                case Al.Z:
                    planeNormalContainer=PeriodicElement<Al.Z,Isotropic>::CrystalStructure::reciprocalPlaneNormals(_covBasis,_contraBasis);
                    break;
                case Ni.Z:
                    planeNormalContainer=PeriodicElement<Ni.Z,Isotropic>::CrystalStructure::reciprocalPlaneNormals(_covBasis,_contraBasis);
//
//                    CrystalOrientation<dim>::template rotate<typename PeriodicElement<Ni.Z,Isotropic>::CrystalStructure>(C2G);
                    break;
                case Cu.Z:
                    planeNormalContainer=PeriodicElement<Cu.Z,Isotropic>::CrystalStructure::reciprocalPlaneNormals(_covBasis,_contraBasis);

//                    CrystalOrientation<dim>::template rotate<typename PeriodicElement<Cu.Z,Isotropic>::CrystalStructure>(C2G);
                    break;
                case W.Z:
                    planeNormalContainer=PeriodicElement<W.Z,Isotropic>::CrystalStructure::reciprocalPlaneNormals(_covBasis,_contraBasis);

//                    CrystalOrientation<dim>::template rotate<typename PeriodicElement<W.Z,Isotropic>::CrystalStructure>(C2G);
                    break;
                    //                case Fe:
                    //                    CrystalOrientation<dim>::template rotate<PeriodicElement<Fe,Isotropic>::CrystalStructure>(C2G);
                    //                    break;
                default:
                    assert(0 && "Material not implemented.");
                    break;
            }
            
            model::cout<<magentaColor<<"Current Crystal Plane Normals are:"<<std::endl;
            for (unsigned int k=0; k<planeNormalContainer.size();++k)
            {
                std::cout<<"    "<<planeNormalContainer[k].cartesian().normalized().transpose()<<std::endl;
            }
            model::cout<<defaultColor<<std::endl;

            
        }
        

        
        /**********************************************************************/
        void rotate(const Eigen::Matrix<double,dim,dim>& C2G_in)
        {/*! Z is atomic number
          */
            model::cout<<greenColor<<"Grain "<<grainID<<", rotating"<<defaultColor<<std::endl;
            
            assert((C2G_in*C2G_in.transpose()-Eigen::Matrix<double,dim,dim>::Identity()).norm()<2.0*DBL_EPSILON*dim*dim && "CRYSTAL TO GLOBAL ROTATION MATRIX IS NOT ORTHOGONAL.");
            // make sure that C2G is proper
            assert(std::fabs(C2G_in.determinant()-1.0) < FLT_EPSILON && "C2G IS NOT PROPER.");
            // store C2G
            C2G=C2G_in;
            setLatticeBasis();
            
            model::cout<<magentaColor<<"Current Crystal Plane Normals are:"<<std::endl;
            for (unsigned int k=0; k<planeNormalContainer.size();++k)
            {
                std::cout<<"    "<<planeNormalContainer[k].cartesian().normalized().transpose()<<std::endl;
            }
            model::cout<<defaultColor<<std::endl;
        }
        
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
        
        /**********************************************************************/
        const MatrixDimD& contraBasis() const
        {
            return _contraBasis;
        }
        
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

