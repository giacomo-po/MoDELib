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
#include <model/DislocationDynamics/Materials/PeriodicElement.h>
#include <model/DislocationDynamics/Materials/FCCcrystal.h>
#include <model/DislocationDynamics/Materials/BCCcrystal.h>
#include <model/DislocationDynamics/Materials/SlipSystem.h>
//#include <model/Math/BestRationalApproximation.h>

namespace model
{
    
    template <int dim>
    class GrainBoundary;
    
    template <int dim>
    class Grain : public Lattice<dim>,
    /* base    */ public std::map<std::pair<size_t,size_t>,const GrainBoundary<dim>* const>
    {
        
        //typedef Simplex<dim,dim> SimplexType;
        typedef Lattice<dim> LatticeType;
        typedef MeshRegion<Simplex<dim,dim> > MeshRegionType;
        typedef MeshRegionObserver<MeshRegionType> MeshRegionObserverType;
        
        
        typedef Eigen::Matrix<long int,dim,1> VectorDimI;
        
        typedef Eigen::Matrix<  double,dim,1> VectorDimD;
        typedef Eigen::Matrix<double,dim,dim> MatrixDimD;
        typedef LatticeVector<dim> LatticeVectorType;
        typedef LatticeDirection<dim> LatticeDirectionType;
        
        typedef ReciprocalLatticeVector<dim> ReciprocalLatticeVectorType;
        typedef ReciprocalLatticeDirection<dim> ReciprocalLatticeDirectionType;
        
        typedef std::vector<LatticePlaneBase> PlaneNormalContainerType;
        typedef std::vector<SlipSystem> SlipSystemContainerType;
        typedef std::vector<unsigned int> PlaneNormalIDContainerType;
        
        
        static constexpr PeriodicElement<13,Isotropic> Al=PeriodicElement<13,Isotropic>();
        static constexpr PeriodicElement<28,Isotropic> Ni=PeriodicElement<28,Isotropic>();
        static constexpr PeriodicElement<29,Isotropic> Cu=PeriodicElement<29,Isotropic>();
        static constexpr PeriodicElement<74,Isotropic>  W=PeriodicElement<74,Isotropic>();
        
//        /**********************************************************************/
//        static Eigen::Matrix<long int,dim,1> rationalApproximation(VectorDimD nd)
//        {
//            const Eigen::Array<double,dim,1> nda(nd.array().abs()); // vector of close-to-integer numbers corresponding to lattice coordinates
//            size_t maxID=0;
//            const double maxVal(nda.maxCoeff(&maxID));
//            nd/=maxVal; // make each value of nd in [-1:1]
//            
//            Eigen::Array<long int,dim,1> nums=Eigen::Matrix<long int,dim,1>::Ones();
//            Eigen::Array<long int,dim,1> dens=Eigen::Matrix<long int,dim,1>::Ones();
//            long int denProd=1;
//            
//            for(int k=0;k<dim;++k)
//            {
//                BestRationalApproximation bra(nd(k),100);
//                
//                nums(k)=bra.num;
//                dens(k)=bra.den;
//                denProd*=bra.den;
//            }
//            
//            for(int k=0;k<dim;++k)
//            {
//                nums(k)*=(denProd/dens(k));
//            }
//            return nums.matrix();
//        }
        
        /**********************************************************************/
        void setLatticeBasis()
        {
//            model::cout<<"Grain"<<grainID<<", material="<<materialZ<<std::endl;
            
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
            
            Lattice<dim>::setLatticeBasis(C2G*A);
//            _covBasis=C2G*A;
//            _contraBasis=_covBasis.inverse().transpose();
//            
//            std::cout<<"Lattice basis (in columns) =\n"<<_covBasis<<std::endl;
//            std::cout<<"Lattice reciprocal basis (in columns) =\n"<<_contraBasis<<std::endl;
        }
        
        PlaneNormalContainerType planeNormalContainer;
        SlipSystemContainerType slipSystemContainer;
        
        //! The static column matrix of lattice vectors
//        MatrixDimD    _covBasis;
//        MatrixDimD _contraBasis;
        
        Eigen::Matrix<double,dim,dim> C2G;
        
        int materialZ;
        
    public:
        static constexpr double roundTol=FLT_EPSILON;
        
        const MeshRegionType& region;
        const int& grainID;
        
        /**********************************************************************/
        Grain(const MeshRegionType& region_in,
              const int& Z,
              const Eigen::Matrix<double,dim,dim>& C2G_in) :
//        /* init */ _covBasis(MatrixDimD::Identity()),
//        /* init */ _contraBasis(MatrixDimD::Identity()),
        /* init */ C2G(C2G_in),
        /* init */ materialZ(Z),
        /* init */ region(region_in),
        /* init */ grainID(region.regionID)
        {
            model::cout<<greenBoldColor<<"Creating Grain "<<grainID<<defaultColor<<std::endl;
            selectMaterial(materialZ);
            rotate(C2G);
        }
        
        const LatticeType& lattice() const
        {
            return *this;
        }
        
        /**********************************************************************/
        const std::map<std::pair<size_t,size_t>,const GrainBoundary<dim>* const> grainBoundaries() const
        {
            return *this;
        }
        
        /**********************************************************************/
        void selectMaterial(const int& Z)
        {
            model::cout<<"  grain "<<grainID<<", selecting material"<<defaultColor<<std::endl;
            
            materialZ=Z;
            setLatticeBasis();
            
            switch (materialZ)
            {
                case Al.Z:
                    planeNormalContainer=PeriodicElement<Al.Z,Isotropic>::CrystalStructure::reciprocalPlaneNormals(lattice());
                    slipSystemContainer=PeriodicElement<Al.Z,Isotropic>::CrystalStructure::slipSystems(lattice());
                    break;
                case Ni.Z:
                    planeNormalContainer=PeriodicElement<Ni.Z,Isotropic>::CrystalStructure::reciprocalPlaneNormals(lattice());
                    slipSystemContainer=PeriodicElement<Ni.Z,Isotropic>::CrystalStructure::slipSystems(lattice());
                    break;
                case Cu.Z:
                    planeNormalContainer=PeriodicElement<Cu.Z,Isotropic>::CrystalStructure::reciprocalPlaneNormals(lattice());
                    slipSystemContainer=PeriodicElement<Cu.Z,Isotropic>::CrystalStructure::slipSystems(lattice());
                    break;
                case W.Z:
                    planeNormalContainer=PeriodicElement<W.Z,Isotropic>::CrystalStructure::reciprocalPlaneNormals(lattice());
                    slipSystemContainer=PeriodicElement<W.Z,Isotropic>::CrystalStructure::slipSystems(lattice());
                    break;
                    //                case Fe:
                    //                    CrystalOrientation<dim>::template rotate<PeriodicElement<Fe,Isotropic>::CrystalStructure>(C2G);
                    //                    break;
                default:
                    assert(0 && "Material not implemented.");
                    break;
            }
            
//            model::cout<<magentaColor<<"Current Crystal Plane Normals are:"<<std::endl;
//            for (unsigned int k=0; k<planeNormalContainer.size();++k)
//            {
//                model::cout<<"    "<<planeNormalContainer[k].cartesian().normalized().transpose()<<std::endl;
//            }
//            model::cout<<defaultColor<<std::endl;
            
        }
        
        /**********************************************************************/
        void rotate(const Eigen::Matrix<double,dim,dim>& C2G_in)
        {/*! Z is atomic number
          */
            model::cout<<"  grain "<<grainID<<", rotating"<<defaultColor<<std::endl;
            
            assert((C2G_in*C2G_in.transpose()-Eigen::Matrix<double,dim,dim>::Identity()).norm()<2.0*DBL_EPSILON*dim*dim && "CRYSTAL TO GLOBAL ROTATION MATRIX IS NOT ORTHOGONAL.");
            // make sure that C2G is proper
            assert(std::fabs(C2G_in.determinant()-1.0) < FLT_EPSILON && "C2G IS NOT PROPER.");
            // store C2G
            C2G=C2G_in;
            setLatticeBasis();
            
//            model::cout<<magentaColor<<"Current Crystal Plane Normals are:"<<std::endl;
//            for (unsigned int k=0; k<planeNormalContainer.size();++k)
//            {
//                model::cout<<"    "<<planeNormalContainer[k].cartesian().normalized().transpose()<<std::endl;
//            }
//            model::cout<<defaultColor<<std::endl;
        }
        
//        /**********************************************************************/
//        LatticeVectorType snapToLattice(const VectorDimD& d) const
//        {
//            VectorDimD nd(_contraBasis.transpose()*d);
//            return LatticeVectorType(RoundEigen<double,dim>::round(nd).template cast<long int>(),_covBasis,_contraBasis);
//        }
        
//        /**********************************************************************/
//        LatticeDirectionType latticeDirection(const VectorDimD& d) const
//        {
//            
//            const VectorDimD nd(_contraBasis.transpose()*d);
//            const LatticeVectorType temp(rationalApproximation(nd),_covBasis,_contraBasis);
//            
//            if(temp.cartesian().normalized().cross(d.normalized()).norm()>FLT_EPSILON)
//            {
//                std::cout<<"input direction="<<d.normalized().transpose()<<std::endl;
//                std::cout<<"lattice direction="<<temp.cartesian().normalized().transpose()<<std::endl;
//                assert(0 && "LATTICE DIRECTION NOT FOUND");
//            }
//            
//            return LatticeDirectionType(temp);
//        }
        
//        /**********************************************************************/
//        ReciprocalLatticeDirectionType reciprocalLatticeDirection(const VectorDimD& d) const
//        {
//            
//            const VectorDimD nd(_covBasis.transpose()*d);
//            const ReciprocalLatticeVectorType temp(rationalApproximation(nd),_covBasis,_contraBasis);
//            
//            if(temp.cartesian().normalized().cross(d.normalized()).norm()>FLT_EPSILON)
//            {
//                std::cout<<"input direction="<<d.normalized().transpose()<<std::endl;
//                std::cout<<"reciprocal lattice direction="<<temp.cartesian().normalized().transpose()<<std::endl;
//                assert(0 && "RECIPROCAL LATTICE DIRECTION NOT FOUND");
//            }
//            
//            return ReciprocalLatticeDirectionType(temp);
//        }
        
        /**********************************************************************/
        std::pair<LatticePlane,LatticePlane> find_confiningPlanes(const LatticeVectorType& sourceL,
                                                                  const LatticeVectorType& sinkL,
                                                                  const LatticeVectorType& Burgers) const
        {
            
            const LatticeVectorType chord(sinkL-sourceL);
            
            bool isOnGB=false;
            const GrainBoundary<dim>* p_GB=NULL;
            for(const auto& gb : grainBoundaries())
            {
                if(   gb.second->latticePlane(grainID).contains(sourceL)
                   && gb.second->latticePlane(grainID).contains(sinkL)
                   )
                {
                    isOnGB=true;
                    p_GB=gb.second;
                    break;
                }
            }
            
            
            // Find the crystallographic planes which may contain chord and Burgers
            PlaneNormalContainerType allowedGLidePlanes;
            for (const auto& planeBase : planeNormals()) // Loop over crystal planes
            {
                if(planeBase.dot(chord)==0
                   && planeBase.dot(Burgers)==0)
                {
                    allowedGLidePlanes.push_back(planeBase);
                }
            }
            if(isOnGB)
            {
                if(   p_GB->latticePlane(grainID).n.dot(chord)==0
                   && p_GB->latticePlane(grainID).n.dot(Burgers)==0
                   )
                {
                    allowedGLidePlanes.push_back(p_GB->latticePlane(grainID).n);
                }
            }
            
            // Find the crystallographic planes which may contain chord and Burgers
            PlaneNormalContainerType allowedSessilePlanes;
            for (const auto& planeBase : planeNormals()) // Loop over crystal planes
            {
                if(planeBase.dot(chord)==0)
                {
                    allowedSessilePlanes.push_back(planeBase);
                }
            }
            if(isOnGB)
            {
                allowedSessilePlanes.push_back(p_GB->latticePlane(grainID).n);
            }
            
            // Return the two planes based on the following conditions:
            if(allowedGLidePlanes.size()) // At least glide plane found
            {
                if(isOnGB) // possibly glissile segmento on GB, or sessile at the intersection of GP and GB
                {
                    return std::make_pair(LatticePlane(sourceL,allowedGLidePlanes[0]),LatticePlane(p_GB->latticePlane(grainID)));
                }
                else // bulk segment, purely glissile
                {
                    return std::make_pair(LatticePlane(sourceL,allowedGLidePlanes[0]),LatticePlane(sourceL,allowedGLidePlanes[0]));
                }
            }
            else // No glide planes found. Check for possibly sessile segment
            {
                if(allowedSessilePlanes.size()>=2)
                {
                    return std::make_pair(LatticePlane(sourceL,allowedSessilePlanes[0]),LatticePlane(sourceL,allowedSessilePlanes[1]));
                }
                else
                {
                    assert(0 && "SESSILE SEGMENTS MUST FORM ON THE INTERSECTION OF TWO PLANES.");
                }
            }
            
        }
        
        
        
        /**********************************************************************/
        std::deque<const LatticePlaneBase*> conjugatePlaneNormal(const LatticeVectorType& B,
                                                                 const ReciprocalLatticeDirectionType& N) const
        {
            std::deque<const LatticePlaneBase*> temp;
            if(B.dot(N)==0) // not sessile
            {
                for (const auto& planeNormal : planeNormalContainer)
                {
                    if(	 B.dot(planeNormal)==0 && N.cross(planeNormal).squaredNorm()>0)
                    {
                        temp.push_back(&planeNormal);
                    }
                }
            }
            return temp;
        }
        
//        /**********************************************************************/
//        const MatrixDimD& covBasis() const
//        {
//            return _covBasis;
//        }
//        
//        /**********************************************************************/
//        const MatrixDimD& contraBasis() const
//        {
//            return _contraBasis;
//        }
        
//        /**********************************************************************/
//        LatticeVectorType latticeVector(const VectorDimD& p) const
//        {
//            return LatticeVectorType(p,_covBasis,_contraBasis);
//        }
//        
//        /**********************************************************************/
//        ReciprocalLatticeVectorType reciprocalLatticeVector(const VectorDimD& p) const
//        {
//            return ReciprocalLatticeVectorType(p,_covBasis,_contraBasis);
//        }
        
        /**********************************************************************/
        const MatrixDimD& get_C2G() const
        {
            return C2G;
        }
        
        /**********************************************************************/
        const PlaneNormalContainerType& planeNormals() const
        {
            return planeNormalContainer;
        }
        
        /**********************************************************************/
        const SlipSystemContainerType& slipSystems() const
        {
            return slipSystemContainer;
        }
        
    };
    
    
} // end namespace
#endif

//        /**********************************************************************/
//        PlaneNormalIDContainerType find_slipSystem(const LatticeVectorType& chord,
//                                                   const LatticeVectorType& Burgers,
//                                                   const bool& isGBsegment) const
//        {/*!
//          */
//            assert(  chord.squaredNorm()>0 && "CHORD HAS ZERO NORM");
//            assert(Burgers.squaredNorm()>0 && "BURGERS HAS ZERO NORM");
//
//
//            PlaneNormalIDContainerType allowedSlipSystems;
//
//            // Try to find a plane which has normal orthogonal to both the chord and the Burgers
//            for (unsigned int k=0;k<planeNormalContainer.size();++k)
//            {
//                if(planeNormalContainer[k].dot(chord)==0
//                   && planeNormalContainer[k].dot(Burgers)==0)
//                {
//                    allowedSlipSystems.push_back(k);
//                }
//            }
//
//            // If no planes are found, check only chord to detect possibly sessile segment
//            if(allowedSlipSystems.size()+isGBsegment==0)
//            {
//                for (unsigned int k=0;k<planeNormalContainer.size();++k)
//                {
//                    if(	planeNormalContainer[k].dot(chord)==0)
//                    {
//                        allowedSlipSystems.push_back(k);
//                    }
//                }
//                if (allowedSlipSystems.size()<2)
//                {
//                    std::cout<<"chord="<<chord.cartesian().transpose()<<std::endl;
//                    std::cout<<"Burgers="<<Burgers.cartesian().transpose()<<std::endl;
//                    for (const auto& planeNormal : planeNormalContainer)
//                    {
//
//                        std::cout<<"n="<<planeNormal.cartesian().normalized().transpose()<<" c*n="<< planeNormal.dot(chord)<<" b*n="<< planeNormal.dot(Burgers) <<std::endl;
//                    }
//                    assert(allowedSlipSystems.size()>=2 && "SESSILE SEGMENTS MUST FORM ON THE INTERSECTION OF TWO CRYSTALLOGRAPHIC PLANES.");
//                }
//            }
//
//            return allowedSlipSystems;
//        }
//
//        /**********************************************************************/
//        const LatticePlaneBase& find_glidePlane(const LatticeVectorType& sourceL,
//                                                const LatticeVectorType& sinkL,
//                                                const LatticeVectorType& Burgers) const
//        {/*!@param[in] sourceL the source node position of a DislocationSegment
//          * @param[in] sourceL the sink node position of a DislocationSegment
//          * @param[in] Burgers the Burgers vector of a DislocationSegment
//          *\returns A const reference to the first vector in planeNormalContainer
//          * which is orthogonal to both chord and Burgers.
//          */
//
//            bool isGBsegment=false;
//            const GrainBoundary<dim>* p_GB=NULL;
//            for(const auto& gb : grainBoundaries())
//            {
//                //                std::cout<<"CHECKING GRAINBOUNDARIES"<<std::endl;
//                //                std::cout<<gb.second->latticePlane(grainID).contains(sourceL)<<std::endl;
//                //                std::cout<<gb.second->latticePlane(grainID).contains(sinkL)<<std::endl;
//                //                std::cout<<gb.second->latticePlane(grainID).n.dot(Burgers)<<std::endl;
//
//                if(   gb.second->latticePlane(grainID).contains(sourceL)
//                   && gb.second->latticePlane(grainID).contains(sinkL)
//                   && gb.second->latticePlane(grainID).n.dot(Burgers)==0
//                   )
//                {
//                    isGBsegment=true;
//                    p_GB=gb.second;
//                    break;
//                }
//            }
//
//            const PlaneNormalIDContainerType allowedSlipSystems=find_slipSystem(sinkL-sourceL,Burgers,isGBsegment);
//
//            return allowedSlipSystems.size()>0? planeNormalContainer[allowedSlipSystems[0]] : p_GB->latticePlane(grainID).n ; // RETURNING THE FIRST PLANE FOUND IS SOMEWHAT ARBITRARY
//        }
//
//        /**********************************************************************/
//        const LatticePlaneBase& find_sessilePlane(const LatticeVectorType& sourceL,
//                                                  const LatticeVectorType& sinkL,
//                                                  const LatticeVectorType& Burgers) const
//        {/*!@param[in] sourceL the source node position of a DislocationSegment
//          * @param[in] sinkL the sink node position of a DislocationSegment
//          * @param[in] Burgers the Burgers vector of a DislocationSegment
//          *\returns A const reference to the a second LatticePlaneBase
//          * which contains the chord.
//
//          * If both sourceL and sinkL are contained in one of the Grainboundaries
//          * of this Grain, then the corresponding LatticePlaneBase is returned.
//          * Otherwise, a search within the crystallographic planes is performed.
//          * If the search returns one plane, then the sessilePlane os the same as
//          * the glidePlane, and the segment is glissile on that plane. If the search
//          * returns two planes, then the glidePlane is returned for a screw segment,
//          * while the other plane is returned for a non-screw segment.
//          */
//
//            // if the segment is on a GB, always return the GB as the sessile plane
//            for(const auto& gb : grainBoundaries())
//            {
//                if(   gb.second->latticePlane(grainID).contains(sourceL)
//                   && gb.second->latticePlane(grainID).contains(sinkL)
//                   //                   && gb.second->latticePlane(grainID).n.dot(Burgers)==0
//                   )
//                {
//                    return gb.second->latticePlane(grainID).n;
//                }
//            }
//
//            const PlaneNormalIDContainerType allowedSlipSystems=find_slipSystem(sinkL-sourceL,Burgers,false);
//
//            int planeID(0);
//            if (allowedSlipSystems.size()>=2)
//            {
//                if ((sinkL-sourceL).cross(Burgers).squaredNorm()!=0) // a sessile segment
//                {
//                    planeID=1;
//                }
//                else // a screw segment
//                {
//                    planeID=0; // allow glide on primary plane
//                }
//            }
//
//            return planeNormalContainer[allowedSlipSystems[planeID]];
//        }


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

//        /**********************************************************************/
//        LatticeDirectionType latticeDirection(const VectorDimD& d) const
//        {
//            bool found=false;
//            VectorDimD rdk(VectorDimD::Zero());
//
//            const VectorDimD nd(_contraBasis.transpose()*d);
//
//
//            for(int k=0;k<dim;++k)
//            {
//                if(fabs(nd(k))>roundTol)
//                {
//                    const VectorDimD ndk(nd/fabs(nd(k)));
//                    rdk=RoundEigen<double,dim>::round(ndk);
//                    if((ndk-rdk).norm()<roundTol)
//                    {
//                        found=true;
//                        break;
//                    }
//                }
//            }
//            assert(found && "Input vector is not on a lattice direction");
//            return LatticeDirectionType(rdk.template cast<long int>(),_covBasis,_contraBasis);
//        }
