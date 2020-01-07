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
//#include <RoundEigen.h>
#include <SimplicialMesh.h>
#include <MeshRegionObserver.h>
#include <MeshRegion.h>
#include <LatticeMath.h>
//#include <PeriodicElement.h>
#include <DislocatedMaterial.h>
#include <SlipSystem.h>
#include <SingleCrystal.h>

//#include <BestRationalApproximation.h>

namespace model
{
    
    template <int dim>
    class GrainBoundary;
    
    template <int dim>
    class Grain : public SingleCrystal<dim>,
    /* base    */ public std::map<std::pair<size_t,size_t>,const GrainBoundary<dim>* const>
    {
        
        typedef Lattice<dim> LatticeType;
        typedef SingleCrystal<dim> SingleCrystalType;
        typedef MeshRegion<Simplex<dim,dim> > MeshRegionType;
        typedef MeshRegionObserver<MeshRegionType> MeshRegionObserverType;
        
        
        typedef Eigen::Matrix<long int,dim,1> VectorDimI;
        
        typedef Eigen::Matrix<  double,dim,1> VectorDimD;
        typedef Eigen::Matrix<double,dim,dim> MatrixDimD;
        typedef LatticeVector<dim> LatticeVectorType;
        typedef LatticeDirection<dim> LatticeDirectionType;
        
        typedef ReciprocalLatticeDirection<dim> ReciprocalLatticeDirectionType;
        typedef std::map<std::pair<size_t,size_t>,const GrainBoundary<dim>* const> GrainBoundaryContainerType;
        

        
        
    public:
        

        
        static constexpr double roundTol=FLT_EPSILON;
        
        const MeshRegionType& region;
        const int& grainID;
        
        /**********************************************************************/
        Grain(const MeshRegionType& region_in,
              const DislocatedMaterial<dim,Isotropic>& material,
              const std::string& polyFile
              ) :
        /* init */ SingleCrystalType(material,TextFileParser(polyFile).readMatrix<double>("C2G"+std::to_string(region_in.regionID),dim,dim,true))
        /* init */,region(region_in)
        /* init */,grainID(region.regionID)
        {
            model::cout<<"  lattice basis="<<this->latticeBasis<<std::endl;
            model::cout<<"  # plane normals="<<this->planeNormals().size()<<std::endl;
            model::cout<<"  # slip systems="<<this->slipSystems().size()<<std::endl;            
        }
        
        /**********************************************************************/
        const GrainBoundaryContainerType& grainBoundaries() const
        {
            return *this;
        }
        
        /**********************************************************************/
        GrainBoundaryContainerType& grainBoundaries()
        {
            return *this;
        }
        
        /**********************************************************************/
        std::deque<const LatticePlaneBase*> conjugatePlaneNormal(const LatticeVectorType& B,
                                                                 const ReciprocalLatticeDirectionType& N) const
        {
            std::deque<const LatticePlaneBase*> temp;
            if(B.dot(N)==0) // not sessile
            {
                for (const auto& planeNormal : this->planeNormals())
                {
                    if(	 B.dot(planeNormal)==0 && N.cross(planeNormal).squaredNorm()>0)
                    {
                        temp.push_back(&planeNormal);
                    }
                }
            }
            return temp;
        }
        
    };
    
}
#endif


//        /**********************************************************************/
//        const int& material() const
//        {
//            return materialZ;
//        }

//        /**********************************************************************/
//        const PlaneNormalContainerType& planeNormals() const
//        {
//            return planeNormalContainer;
//        }
//
//        /**********************************************************************/
//        const SlipSystemContainerType& slipSystems() const
//        {
//            return slipSystemContainer;
//        }

//        /**********************************************************************/
//        const LatticeType& lattice() const
//        {
//            return *this;
//        }

//        /**********************************************************************/
//        void selectMaterial(const int& Z)
//        {
//            model::cout<<"  grain "<<grainID<<", selecting material"<<defaultColor<<std::endl;
//
//            materialZ=Z;
//            setLatticeBasis();
//
//            switch (materialZ)
//            {
//                case Al.Z:
//                    planeNormalContainer=PeriodicElement<Al.Z,Isotropic>::CrystalStructure::reciprocalPlaneNormals(lattice());
//                    slipSystemContainer=PeriodicElement<Al.Z,Isotropic>::CrystalStructure::slipSystems(lattice());
//                    break;
//                case Ni.Z:
//                    planeNormalContainer=PeriodicElement<Ni.Z,Isotropic>::CrystalStructure::reciprocalPlaneNormals(lattice());
//                    slipSystemContainer=PeriodicElement<Ni.Z,Isotropic>::CrystalStructure::slipSystems(lattice());
//                    break;
//                case Cu.Z:
//                    planeNormalContainer=PeriodicElement<Cu.Z,Isotropic>::CrystalStructure::reciprocalPlaneNormals(lattice());
//                    slipSystemContainer=PeriodicElement<Cu.Z,Isotropic>::CrystalStructure::slipSystems(lattice());
//                    break;
//                case W.Z:
//                    planeNormalContainer=PeriodicElement<W.Z,Isotropic>::CrystalStructure::reciprocalPlaneNormals(lattice());
//                    slipSystemContainer=PeriodicElement<W.Z,Isotropic>::CrystalStructure::slipSystems(lattice());
//                    break;
//                    //                case Fe:
//                    //                    CrystalOrientation<dim>::template rotate<PeriodicElement<Fe,Isotropic>::CrystalStructure>(C2G);
//                    //                    break;
//                default:
//                    assert(0 && "Material not implemented.");
//                    break;
//            }
//
//        }

//        /**********************************************************************/
//        void rotate(const Eigen::Matrix<double,dim,dim>& C2G_in)
//        {/*! Z is atomic number
//          */
//            model::cout<<"  grain "<<grainID<<", rotating"<<defaultColor<<std::endl;
//
//            assert((C2G_in*C2G_in.transpose()-Eigen::Matrix<double,dim,dim>::Identity()).norm()<2.0*DBL_EPSILON*dim*dim && "CRYSTAL TO GLOBAL ROTATION MATRIX IS NOT ORTHOGONAL.");
//            // make sure that C2G is proper
//            assert(std::fabs(C2G_in.determinant()-1.0) < FLT_EPSILON && "C2G IS NOT PROPER.");
//            // store C2G
//            C2G=C2G_in;
//            setLatticeBasis();
//        }

//        /**********************************************************************/
//        std::pair<LatticePlane,LatticePlane> find_confiningPlanes(const LatticeVectorType& sourceL,
//                                                                  const LatticeVectorType& sinkL,
//                                                                  const LatticeVectorType& Burgers) const
//        {
//
//            const LatticeVectorType chord(sinkL-sourceL);
//
//            bool isOnGB=false;
//            const GrainBoundary<dim>* p_GB=NULL;
//            for(const auto& gb : grainBoundaries())
//            {
//                if(   gb.second->latticePlane(grainID).contains(sourceL)
//                   && gb.second->latticePlane(grainID).contains(sinkL)
//                   )
//                {
//                    isOnGB=true;
//                    p_GB=gb.second;
//                    break;
//                }
//            }
//
//
//            // Find the crystallographic planes which may contain chord and Burgers
//            PlaneNormalContainerType allowedGLidePlanes;
//            for (const auto& planeBase : planeNormals()) // Loop over crystal planes
//            {
//                if(planeBase.dot(chord)==0
//                   && planeBase.dot(Burgers)==0)
//                {
//                    allowedGLidePlanes.push_back(planeBase);
//                }
//            }
//            if(isOnGB)
//            {
//                if(   p_GB->latticePlane(grainID).n.dot(chord)==0
//                   && p_GB->latticePlane(grainID).n.dot(Burgers)==0
//                   )
//                {
//                    allowedGLidePlanes.push_back(p_GB->latticePlane(grainID).n);
//                }
//            }
//
//            // Find the crystallographic planes which may contain chord and Burgers
//            PlaneNormalContainerType allowedSessilePlanes;
//            for (const auto& planeBase : planeNormals()) // Loop over crystal planes
//            {
//                if(planeBase.dot(chord)==0)
//                {
//                    allowedSessilePlanes.push_back(planeBase);
//                }
//            }
//            if(isOnGB)
//            {
//                allowedSessilePlanes.push_back(p_GB->latticePlane(grainID).n);
//            }
//
//            // Return the two planes based on the following conditions:
//            if(allowedGLidePlanes.size()) // At least glide plane found
//            {
//                if(isOnGB) // possibly glissile segmento on GB, or sessile at the intersection of GP and GB
//                {
//                    return std::make_pair(LatticePlane(sourceL,allowedGLidePlanes[0]),LatticePlane(p_GB->latticePlane(grainID)));
//                }
//                else // bulk segment, purely glissile
//                {
//                    return std::make_pair(LatticePlane(sourceL,allowedGLidePlanes[0]),LatticePlane(sourceL,allowedGLidePlanes[0]));
//                }
//            }
//            else // No glide planes found. Check for possibly sessile segment
//            {
//                if(allowedSessilePlanes.size()>=2)
//                {
//                    return std::make_pair(LatticePlane(sourceL,allowedSessilePlanes[0]),LatticePlane(sourceL,allowedSessilePlanes[1]));
//                }
//                else
//                {
//                    assert(0 && "SESSILE SEGMENTS MUST FORM ON THE INTERSECTION OF TWO PLANES.");
//                }
//            }
//
//        }

