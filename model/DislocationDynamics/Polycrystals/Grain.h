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
        const size_t& grainID;
        
        /**********************************************************************/
        Grain(const MeshRegionType& region_in,
              const DislocatedMaterial<dim,Isotropic>& material,
              const std::string& polyFile
              ) :
        /* init */ SingleCrystalType(material,TextFileParser(polyFile).readMatrix<double>("C2G"+std::to_string(region_in.regionID),dim,dim,true))
        /* init */,region(region_in)
        /* init */,grainID(region.regionID) // remove grain ID, use lattice.sID
        {
//            std::cout<<"region.regionID="<<region.regionID<<std::endl;
//            std::cout<<"this->sID="<<this->sID<<std::endl;
//            std::cout<<"grainID="<<grainID<<std::endl;
            assert(this->sID==region.regionID); // lattice.sID must be same as grainID
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
