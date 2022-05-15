/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_Grain_cpp_
#define model_Grain_cpp_

#include <cfloat> // FLT_EPSILON
#include <Eigen/Core>
//#include <RoundEigen.h>
#include <SimplicialMesh.h>
#include <MeshRegionObserver.h>
#include <MeshRegion.h>
#include <LatticeModule.h>
//#include <PeriodicElement.h>
#include <PolycrystallineMaterial.h>
#include <SlipSystem.h>
#include <SingleCrystal.h>
#include <Grain.h>

//#include <BestRationalApproximation.h>

namespace model
{

    template <int dim>
    Grain<dim>::Grain(const MeshRegionType& region_in,
                      const PolycrystallineMaterial<dim,Isotropic>& material,
                      const std::string& polyFile
                      ) :
    /* init */ SingleCrystalType(material,TextFileParser(polyFile).readMatrix<double>("C2G"+std::to_string(region_in.regionID),dim,dim,true),polyFile)
    /* init */,region(region_in)
    /* init */,grainID(region.regionID) // remove grain ID, use lattice.sID
    {
        assert(this->sID==region.regionID); // lattice.sID must be same as grainID
        std::cout<<"  lattice basis="<<this->latticeBasis<<std::endl;
        std::cout<<"  # plane normals="<<this->planeNormals().size()<<std::endl;
        std::cout<<"  # slip systems="<<this->slipSystems().size()<<std::endl;
    }

    template <int dim>
    const typename Grain<dim>::GrainBoundaryContainerType& Grain<dim>::grainBoundaries() const
    {
        return *this;
    }

    template <int dim>
    typename Grain<dim>::GrainBoundaryContainerType& Grain<dim>::grainBoundaries()
    {
        return *this;
    }

    template <int dim>
    std::deque<const LatticePlaneBase*> Grain<dim>::conjugatePlaneNormal(const LatticeVectorType& B,
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

    template class Grain<3>;

}
#endif



