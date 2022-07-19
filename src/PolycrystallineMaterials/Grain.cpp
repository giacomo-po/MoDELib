/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_Grain_cpp_
#define model_Grain_cpp_

#include <Grain.h>
#include <BCClattice.h>
#include <FCClattice.h>
//#include <HEXlattice.h>

//#include <BestRationalApproximation.h>

namespace model
{

    template <int dim>
    Grain<dim>::Grain(const MeshRegionType& region_in,
                      const PolycrystallineMaterialBase& material,
                      const std::string& polyFile
                      ) :
//    /* init */ SingleCrystalType(material,TextFileParser(polyFile).readMatrix<double>("C2G"+std::to_string(region_in.regionID),dim,dim,true),polyFile)
    /* init */ region(region_in)
    /* init */,grainID(region.regionID) // remove grain ID, use lattice.sID
    /* init */,singleCrystal(getSingleCrystal(region_in,material,polyFile))
    {
        assert(singleCrystal->sID==region.regionID); // lattice.sID must be same as grainID
        std::cout<<"  latticeBasis="<<singleCrystal->latticeBasis<<std::endl;
        std::cout<<"  # planeNormals="<<singleCrystal->planeNormals().size()<<std::endl;
        std::cout<<"  # slipSystems="<<singleCrystal->slipSystems().size()<<std::endl;
        std::cout<<"  # secondPhases="<<singleCrystal->secondPhases().size()<<std::endl;
    }

    template <int dim>
    std::shared_ptr<SingleCrystalBase<dim>> Grain<dim>::getSingleCrystal(const MeshRegionType& region_in,
                                                                               const PolycrystallineMaterialBase& material,
                                                                               const std::string& polyFile)
    {
        
        const MatrixDimD C2G(TextFileParser(polyFile).readMatrix<double>("C2G"+std::to_string(region_in.regionID),dim,dim,true));
        
        if(material.crystalStructure=="BCC")
        {
            return std::shared_ptr<SingleCrystalBase<dim>>(new BCClattice<dim>(C2G,material,polyFile));
        }
        else if(material.crystalStructure=="FCC")
        {
            return std::shared_ptr<SingleCrystalBase<dim>>(new FCClattice<dim>(C2G,material,polyFile));
        }
//        else if(material.crystalStructure=="HEX")
//        {
//            return std::shared_ptr<SingleCrystalBase<dim>>(new HEXlattice<dim>(C2G,material,polyFile));
//        }
        else
        {
            throw std::runtime_error("Grain::getSingleCrystal: unknown crystal structure "+material.crystalStructure);
            return nullptr;
        }
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
            for (const auto& planeNormal : singleCrystal->planeNormals())
            {
                if(	 B.dot(*planeNormal)==0 && N.cross(*planeNormal).squaredNorm()>0)
                {
                    temp.push_back(planeNormal.get());
                }
            }
        }
        return temp;
    }

    template class Grain<3>;

}
#endif



