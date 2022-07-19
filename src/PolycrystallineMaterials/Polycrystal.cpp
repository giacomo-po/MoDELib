/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_Polycrystal_cpp_
#define model_Polycrystal_cpp_

#include <filesystem>
#include <string>       // std::string
#include <iostream>     // std::cout
#include <sstream>
#include <utility>
#include <tuple>
#include <map>
#include <vector>
#include <deque>
#include <tuple>
#include <chrono>
#include <random>
#include <memory>
#include <Eigen/Core>
#include <SimplicialMesh.h>

#include <Grain.h>
#include <GrainBoundary.h>
#include <LatticeVector.h>
//#include <StressStraight.h>
//#include <GrainBoundaryType.h>
//#include <GlidePlane.h>
#include <TextFileParser.h>
//#include <PolycrystallineMaterial.h>
//#include <DislocationMobilityFCC.h>
//#include <DislocationMobilityBCC.h>
#include <Polycrystal.h>

namespace model
{

    template <int dim>
    Polycrystal<dim>::Polycrystal(const std::string& polyFile,
                                  const SimplicialMeshType& mesh_in) :
    /* init */ PolycrystallineMaterialBase(std::filesystem::path(polyFile).parent_path().string()+"/"+TextFileParser(polyFile).readString("materialFile",false),
                            TextFileParser(polyFile).readScalar<double>("absoluteTemperature",true))
    /* init */,mesh(mesh_in)
    /* init */,grains(getGrains(polyFile))
    /* init */,grainBoundaries(getGrainBoundaries())
    /* init */,Omega(getAtomicVolume())
    {
        std::cout<<greenBoldColor<<"Created Polycrystal"<<defaultColor<<std::endl;
        
        // Construct Grains
//        for(const auto& rIter : mesh.regions())
//        {
//            std::cout<<greenBoldColor<<"Creating Grain "<<rIter.second->regionID<<defaultColor<<std::endl;
//
//            //                const auto C2G(TextFileParser(polyFolder+"/polyCrystal.txt").readMatrix<double>("C2G"+std::to_string(rIter.second->regionID),dim,dim,true));
//
//            StaticID<Lattice<dim>>::set_count(rIter.second->regionID);
//            grains().emplace(std::piecewise_construct,
//                             std::forward_as_tuple(rIter.second->regionID),
//                             std::forward_as_tuple(*(rIter.second),
//                                                   *this,
//                                                   polyFile));
//        }
        
        // Construct GrainsBoundaries
        //            grainBoundaryDislocations().clear();
        //            int fileID=1;
//        for(const auto& rgnBnd : mesh.regionBoundaries())
//        {// loop over region boundaries
//            for(const auto& face : rgnBnd.second.faces())
//            {// loop over faces of each region boundary
//                grainBoundaries.emplace(std::piecewise_construct,
//                                          std::forward_as_tuple(rgnBnd.first),
//                                          std::forward_as_tuple(rgnBnd.second,
//                                                                face.second,
//                                                                grain(rgnBnd.first.first),
//                                                                grain(rgnBnd.first.second)
//                                                                )
//                                          );
//            }
//        }
        
        
        //            if(grainBoundaryDislocations().size())
        //            {
        //                model::SequentialOutputFile<'B',1>::set_count(0);
        //                model::SequentialOutputFile<'B',1>::set_increment(1);
        //                model::SequentialOutputFile<'B',true> stressStraightFile;
        //                size_t n=0;
        //                for(const auto& sStraight : grainBoundaryDislocations())
        //                {
        //                    stressStraightFile<<n<<"\t"<<sStraight.P0.transpose()<<"\t"<<sStraight.P1.transpose()<<"\t"<<sStraight.b.transpose()<<std::endl;
        //                    n++;
        //                }
        //            }
        
    }


template <int dim>
std::map<size_t,typename Polycrystal<dim>::GrainType> Polycrystal<dim>::getGrains(const std::string& polyFile) const
{
    std::map<size_t,typename Polycrystal<dim>::GrainType> temp;
    for(const auto& rIter : mesh.regions())
    {
        std::cout<<greenBoldColor<<"Creating Grain "<<rIter.second->regionID<<defaultColor<<std::endl;
        
        //                const auto C2G(TextFileParser(polyFolder+"/polyCrystal.txt").readMatrix<double>("C2G"+std::to_string(rIter.second->regionID),dim,dim,true));
        
        StaticID<Lattice<dim>>::set_count(rIter.second->regionID);
        temp.emplace(std::piecewise_construct,
                         std::forward_as_tuple(rIter.second->regionID),
                         std::forward_as_tuple(*(rIter.second),
                                               *this,
                                               polyFile));
    }
    
    
    
    for(const auto& rgnBnd : mesh.regionBoundaries())
    {// loop over region boundaries
        for(const auto& face : rgnBnd.second.faces())
        {// loop over faces of each region boundary
            
            
            std::shared_ptr<GrainBoundary<dim>> gb(new GrainBoundary<dim>(rgnBnd.second,
                                                                          face.second,
                                                                          grain(rgnBnd.first.first),
                                                                          grain(rgnBnd.first.second)));
                                                   
            temp.at(rgnBnd.first.first ).grainBoundaries().emplace(rgnBnd.first,gb);
            temp.at(rgnBnd.first.second).grainBoundaries().emplace(rgnBnd.first,gb);
//
//            temp.emplace(std::piecewise_construct,
//                                      std::forward_as_tuple(rgnBnd.first),
//                                      std::forward_as_tuple(rgnBnd.second,
//                                                            face.second,
//                                                            grain(rgnBnd.first.first),
//                                                            grain(rgnBnd.first.second)
//                                                            )
//                                      );
        }
    }
    
    return temp;
}

template <int dim>
double Polycrystal<dim>::getAtomicVolume() const
{
    const double omg0(grains.begin()->second.singleCrystal->latticeBasis.determinant());
    for(const auto& grain : grains)
    {
        const double omg(grain.second.singleCrystal->latticeBasis.determinant());
        if(std::fabs(omg-omg0)>FLT_EPSILON)
        {
            throw std::runtime_error("grains have different atomic volume");
        }
    }
    
    return omg0;
}

template <int dim>
std::map<std::pair<size_t,size_t>,const GrainBoundary<dim>* const> Polycrystal<dim>::getGrainBoundaries() const
{
    std::map<std::pair<size_t,size_t>,const GrainBoundaryType* const> temp;
//    for(const auto& rgnBnd : mesh.regionBoundaries())
//    {// loop over region boundaries
//        for(const auto& face : rgnBnd.second.faces())
//        {// loop over faces of each region boundary
//            temp.emplace(std::piecewise_construct,
//                                      std::forward_as_tuple(rgnBnd.first),
//                                      std::forward_as_tuple(rgnBnd.second,
//                                                            face.second,
//                                                            grain(rgnBnd.first.first),
//                                                            grain(rgnBnd.first.second)
//                                                            )
//                                      );
//        }
//    }
    
    for(const auto& grain : grains)
    {
        for(const auto& gb : grain.second.grainBoundaries())
        {
            temp.emplace(gb.first,gb.second.get());
        }
    }
    
    
    return temp;
}



//    template <int dim>
//    typename Polycrystal<dim>::GrainType& Polycrystal<dim>::grain(const size_t& k)
//    {
//        return grains.at(k);
//    }

    template <int dim>
    const typename Polycrystal<dim>::GrainType& Polycrystal<dim>::grain(const size_t& k) const
    {
        return grains.at(k);
    }



//    template <int dim>
//    std::map<size_t,typename Polycrystal<dim>::GrainType>& Polycrystal<dim>::grains()
//    {
//        return *this;
//    }

//    template <int dim>
//    std::map<std::pair<size_t,size_t>,typename Polycrystal<dim>::GrainBoundaryType>& Polycrystal<dim>::grainBoundaries()
//    {
//        return *this;
//    }
//
//    template <int dim>
//    const std::map<std::pair<size_t,size_t>,typename Polycrystal<dim>::GrainBoundaryType>& Polycrystal<dim>::grainBoundaries() const
//    {
//        return *this;
//    }

    template <int dim>
    const typename Polycrystal<dim>::GrainBoundaryType& Polycrystal<dim>::grainBoundary(const size_t& i,
                                                                                        const size_t& j) const
    {
//        assert(i!=j && "GrainBoundary IDs cannot be the same.");
        return (i<j)? *grainBoundaries.at(std::make_pair(i,j)) : *grainBoundaries.at(std::make_pair(j,i));
    }

//    template <int dim>
//    typename Polycrystal<dim>::GrainBoundaryType& Polycrystal<dim>::grainBoundary(const size_t& i,
//                                                                                  const size_t& j)
//    {
//        assert(i!=j && "GrainBoundary IDs cannot be the same.");
//        return (i<j)? grainBoundaries.at(std::make_pair(i,j)) : grainBoundaries.at(std::make_pair(j,i));
//    }

    template <int dim>
    typename Polycrystal<dim>::LatticeVectorType Polycrystal<dim>::latticeVectorFromPosition(const VectorDim& p,
                                                                                             const Simplex<dim,dim>* const guess) const
    {
        const std::pair<bool,const Simplex<dim,dim>*> temp(mesh.searchWithGuess(p,guess));
        assert(temp.first && "Position not found in mesh");
        return grain(temp.second->region->regionID).singleCrystal->latticeVector(p);
    }

    template <int dim>
    typename Polycrystal<dim>::LatticeVectorType Polycrystal<dim>::latticeVectorFromPosition(const VectorDim& p) const
    {
        return latticeVectorFromPosition(p,&(mesh.simplices().begin()->second));
    }

    template <int dim>
    typename Polycrystal<dim>::ReciprocalLatticeVectorType Polycrystal<dim>::reciprocalLatticeVectorFromPosition(const VectorDim& p,
                                                                                                                 const Simplex<dim,dim>* const guess) const
    {
        const std::pair<bool,const Simplex<dim,dim>*> temp(mesh.searchWithGuess(p,guess));
        assert(temp.first && "Position not found in mesh");
        return grain(temp.second->region->regionID).singleCrystal->reciprocalLatticeVector(p);
    }

    template <int dim>
    typename Polycrystal<dim>::ReciprocalLatticeVectorType Polycrystal<dim>::reciprocalLatticeVectorFromPosition(const VectorDim& p) const
    {
        return reciprocalLatticeVectorFromPosition(p,&(mesh.simplices().begin()->second));
    }


    template <int dim>
    typename Polycrystal<dim>::VectorDim Polycrystal<dim>::randomPoint() const
    {
        std::mt19937 generator(std::chrono::system_clock::now().time_since_epoch().count());
        std::uniform_real_distribution<double> distribution(0.0,1.0);
        
        VectorDim P0(VectorDim::Zero());
        
        P0 << mesh.xMin(0)+distribution(generator)*(mesh.xMax(0)-mesh.xMin(0)),
        /* */ mesh.xMin(1)+distribution(generator)*(mesh.xMax(1)-mesh.xMin(1)),
        /* */ mesh.xMin(2)+distribution(generator)*(mesh.xMax(2)-mesh.xMin(2));
        
        return P0;
        
    }

    template <int dim>
    std::pair<LatticeVector<dim>,int> Polycrystal<dim>::randomLatticePointInMesh() const
    {
        const VectorDim P0=randomPoint();
        auto searchResult=mesh.search(P0);
        if(searchResult.first)
        {// point inside
            const LatticeVector<dim> L0 = grain(searchResult.second->region->regionID).singleCrystal->snapToLattice(P0);
            searchResult=mesh.searchRegionWithGuess(L0.cartesian(),searchResult.second);
            if(searchResult.first)
            {// point inside
                return std::make_pair(L0,searchResult.second->region->regionID);
            }
            else
            {
                return randomLatticePointInMesh();
            }
        }
        else
        {
            return randomLatticePointInMesh();
        }
        
    }

template class Polycrystal<3>;

}
#endif

