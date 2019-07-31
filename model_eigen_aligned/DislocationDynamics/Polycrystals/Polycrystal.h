/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_Polycrystal_H_
#define model_Polycrystal_H_

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
#include <MPIcout.h>
#include <Grain.h>
#include <GrainBoundary.h>
#include <LatticeVector.h>
#include <SequentialOutputFile.h>
#include <StressStraight.h>
#include <GrainBoundaryType.h>
//#include <GlidePlane.h>
#include <TextFileParser.h>
#include <DislocatedMaterial.h>
#include <DislocationMobilityFCC.h>
#include <DislocationMobilityBCC.h>

namespace model
{
    
    
    
    template <int dim>
    class Polycrystal : public  DislocatedMaterial<dim,Isotropic>
    /*               */,private std::map<size_t,Grain<dim>>
    /*               */,private std::map<std::pair<size_t,size_t>,GrainBoundary<dim>>
    /*               */,public  std::vector<StressStraight<dim> >
    {
        typedef DislocatedMaterial<dim,Isotropic> MaterialType;
        typedef SimplicialMesh<dim> SimplicialMeshType;
        typedef MeshRegion<Simplex<dim,dim> > MeshRegionType;
        typedef MeshRegionObserver<MeshRegionType> MeshRegionObserverType;
        typedef LatticeVector<dim> LatticeVectorType;
        typedef ReciprocalLatticeVector<dim> ReciprocalLatticeVectorType;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef Grain<dim> GrainType;
        typedef GrainBoundary<dim> GrainBoundaryType;
        
        /**********************************************************************/
        static std::unique_ptr<DislocationMobilityBase> getMobility(const DislocatedMaterial<dim,Isotropic>& material)
        {
            if(material.crystalStructure=="BCC")
            {
                return std::make_unique<DislocationMobilityBCC>(material);
            }
            else if(material.crystalStructure=="FCC")
            {
                return std::make_unique<DislocationMobilityFCC>(material);
            }
            else if(material.crystalStructure=="HEX")
            {
                std::cout<<"FINISH HERE. HEX MOBILITY NOT IMPLEMENTED YET"<<std::endl;
                return std::make_unique<DislocationMobilityFCC>(material);
            }
            else
            {
                std::cout<<"Unknown mobility for crystal structure '"<<material.crystalStructure<<"'. Exiting."<<std::endl;
                exit(EXIT_FAILURE);
            }
        }
        
        
    public:
        
        const SimplicialMeshType& mesh;
        const std::unique_ptr<DislocationMobilityBase> mobility;
        
        /**********************************************************************/
        Polycrystal(const std::string& polyFile,
                    const SimplicialMeshType& mesh_in) :
        /* init */ MaterialType(TextFileParser(polyFile).readString("materialFile",false))
        /* init */,mesh(mesh_in)
        /* init */,mobility(getMobility(*this))
        {
            model::cout<<greenBoldColor<<"Creating Polycrystal"<<defaultColor<<std::endl;
            TextFileParser polyParser(polyFile);
            
            // Construct Grains
            for(const auto& rIter : mesh.regions())
            {
                grains().emplace(std::piecewise_construct,
                                 std::forward_as_tuple(rIter.second->regionID),
                                 std::forward_as_tuple(*(rIter.second),
                                                       *this,
                                                       polyFile));
            }
            
            // Construct GrainsBoundaries
            grainBoundaryDislocations().clear();
            //            int fileID=1;
            for(const auto& rgnBnd : mesh.regionBoundaries())
            {// loop over region boundaries
                for(const auto& face : rgnBnd.second.faces())
                {// loop over faces of each region boundary
                    grainBoundaries().emplace(std::piecewise_construct,
                                              std::forward_as_tuple(rgnBnd.first),
                                              std::forward_as_tuple(rgnBnd.second,
                                                                    face.second,
                                                                    grain(rgnBnd.first.first),
                                                                    grain(rgnBnd.first.second)
                                                                    )
                                              );
                }
            }
            
            
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
        
        /**********************************************************************/
        GrainType& grain(const size_t& k)
        {
            return grains().at(k);
        }
        
        /**********************************************************************/
        const GrainType& grain(const size_t& k) const
        {
            return grains().at(k);
        }
        
        /**********************************************************************/
        const std::map<size_t,GrainType>& grains() const
        {
            return *this;
        }
        
        /**********************************************************************/
        std::map<size_t,GrainType>& grains()
        {
            return *this;
        }
        
        /**********************************************************************/
        std::map<std::pair<size_t,size_t>,GrainBoundaryType>& grainBoundaries()
        {
            return *this;
        }
        
        /**********************************************************************/
        const std::map<std::pair<size_t,size_t>,GrainBoundaryType>& grainBoundaries() const
        {
            return *this;
        }
        
        /**********************************************************************/
        const GrainBoundaryType& grainBoundary(const size_t& i,
                                               const size_t& j) const
        {
            assert(i!=j && "GrainBoundary IDs cannot be the same.");
            return (i<j)? grainBoundaries().at(std::make_pair(i,j)) : grainBoundaries().at(std::make_pair(j,i));
        }
        
        /**********************************************************************/
        GrainBoundaryType& grainBoundary(const size_t& i,
                                         const size_t& j)
        {
            assert(i!=j && "GrainBoundary IDs cannot be the same.");
            return (i<j)? grainBoundaries().at(std::make_pair(i,j)) : grainBoundaries().at(std::make_pair(j,i));
        }
        
        /**********************************************************************/
        LatticeVectorType latticeVectorFromPosition(const VectorDim& p,
                                                    const Simplex<dim,dim>* const guess) const
        {
            const std::pair<bool,const Simplex<dim,dim>*> temp(mesh.searchWithGuess(p,guess));
            assert(temp.first && "Position not found in mesh");
            return grain(temp.second->region->regionID).latticeVector(p);
        }
        
        /**********************************************************************/
        LatticeVectorType latticeVectorFromPosition(const VectorDim& p) const
        {
            return latticeVectorFromPosition(p,&(mesh.simplices().begin()->second));
        }
        
        /**********************************************************************/
        ReciprocalLatticeVectorType reciprocalLatticeVectorFromPosition(const VectorDim& p,
                                                                        const Simplex<dim,dim>* const guess) const
        {
            const std::pair<bool,const Simplex<dim,dim>*> temp(mesh.searchWithGuess(p,guess));
            assert(temp.first && "Position not found in mesh");
            return grain(temp.second->region->regionID).reciprocalLatticeVector(p);
        }
        
        /**********************************************************************/
        ReciprocalLatticeVectorType reciprocalLatticeVectorFromPosition(const VectorDim& p) const
        {
            return reciprocalLatticeVectorFromPosition(p,&(mesh.simplices().begin()->second));
        }
        
        /**********************************************************************/
        std::vector<StressStraight<dim>>& grainBoundaryDislocations()
        {
            return *this;
        }
        
        /**********************************************************************/
        const std::vector<StressStraight<dim>>& grainBoundaryDislocations() const
        {
            return *this;
        }
        
        /**********************************************************************/
        VectorDim randomPoint() const
        {
            std::mt19937 generator(std::chrono::system_clock::now().time_since_epoch().count());
            std::uniform_real_distribution<double> distribution(0.0,1.0);
            
            VectorDim P0(VectorDim::Zero());
            
            P0 << mesh.xMin(0)+distribution(generator)*(mesh.xMax(0)-mesh.xMin(0)),
            /* */ mesh.xMin(1)+distribution(generator)*(mesh.xMax(1)-mesh.xMin(1)),
            /* */ mesh.xMin(2)+distribution(generator)*(mesh.xMax(2)-mesh.xMin(2));
            
            return P0;
            
        }
        
        /**********************************************************************/
        std::pair<LatticeVector<dim>,int> randomLatticePointInMesh() const
        {
            const VectorDim P0=randomPoint();
            auto searchResult=mesh.search(P0);
            if(searchResult.first)
            {// point inside
                const LatticeVector<dim> L0 = grain(searchResult.second->region->regionID).snapToLattice(P0);
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
        
    };
    
}
#endif

