/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_Polycrystal_H_
#define model_Polycrystal_H_

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

namespace model
{
    
    
    
    template <int dim>
    class Polycrystal : public  PolycrystallineMaterialBase
    {
        typedef SimplicialMesh<dim> SimplicialMeshType;
        typedef MeshRegion<dim> MeshRegionType;
        typedef MeshRegionObserver<MeshRegionType> MeshRegionObserverType;
        typedef LatticeVector<dim> LatticeVectorType;
        typedef ReciprocalLatticeVector<dim> ReciprocalLatticeVectorType;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef Grain<dim> GrainType;
        typedef GrainBoundary<dim> GrainBoundaryType;
        
        std::map<size_t,GrainType> getGrains(const std::string& polyFile) const;
        std::map<std::pair<size_t,size_t>,const GrainBoundaryType* const> getGrainBoundaries() const;

//        const std::map<size_t,GrainType>& grains() const;

        double getAtomicVolume() const;
                
    public:
        
        const SimplicialMeshType& mesh;
        const std::map<size_t,Grain<dim>> grains;
        const std::map<std::pair<size_t,size_t>,const GrainBoundaryType* const> grainBoundaries;
        const double Omega;
        
        Polycrystal(const std::string& polyFile,const SimplicialMeshType& mesh_in);
//        GrainType& grain(const size_t& k);
        const GrainType& grain(const size_t& k) const;
//        std::map<std::pair<size_t,size_t>,GrainBoundaryType>& grainBoundaries();
//        const std::map<std::pair<size_t,size_t>,GrainBoundaryType>& grainBoundaries() const;
        const GrainBoundaryType& grainBoundary(const size_t& i,const size_t& j) const;
//        GrainBoundaryType& grainBoundary(const size_t& i,const size_t& j);
        LatticeVectorType latticeVectorFromPosition(const VectorDim& p, const Simplex<dim,dim>* const guess) const;
        LatticeVectorType latticeVectorFromPosition(const VectorDim& p) const;
        ReciprocalLatticeVectorType reciprocalLatticeVectorFromPosition(const VectorDim& p,const Simplex<dim,dim>* const guess) const;
        ReciprocalLatticeVectorType reciprocalLatticeVectorFromPosition(const VectorDim& p) const;
        VectorDim randomPoint() const;
        std::pair<LatticeVector<dim>,int> randomLatticePointInMesh() const;
        
    };
    
}
#endif

