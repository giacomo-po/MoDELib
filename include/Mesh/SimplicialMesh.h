/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SimplicialMesh_H_
#define model_SimplicialMesh_H_

#include <chrono>
#include <sstream>
#include <fstream>
#include <assert.h>
#include <utility>      // std::pair, std::make_pair
#include <set>

#include <map>
//#include <VertexReader.h>
#include <TerminalColors.h>
#include <SimplexTraits.h>
#include <Simplex.h>
#include <SimplexReader.h>
//#include <MeshStats.h>
#include <MeshRegionObserver.h>
 // defines mode::cout
#include <MeshRegionBoundary.h>
#include <SimplexObserver.h>
//#include <SimplicialMeshFace.h>
#include <GmshReader.h>


namespace model
{
    
    
    /**************************************************************************/
    /**************************************************************************/
    template<int _dim>
    class SimplicialMesh : public SimplexObserver<_dim>                                // make sure this is destroyed after map of Simplex<_dim,_dim>
    /*                  */,public MeshRegionObserver<MeshRegion<_dim>>   // make sure this is destroyed after map of Simplex<_dim,_dim>
    /*                  */,public SimplexReader<_dim>
    /*                  */,public std::map<typename SimplexTraits<_dim,_dim>::SimplexIDType, // key
    /*                                */ const Simplex<_dim,_dim>>
    /*                  */,public std::map<std::pair<size_t,size_t>,MeshRegionBoundary<_dim>> // MeshRegionBoundary container
    //    /*                  */,public std::deque<SimplicialMeshFace<_dim>> // MeshRegionBoundary container
    {
     
        typedef Eigen::Matrix<double,_dim,1> VectorDim;
        
        VectorDim _xMin;
        VectorDim _xMax;
        
        double vol0;
        
        void createMesh(const std::set<int>&);
        
        
    public:
        
        
        static constexpr int dim=_dim;
        
        
        typedef std::map<typename SimplexTraits<dim,dim>::SimplexIDType, // key
        /*            */ const Simplex<dim,dim>>  SimplexMapType;
        typedef IDreader<1,dim+2,size_t> ElementReaderType;
        typedef MeshRegion<dim> MeshRegionType;
        typedef MeshRegionObserver<MeshRegionType> MeshRegionObserverType;
        typedef MeshRegionBoundary<dim> MeshRegionBoundaryType;
        typedef std::pair<size_t,size_t> MeshRegionIDType;
        typedef std::map<MeshRegionIDType,MeshRegionBoundaryType> MeshRegionBoundaryContainerType;
        
        SimplicialMesh();
        
        SimplicialMesh(const std::string& meshFileName,const Eigen::Matrix<double,dim,dim>& A,const Eigen::Matrix<double,dim,1>& x0,const std::set<int>&);
        
        void readMesh(const std::string& meshFileName,const Eigen::Matrix<double,dim,dim>& A,const Eigen::Matrix<double,dim,1>& x0,const std::set<int>&);
        
        const SimplexMapType& simplices() const;
        
        SimplexMapType& simplices();
        
        void updateRegions();
        
        void updateRegionBoundaries();
        
        void identifyParallelFaces(const std::set<int>&);
        
        std::vector<VectorDim> periodicBasis() const;
        
        std::vector<VectorDim> periodicShifts(const std::vector<int>& periodicImageSize) const;

        void insertSimplex(const typename SimplexTraits<dim,dim>::SimplexIDType& xIN,const int& regionID);
        
        std::pair<bool,const Simplex<_dim,_dim>*> search(const Eigen::Matrix<double,dim,1>& P) const;
        
        std::pair<bool,const Simplex<_dim,_dim>*> searchWithGuess(const Eigen::Matrix<double,dim,1>& P,
                                                                const Simplex<dim,dim>* const guess) const;
        
        std::pair<bool,const Simplex<_dim,_dim>*> searchRegion(const int& regionID,
                                                             const Eigen::Matrix<double,dim,1>& P) const;
        
        std::pair<bool,const Simplex<_dim,_dim>*> searchRegionWithGuess(const Eigen::Matrix<double,dim,1>& P,
                                                                      const Simplex<dim,dim>* const guess) const;
        
        std::pair<bool,const Simplex<_dim,_dim>*> searchWithGuess(const bool& searchAllRegions,
                                                                const Eigen::Matrix<double,dim,1>& P,
                                                                const Simplex<dim,dim>* const guess,
                                                                std::set<const Simplex<dim,dim>*>& searchSet) const;
        
        
        void checkSearch(const bool& searchAllRegions,
                         const Eigen::Matrix<double,dim,1>& P,
                         const Simplex<dim,dim>* const guess,
                         const std::pair<bool,const Simplex<dim,dim>*>& lastSearched) const;
        
        
        std::pair<bool,const Simplex<_dim,_dim>*> isStrictlyInsideMesh(const Eigen::Matrix<double,dim,1>& P,
                                                                     const Simplex<dim,dim>* const guess,
                                                                     const double& tol) const;
        
        std::pair<bool,const Simplex<_dim,_dim>*> isOnMeshBoundary(const Eigen::Matrix<double,dim,1>& P, const Simplex<dim,dim>* const guess, const double& tol) const;
        
        const Eigen::Matrix<double,_dim,1>& xMin() const;
        
        const double& xMin(const int& k) const;
        
        const Eigen::Matrix<double,_dim,1>& xMax() const;
        
        const double& xMax(const int& k) const;
        
        const double& volume() const;
        
        const MeshRegionBoundaryContainerType& regionBoundaries() const;
        
        MeshRegionBoundaryContainerType& regionBoundaries();
        
        const MeshRegionBoundaryType& regionBoundary(const int& i,const int& j) const;
        
        
        
        
    };
    
}
#endif
