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
//#include <model/IO/VertexReader.h>
#include <model/Utilities/TerminalColors.h>
#include <model/IO/SequentialBinFile.h>
#include <model/Mesh/SimplexTraits.h>
#include <model/Mesh/Simplex.h>
#include <model/Mesh/SimplexReader.h>
//#include <model/Mesh/MeshStats.h>
#include <model/Mesh/MeshRegionObserver.h>
#include <model/MPI/MPIcout.h> // defines mode::cout
#include <model/Mesh/MeshRegionBoundary.h>


namespace model
{
    
    
    /**************************************************************************/
    /**************************************************************************/
    template<int _dim>
    class SimplicialMesh : public SimplexObserver<_dim>,                                // make sure this is destroyed after map of Simplex<_dim,_dim>
    /*                  */ public MeshRegionObserver<MeshRegion<Simplex<_dim,_dim>>>,   // make sure this is destroyed after map of Simplex<_dim,_dim>
    /*                  */ public SimplexReader<_dim>,
    /*                  */ public std::map<typename SimplexTraits<_dim,_dim>::SimplexIDType, // key
    /*                                */ const Simplex<_dim,_dim>, // value
    /*                                */ CompareVectorsByComponent<typename SimplexTraits<_dim,_dim>::ScalarIDType,
    /*                                */ SimplexTraits<_dim,_dim>::nVertices> // key compare
    /*                                */ >,
    /*                  */ public std::map<std::pair<size_t,size_t>,MeshRegionBoundary<Simplex<_dim,_dim-1>>> // MeshRegionBonudary container
    {
        
        Eigen::Matrix<double,_dim,1> _xMin;
        Eigen::Matrix<double,_dim,1> _xMax;
        
        
        double vol0;
        
    public:
        
        static constexpr int dim=_dim;
        
        
        typedef std::map<typename SimplexTraits<dim,dim>::SimplexIDType, // key
        /*            */ const Simplex<dim,dim>, // value
        /*            */ CompareVectorsByComponent<typename SimplexTraits<dim,dim>::ScalarIDType,
        /*                                      */ SimplexTraits<dim,dim>::nVertices> // key compare
        /*            */ >  SimplexMapType;
        
        //        typedef VertexReader<'T',dim+3,size_t> ElementReaderType;
        typedef IDreader<'T',1,dim+2,size_t> ElementReaderType;
        
        typedef MeshRegion<Simplex<dim,dim> > MeshRegionType;
        typedef MeshRegionObserver<MeshRegionType> MeshRegionObserverType;
        
        typedef MeshRegionBoundary<Simplex<dim,dim-1>> MeshRegionBoundaryType;
        typedef std::pair<size_t,size_t> MeshRegionIDType;
        typedef std::map<MeshRegionIDType,MeshRegionBoundaryType> MeshRegionBoundaryContainerType;
        
        /**********************************************************************/
        SimplicialMesh() :
        /* init list */ _xMin(Eigen::Matrix<double,dim,1>::Zero()),
        /* init list */ _xMax(Eigen::Matrix<double,dim,1>::Zero())
        {
        }
        
        /**********************************************************************/
        SimplicialMesh(const int& meshID) :
        /* init */ vol0(0.0)
        {
            readMesh(meshID);
        }
        
        /**********************************************************************/
        const SimplexMapType& simplices() const
        {
            return *this;
        }
        
        /**********************************************************************/
        SimplexMapType& simplices()
        {
            return *this;
        }
        
        /**********************************************************************/
        void readMesh(const int& meshID)
        {/*!
          */
            
            vol0=0.0;
            
            model::cout<<greenColor<<"Reading mesh "<<meshID<<defaultColor<<std::endl;
            simplices().clear();
            
            SimplexReader<dim>::nodeReader.read(meshID,true);
            
            simplices().clear();
            
            //            VertexReader<'T',dim+2,size_t> elementReader;
            ElementReaderType elementReader; // exaple in 3d: [elementID v1 v2 v3 v4 regionID]
            const bool success=elementReader.read(meshID,true);
            
            //            SequentialBinFile<'T',std::pair<int,typename SimplexTraits<dim,dim>::SimplexIDType>,true> binFile;
            
            if (success)
            {
                const auto t0= std::chrono::system_clock::now();
                
                model::cout<<"Creating mesh..."<<std::flush;
                //                for (typename ElementReaderType::const_iterator eIter =elementReader.begin();
                //                     /*                                       */ eIter!=elementReader.end();++eIter)
                for (const auto& eIter : elementReader)
                {
                    //                    insertSimplex(eIter->second);
                    Eigen::Map<const Eigen::Matrix<size_t,1,dim+2>> row(eIter.second.data());
                    insertSimplex(row.template segment<dim+1>(0),row(dim+1));
                    //                    insertSimplex(eIter->second.template segment<dim+1>(0),eIter->second(dim+1));
                    
                    //                    binFile.write(std::make_pair(eIter->first,eIter->second));
                    
                }
                model::cout<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<std::endl;
                
                this->info(); // print mesh info
                SimplexReader<dim>::nodeReader.clear();
                
                if(simplices().size())
                {
                    _xMin=this->template observer<0>().begin()->second->P0;
                    _xMax=this->template observer<0>().begin()->second->P0;
                    
                    for (const auto& nIter : this->template observer<0>())
                    {
                        for(int d=0;d<dim;++d)
                        {
                            if (nIter.second->P0(d)<_xMin(d))
                            {
                                _xMin(d)=nIter.second->P0(d);
                            }
                            if (nIter.second->P0(d)>_xMax(d))
                            {
                                _xMax(d)=nIter.second->P0(d);
                            }
                        }
                    }
                }
                else
                {
                    model::cout<<"Mesh is empty."<<std::endl;
                }
            }
            else
            {
                //                assert(0);
                model::cout<<"Cannot read mesh file T/T_"<<meshID<<".txt . Mesh is empty."<<std::endl;
            }
            
            // Populate MeshRegionBoundaryContainerType
            regionBoundaries().clear();
            for (const auto& simpl : this->template observer<dim-1>())
            {
                if(simpl.second->isRegionBoundarySimplex())
                {
                    std::set<int> regionIDset=simpl.second->regionIDs();
                    std::pair<size_t,size_t> regionIDs(std::make_pair(*regionIDset.begin(),*regionIDset.rbegin()));
                    const auto regionBndIter=regionBoundaries().find(regionIDs);
                    if(regionBndIter!=regionBoundaries().end())
                    {
                        regionBndIter->second.insert(simpl.second);
                    }
                    else
                    {
                        regionBoundaries().emplace(regionIDs,regionIDs).first->second.insert(simpl.second);
                    }
                    
                }
            }
            
            
            
            
            model::cout<<"mesh xMin="<<_xMin.transpose()<<std::endl;
            model::cout<<"mesh xMax="<<_xMax.transpose()<<std::endl;
            //            model::cout<<"mesh volume="<<volume()<<std::endl;
            
            for(auto rIter : MeshRegionObserverType::regions())
            {
                std::cout<<"mesh region "<<rIter.second->regionID<<" contains "<<rIter.second->simplices().size()<<" Simplex<"<<dim<<","<<dim<<">"<<std::endl;
            }
            
            std::cout<<"Mesh contains "<<regionBoundaries().size()<<" mesh region boundaries"<<std::endl;
            for(const auto& rgnBnd : regionBoundaries())
            {
                std::cout<<"    RegionBoundary ("<<rgnBnd.second.regionBndID.first<<","<<rgnBnd.second.regionBndID.second<<") contains "<<rgnBnd.second.size()<<" Simplex<"<<dim<<","<<dim-1<<">"<<std::endl;
            }
            
        }
        
        /**********************************************************************/
        void insertSimplex(const typename SimplexTraits<dim,dim>::SimplexIDType& xIN,const int& regionID)
        {/*!@param[in] xIN the (unsorted) array of mesh node IDs defining a simplex
          *\brief Inserts a Simplex<dim> into the mesh, with ID equal to the
          * sorted array xIN.
          */
            const typename SimplexTraits<dim,dim>::SimplexIDType xID(SimplexTraits<dim,dim>::sortID(xIN));
            const auto pair=simplices().emplace(std::piecewise_construct,
                                                std::make_tuple(xID),
                                                std::make_tuple(this,xID, regionID)
                                                );
            
            assert(pair.second);
            vol0+=pair.first->second.vol0;
        }
        
        /**********************************************************************/
        std::pair<bool,const Simplex<dim,dim>*> search(const Eigen::Matrix<double,dim,1>& P) const
        {/*!@param[in] P position to search for
          *\returns a pair, where:
          * -pair.first is a boolean indicating whether the
          * search succesfully found a Simplex<dim,dim> which includes P.
          * -pair.second is a pointer to the last Simplex<dim,dim> searched.
          *
          * By default the search starts at this->begin()->second
          */
            return searchWithGuess(P,&(simplices().begin()->second));
        }
        
        /**********************************************************************/
        std::pair<bool,const Simplex<dim,dim>*> searchWithGuess(const Eigen::Matrix<double,dim,1>& P,
                                                                const Simplex<dim,dim>* const guess) const
        {/*!@param[in] P position to search for
          * @param[in] guess Simplex* where the search starts
          *\returns a pair, where:
          * -pair.first is a boolean indicating whether the
          * search succesfully found a Simplex<dim,dim> which includes P.
          * -pair.second is a pointer to the last Simplex<dim,dim> searched.
          */
            
            std::set<const Simplex<dim,dim>*> searchSet;
            return searchWithGuess(true,P,guess,searchSet);
        }
        
        /**********************************************************************/
        std::pair<bool,const Simplex<dim,dim>*> searchRegion(const int& regionID,
                                                             const Eigen::Matrix<double,dim,1>& P) const
        {/*!@param[in] P position to search for
          *\returns a pair, where:
          * -pair.first is a boolean indicating whether the
          * search succesfully found a Simplex<dim,dim> which includes P.
          * -pair.second is a pointer to the last Simplex<dim,dim> searched.
          *
          * By default the search starts at this->begin()->second
          */
            return searchRegionWithGuess(P,*this->region(regionID)->begin());
        }
        
        /**********************************************************************/
        std::pair<bool,const Simplex<dim,dim>*> searchRegionWithGuess(const Eigen::Matrix<double,dim,1>& P,
                                                                      const Simplex<dim,dim>* const guess) const
        {/*!@param[in] P position to search for
          * @param[in] guess Simplex* where the search starts
          *\returns a pair, where:
          * -pair.first is a boolean indicating whether the
          * search succesfully found a Simplex<dim,dim> which includes P.
          * -pair.second is a pointer to the last Simplex<dim,dim> searched.
          */
            std::set<const Simplex<dim,dim>*> searchSet;
            return searchWithGuess(false,P,guess,searchSet);
        }
        
        /**********************************************************************/
        std::pair<bool,const Simplex<dim,dim>*> searchWithGuess(const bool& searchAllRegions,
                                                                const Eigen::Matrix<double,dim,1>& P,
                                                                const Simplex<dim,dim>* const guess,
                                                                std::set<const Simplex<dim,dim>*>& searchSet) const
        {/*!@param[in] P position to search for
          * @param[in] guess Simplex* where the search starts. If searchAllRegions=false, only the region of guess is searched
          *\returns a pair, where:
          * -pair.first is a boolean indicating whether the
          * search succesfully found a Simplex<dim,dim> which includes P.
          * -pair.second is a pointer to the last Simplex<dim,dim> searched.
          */
            
            searchSet.clear();
            std::pair<bool,const Simplex<dim,dim>*> lastSearched(false,NULL);
            guess->convexDelaunaynSearch(searchAllRegions,P,lastSearched,searchSet);
            checkSearch(searchAllRegions,P,guess,lastSearched);
            
            if(!lastSearched.first)
            {// if search was not successful, try to search neighbors of last searched Simplex
                const std::pair<bool,const Simplex<dim,dim>*> temp(lastSearched);
                const Eigen::Matrix<double,dim+1,1> bary=temp.second->pos2bary(P);
                for(int k=0;k<dim+1;++k)
                {
                    //                    if(bary(k)<=0.0)
                    if(bary(k)<=FLT_EPSILON)
                    {
                        for(typename Simplex<dim,dim-1>::ParentContainerType::const_iterator pIter=temp.second->child(k).parentBegin();
                            /*                                                            */ pIter!=temp.second->child(k).parentEnd();++pIter)
                        {
                            if((*pIter)->region->regionID==temp.second->region->regionID || searchAllRegions)
                            {
                                (*pIter)->convexDelaunaynSearch(searchAllRegions,P,lastSearched,searchSet);
                                if (lastSearched.first)
                                {
                                    break;
                                }
                            }
                            
                        }
                    }
                    
                }
                
                if(searchAllRegions)
                {// if search is still unsuccessful, reset lastSearched to temp
                    if(!lastSearched.first && !lastSearched.second->isBoundarySimplex())
                    {
                        lastSearched=temp;
                    }
                }
                else
                {// if search is still unsuccessful, reset lastSearched to temp
                    if(!lastSearched.first && !lastSearched.second->isBoundarySimplex() && !lastSearched.second->isRegionBoundarySimplex())
                    {
                        lastSearched=temp;
                    }
                }
                
            }
            
            checkSearch(searchAllRegions,P,guess,lastSearched);
            
            return lastSearched;
        }
        
        /**********************************************************************/
        void checkSearch(const bool& searchAllRegions,
                         const Eigen::Matrix<double,dim,1>& P,
                         const Simplex<dim,dim>* const guess,
                         const std::pair<bool,const Simplex<dim,dim>*>& lastSearched) const
        {
            if(searchAllRegions)
            {// Check that search was successful, or that it ended on a boundary (point outside)
                if(!(lastSearched.first || lastSearched.second->isBoundarySimplex()))
                {
                    std::cout<<"P="<<std::setprecision(15)<<std::scientific<<P.transpose()<<std::endl;
                    std::cout<<"guess="<<guess->xID<<std::endl;
                    std::cout<<"lastSearched="<<lastSearched.second->xID<<std::endl;
                    assert(0 && "SEARCH DID NOT END ON BOUNDARY SIMPLEX");
                }
            }
            else
            {// Check that search was successful, or that it ended on a boundary or region-boundary
                if(!(lastSearched.first || lastSearched.second->isBoundarySimplex() || lastSearched.second->isRegionBoundarySimplex()))
                {
                    std::cout<<"P="<<std::setprecision(15)<<std::scientific<<P.transpose()<<std::endl;
                    std::cout<<"guess="<<guess->xID<<std::endl;
                    std::cout<<"lastSearched="<<lastSearched.second->xID<<std::endl;
                    assert(0 && "SEARCH DID NOT END ON BOUNDARY/REGION-BOUNDARY SIMPLEX");
                }
            }
        }
        
        /**********************************************************************/
        std::pair<bool,const Simplex<dim,dim>*> isStrictlyInsideMesh(const Eigen::Matrix<double,dim,1>& P,
                                                                     const Simplex<dim,dim>* const guess,
                                                                     const double& tol) const
        {
            std::pair<bool,const Simplex<dim,dim>*> temp(searchWithGuess(P,guess));
            if(temp.first)
            {
                const Eigen::Matrix<double,dim+1,1> bary(temp.second->pos2bary(P));
                int kMin;
                const double baryMin(bary.minCoeff(&kMin));
                if (std::fabs(baryMin)<tol && guess->child(kMin).isBoundarySimplex()) // on a boundary face
                {
                    temp.first=false;
                }
            }
            return temp;
        }
        
        /**********************************************************************/
        std::pair<bool,const Simplex<dim,dim>*> isOnMeshBoundary(const Eigen::Matrix<double,dim,1>& P, const Simplex<dim,dim>* const guess, const double& tol) const
        {
            std::pair<bool,const Simplex<dim,dim>*> temp(searchWithGuess(P,guess));
            if(temp.first)
            {
                const Eigen::Matrix<double,dim+1,1> bary(temp.second->pos2bary(P));
                int kMin;
                const double baryMin(bary.minCoeff(&kMin));
                if (!(std::fabs(baryMin)<tol && guess->child(kMin).isBoundarySimplex())) // not on a boundary face
                {
                    temp.first=false;
                }
            }
            return temp;
        }
        
        /**********************************************************************/
        const Eigen::Matrix<double,dim,1>& xMin() const
        {
            return _xMin;
        }
        
        /**********************************************************************/
        const double& xMin(const int& k) const
        {
            return _xMin(k);
        }
        
        /**********************************************************************/
        const Eigen::Matrix<double,dim,1>& xMax() const
        {
            return _xMax;
        }
        
        /**********************************************************************/
        const double& xMax(const int& k) const
        {
            return _xMax(k);
        }
        
        /**********************************************************************/
        const double& volume() const
        {
            return vol0;
        }
        
        /**********************************************************************/
        const MeshRegionBoundaryContainerType& regionBoundaries() const
        {
            return *this;
        }
        
        /**********************************************************************/
        MeshRegionBoundaryContainerType& regionBoundaries()
        {
            return *this;
        }
        
    };
    
}	// close namespace
#endif
