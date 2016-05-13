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
#include <model/Network/Readers/VertexReader.h>
#include <model/Utilities/TerminalColors.h>
#include <model/Utilities/SequentialBinFile.h>
#include <model/Mesh/SimplexTraits.h>
#include <model/Mesh/Simplex.h>
#include <model/Mesh/SimplexReader.h>
#include <model/Mesh/MeshStats.h>
#include <model/Mesh/MeshRegionObserver.h>
#include <model/MPI/MPIcout.h> // defines mode::cout


namespace model
{
    
    
    /**************************************************************************/
    /**************************************************************************/
    template<int _dim>
    class SimplicialMesh : public SimplexReader<_dim>,
    /*                  */ public std::map<typename SimplexTraits<_dim,_dim>::SimplexIDType, // key
    /*                                */ const Simplex<_dim,_dim>, // value
    /*                                */ CompareVectorsByComponent<typename SimplexTraits<_dim,_dim>::ScalarIDType,
    /*                                */ SimplexTraits<_dim,_dim>::nVertices> // key compare
    /*                                */ >
    {
        
        Eigen::Matrix<double,_dim,1> _xMin;
        Eigen::Matrix<double,_dim,1> _xMax;
        
        
        double vol0;
        
    public:
        
        enum {dim=_dim};
        
        
        typedef std::map<typename SimplexTraits<dim,dim>::SimplexIDType, // key
        /*            */ const Simplex<dim,dim>, // value
        /*            */ CompareVectorsByComponent<typename SimplexTraits<dim,dim>::ScalarIDType,
        /*                                      */ SimplexTraits<dim,dim>::nVertices> // key compare
        /*            */ >  SimplexMapType;
        
        typedef VertexReader<'T',dim+3,size_t> ElementReaderType;

        typedef MeshRegion<Simplex<dim,dim> > MeshRegionType;
        typedef MeshRegionObserver<MeshRegionType> MeshRegionObserverType;

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
        void readMesh(const int& meshID)
        {/*!
          */
            
            vol0=0.0;
            
            model::cout<<greenColor<<"Reading mesh "<<meshID<<defaultColor<<std::endl;
            this->clear();
            
            SimplexReader<dim>::nodeReader.read(meshID,true);
            
            this->clear();

//            VertexReader<'T',dim+2,size_t> elementReader;
            ElementReaderType elementReader; // exaple in 3d: [elementID v1 v2 v3 v4 regionID]
            const bool success=elementReader.read(meshID,true);
            
            //            SequentialBinFile<'T',std::pair<int,typename SimplexTraits<dim,dim>::SimplexIDType>,true> binFile;
            
            if (success)
            {
                const auto t0= std::chrono::system_clock::now();

                model::cout<<"Creating mesh..."<<std::flush;
                for (typename ElementReaderType::const_iterator eIter =elementReader.begin();
                     /*                                       */ eIter!=elementReader.end();++eIter)
                {
//                    insertSimplex(eIter->second);
                    insertSimplex(eIter->second.template segment<dim+1>(0),eIter->second(dim+1));
                    
                    //                    binFile.write(std::make_pair(eIter->first,eIter->second));
                    
                }
                model::cout<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<std::endl;

                MeshStats<dim,dim>::stats();
                SimplexReader<dim>::nodeReader.clear();
                
                if(this->size())
                {
                    _xMin=SimplexObserver<dim,0>::simplexBegin()->second->P0;
                    _xMax=SimplexObserver<dim,0>::simplexBegin()->second->P0;
                    
                    for (typename SimplexObserver<dim,0>::const_iterator nIter = SimplexObserver<dim,0>::simplexBegin();
                         /*                                           */ nIter!= SimplexObserver<dim,0>::simplexEnd();
                         /*                                           */ nIter++)
                    {
                        for(int d=0;d<dim;++d)
                        {
                            if (nIter->second->P0(d)<_xMin(d))
                            {
                                _xMin(d)=nIter->second->P0(d);
                            }
                            if (nIter->second->P0(d)>_xMax(d))
                            {
                                _xMax(d)=nIter->second->P0(d);
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

            model::cout<<"mesh xMin="<<_xMin.transpose()<<std::endl;
            model::cout<<"mesh xMax="<<_xMax.transpose()<<std::endl;
//            model::cout<<"mesh volume="<<volume()<<std::endl;
            
            for(auto rIter : MeshRegionObserverType::regions())
            {
                std::cout<<"mesh region "<<rIter.second->regionID<<" contains "<<rIter.second->size()<<" Simplex<"<<dim<<","<<dim<<">"<<std::endl;
            }
            
        }
        
        /**********************************************************************/
        void insertSimplex(const typename SimplexTraits<dim,dim>::SimplexIDType& xIN,const int& regionID)
        {/*!@param[in] xIN the (unsorted) array of mesh node IDs defining a simplex
          *\brief Inserts a Simplex<dim> into the mesh, with ID equal to the
          * sorted array xIN.
          */
            const typename SimplexTraits<dim,dim>::SimplexIDType xID(SimplexTraits<dim,dim>::sortID(xIN));
            const auto pair=this->emplace(std::piecewise_construct,
                                             std::make_tuple(xID),
                                             std::make_tuple(xID, regionID)
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
            return searchWithGuess(P,&(this->begin()->second));
        }

        /**********************************************************************/
        std::pair<bool,const Simplex<dim,dim>*> searchWithGuess(const Eigen::Matrix<double,dim,1>& P,
                                                                const Simplex<dim,dim>* const guess,
                                                                std::set<const Simplex<dim,dim>*>& searchSet) const
        {/*!@param[in] P position to search for
          * @param[in] guess Simplex* where the search starts
          *\returns a pair, where:
          * -pair.first is a boolean indicating whether the
          * search succesfully found a Simplex<dim,dim> which includes P.
          * -pair.second is a pointer to the last Simplex<dim,dim> searched.
          */

            searchSet.clear();
            std::pair<bool,const Simplex<dim,dim>*> lastSearched(false,NULL);
            guess->convexDelaunaynSearch(P,lastSearched,searchSet);
            
            if(!(lastSearched.first || lastSearched.second->isBoundarySimplex()))
            {
                std::cout<<"P="<<std::setprecision(15)<<std::scientific<<P.transpose()<<std::endl;
                std::cout<<"guess="<<guess->xID<<std::endl;
                std::cout<<"lastSearched="<<lastSearched.second->xID<<std::endl;
                assert(0 && "SEARCH DID NOT END ON BOUNDARY SIMPLEX");
            }
            
            if(!lastSearched.first) // search not successful. Force search neighbors of last searched Simplex
            {
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
                            (*pIter)->convexDelaunaynSearch(P,lastSearched,searchSet);
                            if (lastSearched.first)
                            {
                                break;
                            }
                        }
                    }
                    
                }
                
                if(!lastSearched.first && !lastSearched.second->isBoundarySimplex())
                {
                    lastSearched=temp;
                }
                
            }
            
            if(!(lastSearched.first || lastSearched.second->isBoundarySimplex()))
            {
                std::cout<<"P="<<std::setprecision(15)<<std::scientific<<P.transpose()<<std::endl;
                std::cout<<"guess="<<guess->xID<<std::endl;
                std::cout<<"lastSearched="<<lastSearched.second->xID<<std::endl;
                assert(0 && "SEARCH DID NOT END ON BOUNDARY SIMPLEX");
            }

            return lastSearched;
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
            return searchWithGuess(P,guess,searchSet);
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
        
    };
    
}	// close namespace
#endif
