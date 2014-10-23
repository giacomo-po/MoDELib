/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SimplicialMesh_H_
#define model_SimplicialMesh_H_

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
#include <model/Mesh/MeshStats.h>
#include <model/MPI/MPIcout.h> // defines mode::cout


namespace model {
    
    
    /**************************************************************************/
	/**************************************************************************/
	template<int _dim>
	class SimplicialMesh : public std::map<typename SimplexTraits<_dim,_dim>::SimplexIDType, // key
    /*                                */ const Simplex<_dim,_dim>, // value
    /*                                */ CompareVectorsByComponent<typename SimplexTraits<_dim,_dim>::ScalarIDType,
    /*                                */ SimplexTraits<_dim,_dim>::nVertices> // key compare
    /*                                */ >
    {
        
        Eigen::Matrix<double,_dim,1> _xMin;
        Eigen::Matrix<double,_dim,1> _xMax;
        
    public:
        
        enum {dim=_dim};
        
        
        typedef std::map<typename SimplexTraits<dim,dim>::SimplexIDType, // key
        /*            */ const Simplex<dim,dim>, // value
        /*            */ CompareVectorsByComponent<typename SimplexTraits<dim,dim>::ScalarIDType,
        /*                                      */ SimplexTraits<dim,dim>::nVertices> // key compare
        /*            */ >  SimplexMapType;
        
        /**********************************************************************/
        SimplicialMesh() :
        /* init list */ _xMin(Eigen::Matrix<double,dim,1>::Constant( DBL_MAX)),
        /* init list */ _xMax(Eigen::Matrix<double,dim,1>::Constant(-DBL_MAX))

        {
        }
        
        /**********************************************************************/
        SimplicialMesh(const int& meshID)
        {
            readMesh(meshID);
        }
        
        /**********************************************************************/
        void readMesh(const int& meshID)
        {/*!
          */
            
            model::cout<<greenColor<<"Reading mesh "<<meshID<<defaultColor<<std::endl;
            this->clear();
            
            Simplex<dim,0>::nodeReader.read(meshID,true);
            
            this->clear();
            
            VertexReader<'T',dim+2,size_t> elementReader;
            const bool success=elementReader.read(meshID,true);
            
            //            SequentialBinFile<'T',std::pair<int,typename SimplexTraits<dim,dim>::SimplexIDType>,true> binFile;
            
            if (success)
            {
                double t0(clock());
                model::cout<<"Creating mesh..."<<std::flush;
                for (typename VertexReader<'T',dim+2,size_t>::const_iterator eIter =elementReader.begin();
                     /*                                       */ eIter!=elementReader.end();++eIter)
                {
                    insertSimplex(eIter->second);
                    
                    //                    binFile.write(std::make_pair(eIter->first,eIter->second));
                    
                }
                model::cout<<" done.["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]"<<std::endl;
                MeshStats<dim,dim>::stats(true);
                Simplex<dim,0>::nodeReader.clear();
            }
            else
            {
                model::cout<<"Cannot read mesh file T/T_"<<meshID<<".txt . Mesh is empty."<<std::endl;
            }
            
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
        
        /**********************************************************************/
        void insertSimplex(const typename SimplexTraits<dim,dim>::SimplexIDType& xIN)
        {/*!@param[in] xIN the (unsorted) array of mesh node IDs defining a simplex
          *\brief Inserts a Simplex<dim> into the mesh, with ID equal to the 
          * sorted array xIN.
          */
            const typename SimplexTraits<dim,dim>::SimplexIDType xID(SimplexTraits<dim,dim>::sortID(xIN));
            this->emplace(xID,xID); // requires gcc4.8 and above
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
        std::pair<bool,const Simplex<dim,dim>*> searchWithGuess(const Eigen::Matrix<double,dim,1>& P, const Simplex<dim,dim>* const guess) const
        {/*!@param[in] P position to search for
          * @param[in] guess Simplex* where the search starts
          *\returns a pair, where:
          * -pair.first is a boolean indicating whether the
          * search succesfully found a Simplex<dim,dim> which includes P.
          * -pair.second is a pointer to the last Simplex<dim,dim> searched.
          */
            std::set<int> searchSet;
            std::pair<bool,const Simplex<dim,dim>*> lastSearched(false,NULL);
            guess->convexDelaunaynSearch(P,lastSearched,searchSet);
            return lastSearched;
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
        const Eigen::Matrix<double,dim,1>& xMax() const
        {
            return _xMax;
        }
        
	};
    
}	// close namespace
#endif
