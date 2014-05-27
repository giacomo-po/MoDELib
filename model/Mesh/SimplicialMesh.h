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


namespace model {
    
    
    /**************************************************************************/
	/**************************************************************************/
	template<int _dim>
	struct SimplicialMesh : public std::map<typename SimplexTraits<_dim,_dim>::SimplexIDType, // key
    /*                                */ const Simplex<_dim,_dim>, // value
    /*                                */ CompareVectorsByComponent<typename SimplexTraits<_dim,_dim>::ScalarIDType,
    /*                                */ SimplexTraits<_dim,_dim>::nVertices> // key compare
    /*                                */ >
    {
        
        enum {dim=_dim};
        
        
        typedef std::map<typename SimplexTraits<dim,dim>::SimplexIDType, // key
        /*            */ const Simplex<dim,dim>, // value
        /*            */ CompareVectorsByComponent<typename SimplexTraits<dim,dim>::ScalarIDType,
        /*                                      */ SimplexTraits<dim,dim>::nVertices> // key compare
        /*            */ >  SimplexMapType;
        
        
        
        
        /**********************************************************************/
        SimplicialMesh()
        {
            std::cout<<greenColor<<"Creating SimplicialMesh"<<defaultColor<<std::endl;
        }
        
        /**********************************************************************/
        SimplicialMesh(const int& meshID)
        {
            std::cout<<greenColor<<"Creating SimplicialMesh"<<defaultColor<<std::endl;
            readMesh(meshID);
        }
        
        /**********************************************************************/
        void readMesh(const int& meshID)
        {/*!
          */
            
            Simplex<dim,0>::nodeReader.read(meshID,true);
            
            this->clear();
            
            //            std::stringstream filestream;
            //            filestream << "T/T_" << meshID << ".txt";
            //            std::string filename(filestream.str());
            
            
            VertexReader<'T',dim+2,size_t> elementReader;
            const bool success=elementReader.read(meshID,true);
            
            //            SequentialBinFile<'T',std::pair<int,typename SimplexTraits<dim,dim>::SimplexIDType>,true> binFile;
            
            if (success)
            {
                double t0(clock());
                std::cout<<"Creating mesh..."<<std::flush;
                for (typename VertexReader<'T',dim+2,size_t>::const_iterator eIter =elementReader.begin();
                     /*                                       */ eIter!=elementReader.end();++eIter)
                {
                    insertSimplex(eIter->second);
                    
                    //                    binFile.write(std::make_pair(eIter->first,eIter->second));
                    
                }
                std::cout<<" done.["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]"<<std::endl;
                MeshStats<dim,dim>::stats(true);
                Simplex<dim,0>::nodeReader.clear();
            }
            else
            {
                std::cout<<"Cannot read mesh file T/T_"<<meshID<<".txt . Mesh is empty."<<std::endl;
            }
            
        }
        
        /**********************************************************************/
        void insertSimplex(const typename SimplexTraits<dim,dim>::SimplexIDType& xIN)
        {
            const typename SimplexTraits<dim,dim>::SimplexIDType xID(SimplexTraits<dim,dim>::sortID(xIN));
            this->emplace(xID,xID); // requires gcc4.8 and above
        }
        
        
        /**********************************************************************/
        std::pair<bool,const Simplex<dim,dim>*> search(const Eigen::Matrix<double,dim,1>& P) const
        {/*!
          *\returns a pair, where:
          * -pair.first is a boolean indicating whether the
          * search succesfully found a Simplex<dim,dim> which includes P.
          * -pair.second is a pointer to the last Simplex<dim,dim> searched.
          */

            return searchWithGuess(P,(const Simplex<dim,dim>*) NULL);
            //            std::set<int> searchSet;
//            std::pair<bool,const Simplex<dim,dim>*> lastSearched(false,NULL);
//            this->begin()->second.search(P,lastSearched,searchSet);
//            return lastSearched;
//            return this->begin()->second.search(P,searchSet);
        }
        
        /**********************************************************************/
        std::pair<bool,const Simplex<dim,dim>*> searchWithGuess(const Eigen::Matrix<double,dim,1>& P, const Simplex<dim,dim>* const guess) const
        {
            
            std::set<int> searchSet;
            std::pair<bool,const Simplex<dim,dim>*> lastSearched(false,guess);
            this->begin()->second.search(P,lastSearched,searchSet);
            return lastSearched;
//            std::set<int> searchSet;
//            return guess->search(P,searchSet);
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
        
	};
    
}	// close namespace
#endif


//            std::ifstream ifs ( filename.c_str() , std::ifstream::in );
//
//            if (ifs.is_open())
//            {
//                std::cout<<"Reading mesh elements from file "<<filename<<"..."<<std::flush;
//
//
//                std::string line;
//
//                while (std::getline(ifs, line))
//                {
//                    std::stringstream ss(line);
//                    typename SimplexTraits<dim,dim>::SimplexIDType xID;
//                    int i(0);
//                    int temp;
//
//                    while (ss >> temp)
//                    {
//                        xID(i) =temp;
//                        i++;
//                    }
//
//                    insertSimplex(xID);
//                }
//
//                std::cout<<" done."<<std::endl;
//                MeshStats<dim,dim>::stats(true);
//                Simplex<dim,0>::nodeReader.clear();
//            }
//            else
//            {
//                std::cout<<"Cannot read mesh file "<<filename<<". Mesh is empty."<<std::endl;
//                //assert(0);
//            }
