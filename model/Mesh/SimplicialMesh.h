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


#include <map>
#include <model/Network/Readers/VertexReader.h>

#include <model/Utilities/TerminalColors.h>
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
        void readMesh(const int& meshID)
        {/*!
          */
            
            Simplex<dim,0>::nodeReader.read(meshID,true);
            
            
            this->clear();
            
            std::stringstream filestream;
            filestream << "T/T_" << meshID << ".txt";
            std::string filename(filestream.str());

            
            
            std::ifstream ifs ( filename.c_str() , std::ifstream::in );

            if (ifs.is_open())
            {
                std::cout<<"Reading mesh from file "<<filename<<"..."<<std::flush;

                
                std::string line;

                while (std::getline(ifs, line))
                {
                    std::stringstream ss(line);
                    typename SimplexTraits<dim,dim>::SimplexIDType xID;
                    int i(0);
                    int temp;

                    while (ss >> temp)
                    {
                        xID(i) =temp;
                        i++;
                    }
                    
                    insertSimplex(xID);
                }
                
                std::cout<<" done."<<std::endl;

            }
            else
            {
                std::cout<<"Cannot read mesh file "<<filename<<std::endl;
                assert(0);
            }
        
        
            MeshStats<dim,dim>::stats(true);
        
            Simplex<dim,0>::nodeReader.clear();
        }
        
        
        
        /**********************************************************************/
        void insertSimplex(const typename SimplexTraits<dim,dim>::SimplexIDType& xIN)
        {
            const typename SimplexTraits<dim,dim>::SimplexIDType xID(SimplexTraits<dim,dim>::sortID(xIN));
            this->emplace(xID,xID); // requires gcc4.8 and above
        }
        
        
        
        
	};

    
}	// close namespace
#endif