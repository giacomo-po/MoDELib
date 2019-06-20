/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SimplexReader_H_
#define model_SimplexReader_H_

#include <map>
#include <iostream>
#include <Eigen/Dense>
#include <SimplexTraits.h>
//#include <VertexReader.h>
#include <IDreader.h>
#include <GmshReader.h>
#include <SimplexTraits.h>

namespace model
{
    
    
    /**************************************************************************/
    /**************************************************************************/
    template<int dim>
    struct SimplexReader : private std::map<size_t,Eigen::Matrix<double,dim,1>>
    /*                  */,private std::map<size_t,std::pair<typename SimplexTraits<dim,dim>::SimplexIDType,size_t>>
    {
        
        typedef std::map<size_t,Eigen::Matrix<double,dim,1>> NodeContainerType;
        typedef std::map<size_t,std::pair<typename SimplexTraits<dim,dim>::SimplexIDType,size_t>> ElementContainerType;
        
        /**********************************************************************/
        const NodeContainerType& nodes() const
        {
            return *this;
        }
        
        /**********************************************************************/
        NodeContainerType& nodes()
        {
            return *this;
        }
        
        /**********************************************************************/
        const ElementContainerType& elements() const
        {
            return *this;
        }
        
        /**********************************************************************/
        ElementContainerType& elements()
        {
            return *this;
        }
        
//        /**********************************************************************/
//        void readGmsh(const std::string& meshFileName)
//        {
//            std::cout<<"Reading mesh file "<<meshFileName<<std::endl;
//            GmshReader gmsh(meshFileName);
//            this->nodes()=gmsh.nodes();
//            f//or()
//        }
        
        /**********************************************************************/
        void readTN(const int& meshID)
        {
            clear();
            IDreader<'N',1,dim,double> nodeReader;
            nodeReader.read(meshID,true);
            for(const auto& node : nodeReader)
            {
                nodes().emplace(node.first,Eigen::Map<const Eigen::Matrix<double,dim,1>>(node.second.data()));
            }
            assert(nodeReader.size()==nodes().size());
            
            IDreader<'T',1,dim+2,size_t> elementReader;
            elementReader.read(meshID,true);
            for(const auto& eIter : elementReader)
            {
                typename SimplexTraits<dim,dim>::SimplexIDType key;
                for(int d=0;d<dim+1;++d)
                {
                    key[d]=eIter.second[d];
                }
                size_t regionID(eIter.second[dim+1]);
                elements().emplace(eIter.first,std::make_pair(key,regionID));
            }
            assert(elementReader.size()==elements().size());
        }
        
        /**********************************************************************/
        void clear()
        {
            nodes().clear();
            elements().clear();
        }
        
        SimplexReader& simplexReader()
        {
            return *this;
        }
        
        const SimplexReader& simplexReader() const
        {
            return *this;
        }
        
//        static IDreader<'N',1,dim,double> nodeReader;
        
        /**********************************************************************/
        const Eigen::Matrix<double,dim,1>& get_P0(const typename SimplexTraits<dim,0>::SimplexIDType& xID)
        {
            //            const typename VertexReader<'N',dim+1,double>::const_iterator nIter(nodeReader.find((xID)(0)));
            const auto nIter(nodes().find((xID)(0)));
            assert((nIter!=nodes().end()) && "MESH VERTEX NOT FOUND IN N/N_x.txt.");
            return nIter->second;
        }
        
//        /**********************************************************************/
//        static Eigen::Matrix<double,dim,1> get_P0(const typename SimplexTraits<dim,0>::SimplexIDType& xID)
//        {
//            //            const typename VertexReader<'N',dim+1,double>::const_iterator nIter(nodeReader.find((xID)(0)));
//            const typename IDreader<'N',1,dim,double>::const_iterator nIter(nodeReader.find((xID)(0)));
//            assert((nIter!=nodeReader.end()) && "MESH VERTEX NOT FOUND IN N/N_x.txt.");
//            return Eigen::Map<const Eigen::Matrix<double,dim,1>>(nIter->second.data());
//        }
        
    };
    
//    template<int dim>
//    //    VertexReader<'N',dim+1,double> SimplexReader<dim>::nodeReader;
//    IDreader<'N',1,dim,double> SimplexReader<dim>::nodeReader;
    
    
}	// close namespace
#endif
