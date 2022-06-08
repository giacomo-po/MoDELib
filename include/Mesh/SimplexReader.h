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
    struct SimplexReaderBase : private std::map<size_t,Eigen::Matrix<double,dim,1>>
    /*                      */,private std::map<size_t,std::pair<typename SimplexTraits<dim,dim>::SimplexIDType,size_t>>
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
        
        /**********************************************************************/
        void clear()
        {
            nodes().clear();
            elements().clear();
        }
        
        SimplexReaderBase<dim>& simplexReader()
        {
            return *this;
        }
        
        const SimplexReaderBase<dim>& simplexReader() const
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
    };
    
    /**************************************************************************/
    /**************************************************************************/
    template<int dim>
    struct SimplexReader : public SimplexReaderBase<dim>
    {
        
        /**********************************************************************/
        void read(const int& meshID)
        {
            
            this->clear();
            IDreader<1,dim,double> nodeReader("N/N");
            nodeReader.read(meshID,true);
            for(const auto& node : nodeReader)
            {
                this->nodes().emplace(node.first,Eigen::Map<const Eigen::Matrix<double,dim,1>>(node.second.data()));
            }
            assert(nodeReader.size()==this->nodes().size());
            
            IDreader<1,dim+2,size_t> elementReader("T/T");
            elementReader.read(meshID,true);
            for(const auto& eIter : elementReader)
            {
                typename SimplexTraits<dim,dim>::SimplexIDType key;
                for(int d=0;d<dim+1;++d)
                {
                    key[d]=eIter.second[d];
                }
                size_t regionID(eIter.second[dim+1]);
                this->elements().emplace(eIter.first,std::make_pair(key,regionID));
            }
            //            assert(elementReader.size()==this->elements().size());
        }
        
    };
    
    template<>
    struct SimplexReader<3> : public SimplexReaderBase<3>
    {
        
        static constexpr int dim=3;
        
        /**********************************************************************/
        void read(const std::string& meshFileName,const Eigen::Matrix<double,3,3>& A,const Eigen::Matrix<double,3,1>& x0)
        {
            this->clear();
            GmshReader gmshReader(meshFileName,A,x0);
            this->nodes()=gmshReader.nodes();
            //assert(nodeReader.size()==this->nodes().size());
            
            for(const auto& eIter : gmshReader.elements())
            {
                if(   eIter.second.type==4  // 4-node (linear) tetrahedron
                   || eIter.second.type==11 //11-node (quadratic) tetrahedron
                   )
                {
//                    typename SimplexTraits<dim,dim>::SimplexIDType key;
                    std::set<size_t> set;
                    for(int d=0;d<dim+1;++d)
                    {
                        set.insert(eIter.second.nodeIDs[d]);
//                        key[d]=eIter.second.nodeIDs[d]; // first dim+1 nodes are the principa vertex nodes. See Node order in http://gmsh.info/doc/texinfo/gmsh.html#Node-ordering
                    }
                    size_t regionID(eIter.second.tags[1]); // FINISH HERE, WHAT IS THE RIGHT TAG FOR A POLYCRYSTAL
                    this->elements().emplace(eIter.first,std::make_pair(typename SimplexTraits<dim,dim>::SimplexIDType(set),regionID));
                    
                }
            }
        }
        
    };
    
}	// close namespace
#endif
