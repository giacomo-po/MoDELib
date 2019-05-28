/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_GmshReader_H_
#define model_GmshReader_H_

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <map>
#include <assert.h>

#include <Eigen/Dense>


namespace model
{
    
    struct GmshElement
    {
    
        int type;
        std::vector<int> tags;
        std::vector<size_t> nodeIDs;
    
    };
    
    class GmshReader : private std::map<size_t,Eigen::Vector3d>
    /*              */,private std::map<size_t,GmshElement>
    {// based on documentation here http://www.manpagez.com/info/gmsh/gmsh-2.2.6/gmsh_63.php
        
        /**********************************************************************/
        static int numberOfNodes(const int& type)
        {
        
            switch (type)
            {
                case 2: // 3-node triangle
                    return 3;
                    break;
                    
                case 4: // 4-node tetrahedron
                    return 4;
                    break;
                    
                default:
                    assert(0 && "FINISH HERE");
                    return 0;
                    break;
            }
        }
        
    public:
        
        /**********************************************************************/
        GmshReader(const std::string& meshFileName)
        {
            read(meshFileName);
        }
        
        /**********************************************************************/
        const std::map<size_t,Eigen::Vector3d>& nodes() const
        {
            return *this;
        }
        
        /**********************************************************************/
        std::map<size_t,Eigen::Vector3d>& nodes()
        {
            return *this;
        }
        
        /**********************************************************************/
        const std::map<size_t,GmshElement>& elements() const
        {
            return *this;
        }
        
        /**********************************************************************/
        std::map<size_t,GmshElement>& elements()
        {
            return *this;
        }
        
        /**********************************************************************/
        void read(const std::string& meshFileName)
        {
            nodes().clear();
            elements().clear();
            
            std::ifstream meshFile(meshFileName);
            
            if(meshFile.is_open())
            {
                std::string line;
                while (std::getline(meshFile, line))
                {
                    
                    if(line=="$MeshFormat")
                    {
                        //                        std::cout<<line<<std::endl;
                        // FINISH HERE
                        
                    }
                    
                    if(line=="$PhysicalNames")
                    {
                        std::getline(meshFile, line);
                        const int numberPhysicalNames=std::atoi(line.c_str());
                        std::cout<<"numberPhysicalNames="<<numberPhysicalNames<<std::endl;
                        for(int k=0;k<numberPhysicalNames;++k)
                        {
                            std::getline(meshFile, line);
                            // FINISH HERE
                        }
                    }
                    
                    if(line=="$Nodes")
                    {
                        
                        size_t nodeID(0);
                        Eigen::Vector3d nodePos;
                        
                        std::getline(meshFile, line);
                        const int numberNodes=std::atoi(line.c_str());
                        std::cout<<"number of Nodes="<<numberNodes<<std::endl;
                        for(int k=0;k<numberNodes;++k)
                        {
                            std::getline(meshFile, line);
                            std::stringstream ss(line);
                            ss>>nodeID;
                            ss>>nodePos(0);
                            ss>>nodePos(1);
                            ss>>nodePos(2);
                            const bool success=nodes().emplace(nodeID,nodePos).second;
                            if(!success)
                            {
                                std::cout<<"Could not insert node "<<nodeID<<" "<<nodePos.transpose()<<". Exiting."<<std::endl;
                                exit(EXIT_FAILURE);
                            }
                        }
                    }
                    
                    if(line=="$Elements")
                    {
                        
                        size_t elementID(0);
                        int type(0);
                        int tags(0);
                        int tag(0);
                        int nodeID(0);

                        std::getline(meshFile, line);
                        const int numberElements=std::atoi(line.c_str());
                        std::cout<<"number of Elements="<<numberElements<<std::endl;
                        for(int k=0;k<numberElements;++k)
                        {
                            GmshElement element;

                            
                            std::getline(meshFile, line);
                            std::stringstream ss(line);
                            ss>>elementID;
                            ss>>type;
                            ss>>tags;
                            for(int n=0;n<tags;++n)
                            {
                                ss>>tag;
                                element.tags.push_back(tag);
                            }
                            for(int n=0;n<numberOfNodes(type);++n)
                            {
                                ss>>nodeID;
                                element.nodeIDs.push_back(nodeID);
                            }
                            
                            const bool success=elements().emplace(elementID,element).second;
                            if(!success)
                            {
                                std::cout<<"Could not insert element "<<elementID<<". Exiting."<<std::endl;
                                exit(EXIT_FAILURE);
                            }
                        }
                    }
                }
                
            }
            else
            {
                std::cout<<"File "<<meshFileName<<" cannot be opened. Exiting."<<std::endl;
                exit(EXIT_FAILURE);
            }
            
            
        }
        
    };
}
#endif

