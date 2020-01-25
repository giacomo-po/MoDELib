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
#include <chrono>

#include <Eigen/Dense>

#include <TerminalColors.h>


namespace model
{
    
    struct GmshPhysicalName
    {
        
        const int dim;
        const int tag;
        const std::string name;
        
    };

    struct GmshElement
    {
        int ID;
        int type;
        std::vector<int> tags;
        std::vector<size_t> nodeIDs;
    
    };
    
    class GmshReader : private std::map<size_t,Eigen::Vector3d>
    /*              */,private std::map<size_t,GmshElement>
    /*              */,private std::map<int,GmshPhysicalName>
    {// based on documentation here http://www.manpagez.com/info/gmsh/gmsh-2.2.6/gmsh_63.php
        
        
        static std::map<int,std::string> getElementTypes()
        {
            std::map<int,std::string> temp;
            temp.emplace(1,"2-node line");
            temp.emplace(2,"3-node triangle");
            temp.emplace(3,"4-node quadrangle");
            temp.emplace(4,"4-node tetrahedron");
            temp.emplace(5,"8-node hexahedron");
            temp.emplace(6,"6-node prism");
            temp.emplace(7,"5-node pyramid");
            temp.emplace(8,"3-node second-order line");
            temp.emplace(9,"6-node second-order triangle");
            temp.emplace(10,"9-node second-order quadrangle");
            temp.emplace(11,"10-node tetrahedron");
            temp.emplace(12,"27-node second order hexahedron");
            temp.emplace(13,"18-node second order prism");
            temp.emplace(14,"14-node second order pyramid");
            temp.emplace(15,"1-node point");
            return temp;
        }
        
        /**********************************************************************/
        static int numberOfNodes(const int& type)
        {
        
            switch (type)
            {
                case 1: // 2-node line
                    return 2;
                    break;

                case 2: // 3-node triangle
                    return 3;
                    break;
                    
                case 3: // 4-node quadrangle
                    return 4;
                    break;
                    
                case 4: // 4-node tetrahedron
                    return 4;
                    break;
                    
                case 8: // 3-node second order line (2 nodes associated with the vertices and 1 with the edge).
                    return 3;
                    break;
                    
                case 9: // 6-node second order triangle (3 nodes associated with the vertices and 3 with the edges).
                    return 6;
                    break;
                    
                case 11: // 10-node tetrahedron
                    return 10;
                    break;
                    
                case 12: // 27-node second order hexahedron (8 nodes associated with the vertices, 12 with the edges, 6 with the faces and 1 with the volume).
                    return 27;
                    break;
                    
                case 13: // 18-node second order prism (6 nodes associated with the vertices, 9 with the edges and 3 with the quadrangular faces).
                    return 18;
                    break;
                    
                case 14: // 14-node second order pyramid (5 nodes associated with the vertices, 8 with the edges and 1 with the quadrangular face).
                    return 14;
                    break;
                    
                case 15: // 1-node point.
                    return 1;
                    break;
                    
                default:
                {
                    std::cout<<"Element type="<<type<<std::endl;
                    assert(0 && "FINISH HERE");
                    return 0;
                    break;
                }
            }
        }
        
    public:

        const std::map<int,std::string> elementTypes;
        float version;
        int format;
        int size;
        int numberPhysicalNames;
        int numberNodes;
        int numberElements;
        std::map<int,size_t> numberElementsByType;
        
        /**********************************************************************/
        GmshReader(const std::string& meshFileName) :
        /* init */elementTypes(getElementTypes())
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
        
        std::map<int,GmshPhysicalName>& physicalNames()
        {
            return *this;
        }

        const std::map<int,GmshPhysicalName>& physicalNames() const
        {
            return *this;
        }

        
        void readVersion2(std::ifstream& meshFile)
        {
            std::string line;
            while (std::getline(meshFile, line))
            {
                
                if(line=="$PhysicalNames")
                {
                    std::getline(meshFile, line);
                    numberPhysicalNames=std::atoi(line.c_str());
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
                    numberNodes=std::atoi(line.c_str());
                    for(int k=0;k<numberNodes;++k)
                    {// Note Gmsh is only 3-dimensional
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
                    
                    std::getline(meshFile, line);
                    numberElements=std::atoi(line.c_str());
                    for(int k=0;k<numberElements;++k)
                    {
                        GmshElement element;
                        
                        
                        std::getline(meshFile, line);
                        std::stringstream ss(line);
                        ss>>element.ID;
                        ss>>element.type;
                        int tags(0);
                        ss>>tags;
                        for(int n=0;n<tags;++n)
                        {
                            element.tags.push_back(0);
                            ss>>element.tags.back();
                        }
                        for(int n=0;n<numberOfNodes(element.type);++n)
                        {
                            element.nodeIDs.push_back(0);
                            ss>>element.nodeIDs.back();
                        }
                        
                        const bool success=elements().emplace(element.ID,element).second;
                        if(!success)
                        {
                            std::cout<<"Could not insert element "<<element.ID<<". Exiting."<<std::endl;
                            exit(EXIT_FAILURE);
                        }

                        numberElementsByType[element.type]++;

                    }
                }
            }
        }
        
        /**********************************************************************/
        void read(const std::string& meshFileName)
        {
            nodes().clear();
            elements().clear();
            
            std::ifstream meshFile(meshFileName);
            
            if(meshFile.is_open())
            {
                const auto t0= std::chrono::system_clock::now();
                std::cout<<greenBoldColor<<"Reading msh file "<<meshFileName<<defaultColor<<std::flush;
                std::string line;
                while (std::getline(meshFile, line))
                {
                    
                    if(line=="$MeshFormat")
                    {
                        meshFile >> version >> format >> size;
                        

                        if(format!=0)
                        {
                            std::cout<<"Unsupported format. Only format=0 (ASCII) is supported. "<<std::endl;
                            exit (EXIT_FAILURE);
                        }
                        
                        if(version>=2.0 && version<3.0)
                        {
                            readVersion2(meshFile);
                        }
                        else
                        {
                                std::cout<<"Unsupported msh version"<<std::endl;
                                exit (EXIT_FAILURE);
                        }
                        break;
                    }
                }
                std::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
                std::cout<<"version="<<version<<std::endl;
                std::cout<<"format="<<format<<std::endl;
                std::cout<<"size="<<size<<std::endl;
                std::cout<<"numberPhysicalNames="<<numberPhysicalNames<<std::endl;
                std::cout<<"number of Nodes="<<numberNodes<<std::endl;
                std::cout<<"number of Elements="<<numberElements<<std::endl;
                for(const auto& pair : numberElementsByType)
                {
                    std::cout<<"    "<<pair.second<<" of type "<<pair.first<<": "<<elementTypes.at(pair.first)<<std::endl;
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

