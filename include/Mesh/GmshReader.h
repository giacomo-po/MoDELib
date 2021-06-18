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
        
        
        static std::map<int,std::string> getElementTypes();
        static int numberOfNodes(const int& type);
        
    public:

        const std::map<int,std::string> elementTypes;
        float version;
        int format;
        int size;
        int numberPhysicalNames;
        int numberNodes;
        int numberElements;
        std::map<int,size_t> numberElementsByType;

        GmshReader(const std::string& meshFileName);
        GmshReader(const std::string& meshFileName,const Eigen::Matrix<double,3,3>& A,const Eigen::Matrix<double,3,1>& x0);
        const std::map<size_t,Eigen::Vector3d>& nodes() const;
        std::map<size_t,Eigen::Vector3d>& nodes();
        const std::map<size_t,GmshElement>& elements() const;
        std::map<size_t,GmshElement>& elements();
        std::map<int,GmshPhysicalName>& physicalNames();
        const std::map<int,GmshPhysicalName>& physicalNames() const;
        void readVersion2(std::ifstream& meshFile,const Eigen::Matrix<double,3,3>& A,const Eigen::Matrix<double,3,1>& x0);
        void read(const std::string& meshFileName,const Eigen::Matrix<double,3,3>& A,const Eigen::Matrix<double,3,1>& x0);
        
    };
}
#endif

