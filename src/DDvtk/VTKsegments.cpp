/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_VTKsegments_cpp_
#define model_VTKsegments_cpp_


#include <VTKsegments.h>

namespace model
{


    VTKsegments::VTKsegments(const std::string& folderName) :
    ///* init */ traitsIO(folderName)
    /* init */ perser(folderName+"/inputFiles/vtkSegments.txt")
    /* init */,material(perser.readString("materialFile",true),perser.readScalar<double>("absoluteTemperature",true))
    /* init */,quadPerLength(perser.readScalar<double>("quadPerLength",true))
    {
        
    }

//    typename VTKsegments::VectorDimI VTKsegments::getPbcFlags(const std::string& filename) const
//    {
//
//    }


    const std::vector<DislocationQuadraturePoint<3,0>>& VTKsegments::quadraturePoints() const
    {
        return *this;
    }

    std::vector<DislocationQuadraturePoint<3,0>>& VTKsegments::quadraturePoints()
    {
        return *this;
    }

    const std::vector<StressStraight<3>>& VTKsegments::segments() const
    {
        return *this;
    }

    std::vector<StressStraight<3>>& VTKsegments::segments()
    {
        return *this;
    }

const std::vector<typename VTKsegments::VectorDim>& VTKsegments::nodes() const
{
    return *this;
}

std::vector<typename VTKsegments::VectorDim>& VTKsegments::nodes()
{
    return *this;
}


    void VTKsegments::writeVTK(const std::string& vtkFilePrefix) const
    {
        const std::string quadFileName(vtkFilePrefix+"_quadrature.vtk");
        std::ofstream quadFile(quadFileName);
        quadFile<<"# vtk DataFile Version 3.0\n";
        quadFile<<"# Dislocation lines converted from MoDELib file\n";
        quadFile<<"ASCII\n";
        quadFile<<"DATASET UNSTRUCTURED_GRID\n";
        quadFile<<"POINTS "+std::to_string(quadraturePoints().size()+nodes().size())+" double\n";

        for(const auto& node : nodes())
        {// write all nodes positions
            quadFile<<std::setprecision(15)<<std::scientific<<node.transpose()*material.b_SI*1.0e10<<"\n";
        }
        
        for(const auto& qp : quadraturePoints())
        {// write all nodes positions
            quadFile<<std::setprecision(15)<<std::scientific<<qp.r.transpose()*material.b_SI*1.0e10<<"\n";
        }
        
        // Create map of quadraturePoints by segments
        std::map<std::pair<size_t,size_t>,std::set<size_t>> qPointMap;
        for(size_t q=0;q<quadraturePoints().size();++q)
        {
            const auto& qp(quadraturePoints()[q]);
            const std::pair<size_t,size_t> key(std::make_pair(qp.sourceID,qp.sinkID));
            qPointMap[key].insert(q);
        }
        
        if(qPointMap.size()!=segments().size())
        {
            std::cout<<"qPointMap.size()="<<qPointMap.size()<<std::endl;
            std::cout<<"segments().size()="<<segments().size()<<std::endl;
            throw std::runtime_error("qPointMap.size()!=segments().size()");
        }
        
        quadFile<<"\nCELLS "+std::to_string(qPointMap.size())+" "+std::to_string(quadraturePoints().size()+3*qPointMap.size())+"\n";
        for(const auto& pair : qPointMap)
        {
            quadFile<<pair.second.size()+2<<" "<<pair.first.first<<" ";
            for(const auto& q : pair.second)
            {
                quadFile<<q+nodes().size()<<" ";
            }
            quadFile<<pair.first.second<<"\n";
        }
        
    }


    void VTKsegments::updateQuadraturePoints(const std::string& vtkFilePrefix)
    {
        
        VectorDimI pbcFlags;
        const std::string caFileName(vtkFilePrefix+".ca");
        std::ifstream caFile(caFileName); //access vtk file
        
        
        
        if(caFile.is_open())
        {
            
            std::string line;
            while (std::getline(caFile, line)) //begin parsing vtk file for lines
            {
                if(line.find("SIMULATION_CELL_MATRIX")!=std::string::npos) //if POINTS is read, read npos
                {
                    for(int d=0;d<dim;++d)
                    {
                        std::getline(caFile, line);
                        std::stringstream ss(line);//store lines of vtk file in ss
                        ss>>cellMatrix(d,0)>>cellMatrix(d,1)>>cellMatrix(d,2);//store nodal positions to x,y,z
                    }
                    std::cout<<"SIMULATION_CELL_MATRIX=\n"<<cellMatrix<<std::endl;
                }
                
                if(line.find("PBC_FLAGS")!=std::string::npos) //if POINTS is read, read npos
                {
                    std::string temp;
                    std::stringstream ss(line);//store lines of vtk file in ss
                    ss>>temp;
                    ss>>pbcFlags(0)>>pbcFlags(1)>>pbcFlags(2);
                    std::cout<<"PBC_FLAGS="<<pbcFlags.transpose()<<std::endl;

                }
                
            }
            
            periodicShifts.clear();
            for(int i=-pbcFlags(0);i<pbcFlags(0)+1;++i)
            {
                for(int j=-pbcFlags(1);j<pbcFlags(1)+1;++j)
                {
                    for(int k=-pbcFlags(2);k<pbcFlags(2)+1;++k)
                    {
                        periodicShifts.push_back( i*cellMatrix.row(0)*1.0e-10/material.b_SI
                                                 +j*cellMatrix.row(1)*1.0e-10/material.b_SI
                                                 +k*cellMatrix.row(2)*1.0e-10/material.b_SI);
                    }
                }
            }
            std::cout<<"periodicShifts.size()="<<periodicShifts.size()<<std::endl;
        }
        else
        {
            throw std::runtime_error("Cannot open file "+ caFileName);
        }
        
        
        
        
    #ifdef _OPENMP
        const size_t nThreads = omp_get_max_threads();
    #else
        const size_t nThreads = 1;
    #endif
        std::cout <<"updating QuadraturePoints (" << nThreads << " threads) " << std::flush;
        const auto t0= std::chrono::system_clock::now();

    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
        for(size_t q=0;q<quadraturePoints().size();++q)
        {
            auto& qPoint(quadraturePoints()[q]);
            for(const auto& seg : segments())
            {
                for(const auto& shift : periodicShifts)
                {
                    qPoint.stress+=seg.stress(qPoint.r+shift);
                }
            }
        }
        std::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<defaultColor<<std::endl;
    
        writeVTK(vtkFilePrefix);
    }

    void VTKsegments::readVTK(const std::string& vtkFilePrefix)
    {
        nodes().clear();
        segments().clear();
        quadraturePoints().clear();
        
        const std::string vtkFileName(vtkFilePrefix+".vtk");
        std::ifstream vtkFile(vtkFileName); //access vtk file
        
        if(vtkFile.is_open())
        {
//            std::vector<VectorDim> rawPoints;
            std::vector<std::vector<size_t>> cells;
            std::vector<VectorDim> rawBurgers;
            const auto sfCoeffs(SplineSegmentBase<3,0>::sfCoeffs(0.0)); // shape function coefficients for linear segments
            
            int loopsSize(0);
            
            std::string line;
            
            
            while (std::getline(vtkFile, line)) //begin parsing vtk file for lines
            {
                if(line.find("POINTS")!=std::string::npos) //if POINTS is read, read npos
                {
                    const size_t firstSpace(line.find(' ')); //vtk file formatting const
                    const size_t secondSpace(line.find(' ',firstSpace+1)); //vtk file formatting const
                    
                    if(secondSpace>firstSpace) // if the line indentation is longer...
                    {
                        const size_t nodesSize(std::atoi(line.substr(firstSpace+1,secondSpace-firstSpace-1).c_str())); //read the number of nodes
                        
                        std::cout<<"Reading "<<nodesSize<<" nodes"<<std::endl; //print the number of nodes
                        
                        double x,y,z; //create doubles to store node positions
                        
                        for(size_t n=0;n<nodesSize;++n) //index over all nodes
                        {
                            std::getline(vtkFile, line); //access lines of vtk file
                            std::stringstream ss(line);//store lines of vtk file in ss
                            ss>>x>>y>>z;//store nodal positions to x,y,z
                            nodes().push_back((VectorDim()<<x,y,z).finished()*1.0e-10/material.b_SI);//construct nodes()  [<x_1,x_2,x_3>,loopNo]
                        }
                    }
                    else
                    {//POINTS was not read
                        throw std::runtime_error("Unable to extract number of nodes from line: "+line);
                    }
                }
                
                if(line.find("CELLS")!=std::string::npos) //if CELLS is read, read npos
                {
                    const size_t firstSpace(line.find(' ')); //vtk file formatting const
                    const size_t secondSpace(line.find(' ',firstSpace+1)); //vtk file formatting const
                    
                    if(secondSpace>firstSpace) // if the line indentation is longer...
                    {
                        
                        loopsSize=std::atoi(line.substr(firstSpace+1,secondSpace-firstSpace-1).c_str()); //read the total number of loops loopsSize
                        
                        std::cout<<"Reading "<<loopsSize<<" cells"<<std::endl; // print the total number of loops
                        //                    cells.emplace_back();
                        for(int n=0;n<loopsSize;++n) // index over the number of loops
                        {
                            std::getline(vtkFile, line); // read the contents of cells
                            std::stringstream ss(line); // store the contents of cells in ss
                            cells.emplace_back(); //populate cells
                            int temp; //creante int temp
                            int count =0; //create a count initialized at 0
                            
                            while (ss >> temp) // while the nodal index is less than temp...
                            {
                                if(count==0) //if count = 0
                                {
                                    count++; // increase count by 1
                                    continue; //rerun the loop
                                }
                                cells.back().push_back(temp); //populate cells wih node connections <node pts, Burgers vector>
                            }
                        }
                    }
                    else
                    {//CELLS was not read
                        throw std::runtime_error("Unable to extract number of loops from line: "+line);
                    }
                }
                
                
                if(line.find("burgers_vector_world")!=std::string::npos) //if burgers_vector_world is read, read npos
                {
                    double x,y,z; //initialize x,y,z, for burgers vector components
                    
                    for(int n=0;n<loopsSize;++n) //index over the number of loops
                    {
                        std::getline(vtkFile, line);//read the contents of burgers_vector_world
                        std::stringstream ss(line);// store the contents of burgers_vector_world in ss
                        ss>>x>>y>>z;// pass ss to x,y,z
                        rawBurgers.push_back((VectorDim()<<x,y,z).finished()*1.0e-10/material.b_SI);
                    }
                }
            }
            
            if(rawBurgers.size()==cells.size())
            {
                for(size_t k=0;k<cells.size();++k)
                {
                    const auto Burgers(rawBurgers[k]);
                    
                    for(size_t n=0;n<cells[k].size()-1;++n)
                    {
                        const size_t sourceID(cells[k][n]);
                        const size_t   sinkID(cells[k][n+1]);
                        const VectorDim& sourceP(nodes()[sourceID]);
                        const VectorDim& sinkP(nodes()[sinkID]);
                        
                        segments().emplace_back(material,sourceP,sinkP,Burgers);
                        const VectorDim chord(sinkP-sourceP);
                        const double chordLength(chord.norm());
                        const int qOrder=QuadratureDynamicType::lowerOrder(quadPerLength*chordLength);
                        const MatrixNcoeffDim dofs((MatrixNcoeffDim()<< sourceP.transpose(),sinkP.transpose()).finished());
                        
                        for(int q=0;q<qOrder;++q)
                        {
                            quadraturePoints().emplace_back(sourceID,sinkID,q,qOrder,sfCoeffs,dofs);
                        }
                    }
                    
                }
                std::cout<<"segments().size()="<<segments().size()<<std::endl;
                std::cout<<"quadraturePoints().size()="<<quadraturePoints().size()<<std::endl;
            }
            else
            {
                std::cout<<"rawBurgers.size()="<<rawBurgers.size()<<std::endl;
                std::cout<<"cells.size()="<<cells.size()<<std::endl;
                throw std::runtime_error("rawBurgers.size() NOT EQUAL TO cells.size()");
            }
            
        }
        else
        {
            std::cout<<"Cannot open "<<vtkFileName<<std::endl;
        }
        
        updateQuadraturePoints(vtkFilePrefix);
        
    }

}
#endif
