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

    void VTKsegments::updateQuadraturePoints()
    {
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
                qPoint.stress+=seg.stress(qPoint.r);
            }
        }
        std::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<defaultColor<<std::endl;

    }

    void VTKsegments::readVTK(const std::string& vtkFileName)
    {
        segments().clear();
        quadraturePoints().clear();
        
        std::ifstream vtkFile(vtkFileName); //access vtk file
        
        if(vtkFile.is_open())
        {
            std::vector<VectorDim> rawPoints;
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
                            rawPoints.push_back((VectorDim()<<x,y,z).finished()*1.0e-10/material.b_SI);//construct rawPoints  [<x_1,x_2,x_3>,loopNo]
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
                        const VectorDim& sourceP(rawPoints[sourceID]);
                        const VectorDim& sinkP(rawPoints[sinkID]);
                        
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
        
        updateQuadraturePoints();
        
    }

}
#endif
