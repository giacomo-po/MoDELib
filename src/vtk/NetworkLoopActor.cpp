/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_NetworkLoopActor_cpp_
#define model_NetworkLoopActor_cpp_

#include <iostream>
#include <deque>
#include <string>

#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkMath.h>
#include <vtkProperty.h>
#include <vtkTubeFilter.h>
#include <vtkPolyLine.h>
#include <vtkSphereSource.h>
#include <vtkArrowSource.h>
#include <vtkGlyph3D.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkLabeledDataMapper.h>
#include <vtkFloatArray.h>
#include <vtkTextProperty.h>
#include <vtkActor2D.h>
#include <vtkProperty2D.h>
#include <vtkRenderer.h>
#include <vtkLine.h>
#include <vtkCellData.h>
#include <vtkPolygon.h>
#include <vtkCellArray.h>

//#include <IDreader.h>
//#include <PlanarPolygon.h>
#include <DDconfigIO.h>
#include <MeshPlane.h>
#include <NetworkLoopActor.h>


namespace model
{
    
        
        /**********************************************************************/
        NetworkLoopActor::NetworkLoopActor(vtkGenericOpenGLRenderWindow* const renWin,vtkRenderer* const renderer) :
        /* init */ renderWindow(renWin)
        /* init */,mainLayout(new QGridLayout(this))
        /* init */,showLoops(new QCheckBox(this))
        /* init */,loopPolyData(vtkSmartPointer<vtkPolyData>::New())
        /* init */,loopMapper(vtkSmartPointer<vtkPolyDataMapper>::New())
        /* init */,loopActor(vtkSmartPointer<vtkActor>::New())
        {
            showLoops->setChecked(false);
            showLoops->setText("loops");

            mainLayout->addWidget(showLoops,0,0,1,1);
            this->setLayout(mainLayout);

            connect(showLoops,SIGNAL(stateChanged(int)), this, SLOT(modify()));

            loopMapper->SetInputData(loopPolyData);
            loopActor->SetMapper(loopMapper);
            loopActor->GetProperty()->SetColor(0.0, 0.0, 1.0); //(R,G,B)
            loopActor->SetVisibility(false);

            renderer->AddActor(loopActor);

        }
        

        
        /**********************************************************************/
        void NetworkLoopActor::updateConfiguration(const DDconfigIO<3>& configIO,vtkPolyData* const nodePolyData)
        {
            std::cout<<"Updating loops..."<<std::flush;
            const auto t0= std::chrono::system_clock::now();
            

            vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
            std::map<size_t,size_t> loopNodesMap; // nodeID,nodePositionInDDauxIO
            for(size_t k=0;k<configIO.loopNodes().size();++k)
            {
                const auto& node(configIO.loopNodes()[k]);
                loopNodesMap.emplace(node.sID,k);
                
                points->InsertNextPoint(node.P(0),
                                        node.P(1),
                                        node.P(2));
            }

            vtkSmartPointer<vtkCellArray> cells(vtkSmartPointer<vtkCellArray>::New());
            for(const auto& loopLink : configIO.loopLinks())
            {
                
                const auto sourceIter(loopNodesMap.find(loopLink.sourceID));
                if(sourceIter!=loopNodesMap.end())
                {
                    const auto sinkIter(loopNodesMap.find(loopLink.sinkID));
                    if(sinkIter!=loopNodesMap.end())
                    {
                        vtkSmartPointer<vtkLine> line(vtkSmartPointer<vtkLine>::New());
                        line->GetPointIds()->SetId(0, sourceIter->second); // the second 0 is the index of the Origin in linesPolyData's points
                        line->GetPointIds()->SetId(1, sinkIter->second);
                        cells->InsertNextCell(line);
                    }
                    else
                    {
                        throw std::runtime_error("Sink vertex not found in nodeMap");
                    }
                }
                else
                {
                    throw std::runtime_error("Source vertex not found in nodeMap");
                }
                
            }
            
//            // Create the polygon
//            vtkSmartPointer<vtkPolygon> polygon = vtkSmartPointer<vtkPolygon>::New();
//            polygon->GetPointIds()->SetNumberOfIds(n);
//            for (int j = 0; j < n; j++)
//            {
//                polygon->GetPointIds()->SetId(j, j);
//            }
//            // Add the polygon to a list of polygons
//            vtkSmartPointer<vtkCellArray> polygons = vtkSmartPointer<vtkCellArray>::New();
//            polygons->InsertNextCell(polygon);

            // Create a PolyData
//            vtkPolyData* polygonPolyData =  vtkPolyData::New();
            loopPolyData->SetPoints(points);
            loopPolyData->SetLines(cells);
            loopPolyData->Modified();

//            polygonPolyData->SetPolys(polygons);

            // create mapper and actor using this polydata  - the usual stuff
            
            

            std::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
        }
        
        /**********************************************************************/
        void NetworkLoopActor::modify()
        {
            
            loopActor->SetVisibility(showLoops->isChecked());


            renderWindow->Render();
        }

        
    
} // namespace model
#endif
