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

#include <DDconfigIO.h>
#include <MeshPlane.h>
#include <NetworkLoopActor.h>
#include <DislocationLoopPatches.h>

namespace model
{


/**********************************************************************/
NetworkLoopActor::NetworkLoopActor(vtkGenericOpenGLRenderWindow* const renWin,vtkRenderer* const renderer,
                                   const DDconfigFields<3>& configFields_in) :
/* init */ renderWindow(renWin)
/* init */,mainLayout(new QGridLayout(this))
/* init */,showLoops(new QCheckBox(this))
/* init */,slippedAreaBox(new QGroupBox(tr("&Slip Area")))
/* init */,sliderSlippedArea(new QSlider(this))
/* init */,loopPolyData(vtkSmartPointer<vtkPolyData>::New())
/* init */,loopMapper(vtkSmartPointer<vtkPolyDataMapper>::New())
/* init */,loopActor(vtkSmartPointer<vtkActor>::New())
/* init */,areaPolyData(vtkSmartPointer<vtkPolyData>::New())
/* init */,areaTriangleFilter(vtkSmartPointer<vtkTriangleFilter>::New())
/* init */,areaMapper(vtkSmartPointer<vtkPolyDataMapper>::New())
/* init */,areaActor(vtkSmartPointer<vtkActor>::New())
/* init */,configFields(configFields_in)
{
    showLoops->setChecked(false);
    showLoops->setText("Loops");
    
    slippedAreaBox->setCheckable(true);
    slippedAreaBox->setChecked(false);
    
    
    
    
    //            showSlippedArea->setChecked(false);
    //            showSlippedArea->setText("SlippedArea");
    //            sliderSlippedArea->setEnabled(false);
    sliderSlippedArea->setMinimum(0);
    sliderSlippedArea->setMaximum(10);
    sliderSlippedArea->setValue(5);
    sliderSlippedArea->setOrientation(Qt::Horizontal);
    
    //            colorSlippedArea->setWindowFlags(Qt::Widget );
    
    QVBoxLayout *slippedAreaLayout = new QVBoxLayout();
    slippedAreaLayout->addWidget(sliderSlippedArea);
    //            slippedAreaLayout->addWidget(colorSlippedArea);
    slippedAreaBox->setLayout(slippedAreaLayout);
    
    /* a few options that we must set for it to work nicely */
    //            colorSlippedArea->setOptions(
    //                            /* do not use native dialog */
    //                            QColorDialog::DontUseNativeDialog
    //                            /* you don't need to set it, but if you don't set this
    //                                the "OK" and "Cancel" buttons will show up, I don't
    //                                think you'd want that. */
    ////                            | QColorDialog::NoButtons
    //                );
    
    //            QGroupBox *groupBox = new QGroupBox(tr("&Push Buttons"));
    //                groupBox->setCheckable(true);
    //                groupBox->setChecked(true);
    
    
    mainLayout->addWidget(showLoops,0,0,1,1);
    mainLayout->addWidget(slippedAreaBox,1,0,1,1);
    //            mainLayout->addWidget(sliderSlippedArea,1,1,1,1);
    //            mainLayout->addWidget(colorSlippedArea,1,2,1,1);
    
    
    this->setLayout(mainLayout);
    
    connect(showLoops,SIGNAL(stateChanged(int)), this, SLOT(modify()));
    //            connect(slippedAreaBox,SIGNAL(stateChanged(int)), this, SLOT(modify()));
    connect(slippedAreaBox,SIGNAL(toggled(bool)), this, SLOT(modify()));
    
    connect(sliderSlippedArea,SIGNAL(valueChanged(int)), this, SLOT(modify()));
    
    loopMapper->SetInputData(loopPolyData);
    loopActor->SetMapper(loopMapper);
    loopActor->GetProperty()->SetColor(0.0, 0.0, 1.0); //(R,G,B)
    loopActor->SetVisibility(false);
    
    areaTriangleFilter->SetInputData(areaPolyData);
    //            filter->SetInputData(polygonPolyData);
    //            areaTriangleFilter->Update();
    
    //            areaMapper->SetInputData(areaPolyData);
    areaMapper->SetInputData(areaTriangleFilter->GetOutput());
    areaActor->SetMapper(areaMapper);
    areaActor->GetProperty()->SetColor(0.0, 1.0, 1.0); //(R,G,B)
    areaActor->SetVisibility(false);
    
    
    renderer->AddActor(loopActor);
    renderer->AddActor(areaActor);
    
    //            std::vector<Eigen::Vector2d> clippedPolygonTemp(SutherlandHodgman::clip(polygon2d,box2d));
    
}



void NetworkLoopActor::updateConfiguration()
{
    std::cout<<"Updating loops..."<<std::flush;
    const auto t0= std::chrono::system_clock::now();
    
    //            std::map<size_t,std::map<size_t,size_t>> loop2linkMap; // loopID->pair(source,sink)
    
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    std::map<size_t,size_t> loopNodesMap; // nodeID,nodePositionInDDauxIO
    for(size_t k=0;k<configFields.configIO.loopNodes().size();++k)
    {
        const auto& node(configFields.configIO.loopNodes()[k]);
        loopNodesMap.emplace(node.sID,k);
        
        points->InsertNextPoint(node.P(0),
                                node.P(1),
                                node.P(2));
    }
    
    vtkSmartPointer<vtkCellArray> cells(vtkSmartPointer<vtkCellArray>::New());
    for(const auto& loopLink : configFields.configIO.loopLinks())
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
                
                //                        loop2linkMap[loopLink.loopID].emplace(loopLink.sourceID,loopLink.sinkID);
                
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
    loopPolyData->SetPoints(points);
    loopPolyData->SetLines(cells);
    loopPolyData->Modified();
    
    // Slipped area
    vtkNew<vtkPoints> areaPoints;
    vtkNew<vtkCellArray> areaPolygons;
    
    size_t areaPointID(0);
    for(const auto& currentPatches : configFields.loopPatches())
    {
        for(const auto& pair : currentPatches.second.globalPatches())
        {
            vtkNew<vtkPolygon> polygon;
            for(const auto& globalPos : pair.second)
            {
                areaPoints->InsertNextPoint(globalPos(0),
                                            globalPos(1),
                                            globalPos(2));
                polygon->GetPointIds()->InsertNextId(areaPointID);
                areaPointID++;
            }
            areaPolygons->InsertNextCell(polygon);
        }
    }
    areaPolyData->SetPoints(areaPoints);
    areaPolyData->SetPolys(areaPolygons);
    areaTriangleFilter->Update();
    
    std::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
}

void NetworkLoopActor::modify()
{
    
    loopActor->SetVisibility(showLoops->isChecked());
    areaActor->SetVisibility(slippedAreaBox->isChecked());
    //            sliderSlippedArea->setEnabled(slippedAreaBox->isChecked());
    areaActor->GetProperty()->SetOpacity(sliderSlippedArea->value()/10.0);
    
    renderWindow->Render();
}

} // namespace model
#endif
