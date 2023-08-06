/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_NetworkLinkActor_cpp_
#define model_NetworkLinkActor_cpp_

#include <numbers>
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

//#include <IDreader.h>
//#include <PlanarPolygon.h>
#include <DDconfigIO.h>
#include <MeshPlane.h>
#include <NetworkLinkActor.h>


namespace model
{
    
        
        /**********************************************************************/
        NetworkLinkActor::NetworkLinkActor(vtkGenericOpenGLRenderWindow* const renWin,vtkRenderer* const renderer) :
        /* init */ renderWindow(renWin)
        /* init */,mainLayout(new QGridLayout(this))
        /* init */,showLinks(new QCheckBox(this))
        /* init */,sliderLinksRadius(new QSlider(this))
        /* init */,showZeroLinks(new QCheckBox(this))
        /* init */,linksColorBox(new QComboBox(this))
//        /* init */,tubeRadius(1.0)
        /* init */,clr(ColorScheme::colorBurgers)
        /* init */,polyData(vtkSmartPointer<vtkPolyData>::New())
        /* init */,polyDataBnd(vtkSmartPointer<vtkPolyData>::New())
        /* init */,polyData0(vtkSmartPointer<vtkPolyData>::New())
        /* init */,tubeFilter(vtkSmartPointer<vtkTubeFilter>::New())
        /* init */,tubeFilterBnd(vtkSmartPointer<vtkTubeFilter>::New())
        /* init */,tubeFilter0(vtkSmartPointer<vtkTubeFilter>::New())
        /* init */,tubeMapper(vtkSmartPointer<vtkPolyDataMapper>::New())
        /* init */,tubeActor(vtkSmartPointer<vtkActor>::New())
        /* init */,tubeMapperBnd(vtkSmartPointer<vtkPolyDataMapper>::New())
        /* init */,tubeActorBnd(vtkSmartPointer<vtkActor>::New())
        /* init */,tubeMapper0(vtkSmartPointer<vtkPolyDataMapper>::New())
        /* init */,tubeActor0(vtkSmartPointer<vtkActor>::New())
        {
            showLinks->setChecked(true);
            showLinks->setText("links");
            sliderLinksRadius->setEnabled(true);
            sliderLinksRadius->setMinimum(0);
            sliderLinksRadius->setMaximum(30);
            sliderLinksRadius->setValue(5);
            sliderLinksRadius->setOrientation(Qt::Horizontal);

            
            showZeroLinks->setChecked(false);
            showZeroLinks->setText("0-Burgers links");
            tubeActor0->SetVisibility(false);

            linksColorBox->addItem("Burgers vector");
            linksColorBox->addItem("plane normal");
            linksColorBox->addItem("glissile/sessile");
            linksColorBox->addItem("edge/screw");
            
            mainLayout->addWidget(showLinks,0,0,1,1);
            mainLayout->addWidget(sliderLinksRadius,0,1,1,1);
            mainLayout->addWidget(linksColorBox,1,0,1,1);
            mainLayout->addWidget(showZeroLinks,2,0,1,1);

            this->setLayout(mainLayout);

            connect(showLinks,SIGNAL(stateChanged(int)), this, SLOT(modify()));
            connect(showZeroLinks,SIGNAL(stateChanged(int)), this, SLOT(modify()));
            connect(sliderLinksRadius,SIGNAL(valueChanged(int)), this, SLOT(modify()));


            tubeFilter->SetInputData(polyData);
            tubeFilter->SetRadius(5.0); // this must be a function similar to setColor
//            if(scaleRadiusByBurgers)
//            {
//                tubeFilter->SetVaryRadiusToVaryRadiusByAbsoluteScalar();
//            }
            tubeFilter->SetNumberOfSides(10);
            tubeFilter->Update();
            tubeMapper->SetInputConnection(tubeFilter->GetOutputPort());
            tubeMapper->ScalarVisibilityOn();
            tubeMapper->SetScalarModeToUseCellFieldData();
            tubeMapper->SelectColorArray("Colors");
            tubeActor->SetMapper(tubeMapper);
            
            // Boundary segments
            tubeFilterBnd->SetInputData(polyDataBnd);
            tubeFilterBnd->SetRadius(5.0); // this must be a function similar to setColor
            tubeFilterBnd->SetNumberOfSides(10);
            tubeFilterBnd->Update();
            tubeMapperBnd->SetInputConnection(tubeFilterBnd->GetOutputPort());
            tubeMapperBnd->ScalarVisibilityOn();
            tubeActorBnd->SetMapper(tubeMapperBnd);
            tubeActorBnd->GetProperty()->SetOpacity(0.3); //(R,G,B)
            
            // Zero-Burgers Segments
            tubeFilter0->SetInputData(polyData0);
            tubeFilter0->SetRadius(5.0); // this must be a function similar to setColor
            tubeFilter0->SetNumberOfSides(10);
            tubeFilter0->Update();
            tubeMapper0->SetInputConnection(tubeFilter0->GetOutputPort());
            tubeMapper0->ScalarVisibilityOn();
            tubeActor0->SetMapper(tubeMapper0);
            tubeActor0->GetProperty()->SetColor(0.5, 0.5, 0.5); //(R,G,B)
            tubeActor0->GetProperty()->SetOpacity(0.3); //(R,G,B)

            renderer->AddActor(tubeActor);
            renderer->AddActor(tubeActorBnd);
            renderer->AddActor(tubeActor0);


        }
        
const std::map<std::pair<size_t,size_t>,DislocationSegmentIO<3>>& NetworkLinkActor::segments() const
{
    return _segments;
}

        
        /**********************************************************************/
        void NetworkLinkActor::updateConfiguration(const DDconfigIO<3>& configIO,vtkPolyData* const nodePolyData)
        {// https://stackoverflow.com/questions/6878263/remove-individual-points-from-vtkpoints
            std::cout<<"Updating links..."<<std::flush;
            const auto t0= std::chrono::system_clock::now();
            
//            std::cout<<"Updating segments..."<<std::flush;
//            const auto t2= std::chrono::system_clock::now();
            _segments=configIO.segments(); // compute and store segments from configIO

            vtkSmartPointer<vtkFloatArray> radii(vtkSmartPointer<vtkFloatArray>::New());
            vtkSmartPointer<vtkUnsignedCharArray> colors(vtkSmartPointer<vtkUnsignedCharArray>::New());
            vtkSmartPointer<vtkUnsignedCharArray> colorsBnd(vtkSmartPointer<vtkUnsignedCharArray>::New());
            radii->SetName("TubeRadius");
            colors->SetName("Colors");
            colorsBnd->SetNumberOfComponents(3);
            colors->SetNumberOfComponents(3);
            vtkSmartPointer<vtkCellArray> cells(vtkSmartPointer<vtkCellArray>::New());
            vtkSmartPointer<vtkCellArray> cellsBnd(vtkSmartPointer<vtkCellArray>::New());
            vtkSmartPointer<vtkCellArray> cells0(vtkSmartPointer<vtkCellArray>::New());

            for (const auto& segment : segments())
            {
                
                auto itSource(configIO.nodeMap().find(segment.second.sourceID)); //source
                if(itSource!=configIO.nodeMap().end())
                {
                    auto   itSink(configIO.nodeMap().find(segment.second.sinkID)); //sink
                    if(itSink!=configIO.nodeMap().end())
                    {
                        vtkSmartPointer<vtkLine> line(vtkSmartPointer<vtkLine>::New());
                        line->GetPointIds()->SetId(0, std::distance(configIO.nodeMap().begin(),itSource)); // the second 0 is the index of the Origin in linesPolyData's points
                        line->GetPointIds()->SetId(1, std::distance(configIO.nodeMap().begin(),itSink));
                        
                        const auto chord(configIO.nodes()[itSink->second].P-configIO.nodes()[itSource->second].P);
                        const double burgersNorm(segment.second.b.norm());
                        if(burgersNorm>FLT_EPSILON)
                        {
                            Eigen::Matrix<int,3,1> colorVector=computeColor(segment.second.b,chord,segment.second.n);
                            unsigned char lineClr[3]={(unsigned char) colorVector(0),(unsigned char) colorVector(1),(unsigned char) colorVector(2)};
                            
                            if(segment.second.meshLocation!=0)
                            {// 0=MeshLocation::insideMesh
                                cellsBnd->InsertNextCell(line);
        //                        polyDataBnd->InsertNextCell(VTK_LINE,2,connectivity); //Connects the first and fourth point we inserted into a line
                                colorsBnd->InsertNextTypedTuple(lineClr);
                            }
                            else
                            {
        //                        polyData->InsertNextCell(VTK_LINE,2,connectivity); //Connects the first and fourth point we inserted into a line
                                cells->InsertNextCell(line);
                                if(segment.second.meshLocation==2 /*&& blackGrainBoundarySegments*/)
                                {
                                    unsigned char lineClr1[3]={1,1,1};
                                    colors->InsertNextTypedTuple(lineClr1);
                                }
                                else
                                {
                                    colors->InsertNextTypedTuple(lineClr);
                                }

                                radii->InsertNextValue(burgersNorm*sliderLinksRadius->value());
                            }
                        }
                        else
                        {
                            cells0->InsertNextCell(line);
        //                    polyData0->InsertNextCell(VTK_LINE,2,connectivity); //Connects the first and fourth point we inserted into a line
                        }
                    }
                    else
                    {
                        throw std::runtime_error("SINK VERTEX NOT FOUND IN nodeMap");
                    }
                }
                else
                {
                    throw std::runtime_error("SOURCE VERTEX NOT FOUND IN nodeMap");
                }
            }
            
            polyData->SetPoints(nodePolyData->GetPoints());
            polyData->SetLines(cells);
            polyData->GetCellData()->AddArray(colors);
            polyData->GetPointData()->AddArray(radii);
            polyData->GetPointData()->SetActiveScalars("TubeRadius");
            polyData->Modified();
            
            polyDataBnd->SetPoints(nodePolyData->GetPoints());
            polyDataBnd->SetLines(cellsBnd);
            polyDataBnd->GetCellData()->SetScalars(colorsBnd);
            polyDataBnd->Modified();
            
            polyData0->SetPoints(nodePolyData->GetPoints());
            polyData0->SetLines(cells0);
            polyData0->Modified();


            std::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
        }
        
        /**********************************************************************/
        void NetworkLinkActor::modify()
        {
            sliderLinksRadius->setEnabled(showLinks->isChecked());
            tubeActor->SetVisibility(showLinks->isChecked());
            linksColorBox->setEnabled(showLinks->isChecked());
            tubeActor0->SetVisibility(showZeroLinks->isChecked());
            tubeFilter->SetRadius(sliderLinksRadius->value()); // this must be a function similar to setColor
            tubeFilterBnd->SetRadius(sliderLinksRadius->value()); // this must be a function similar to setColor
            tubeFilter0->SetRadius(sliderLinksRadius->value()); // this must be a function similar to setColor


//            labelActor->SetVisibility(showLinks->isChecked());
//            velocityActor->SetVisibility(showLinks->isChecked());
            
//            nodeGlyphs->SetScaleFactor(2.0*this->tubeRadius*1.2);
//            
//            if(this->showVelocities)
//            {
//                velocityActor->VisibilityOn();
//            }
//            else
//            {
//                velocityActor->VisibilityOff();
//            }
//            
//            if(this->showNodeIDs)
//            {
//                labelActor->VisibilityOn();
//                
//            }
//            else
//            {
//                labelActor->VisibilityOff();
//            }
//            
//            velocityGlyphs->SetScaleFactor(this->velocityFactor);
//            
//            
//            if(this->showSingleNode)
//            {
//                // HERE WE SHOULD CHANGE THE NODE POSITION BASED ON NODE ID
//                // OTHERWISE THE SELECTED NODE WILL BE VISIBLE ONLY UPON LOADING A NEW FRAME
//                std::cout<<"RELOAD FRAME TO SHOW SELECTED NODE"<<std::endl;
//                singleNodeLabelActor->VisibilityOn();
//            }
//            else
//            {
//                singleNodeLabelActor->VisibilityOff();
//            }
//            
//            if(this->showLinks)
//            {
//                nodeActor->VisibilityOn();
//            }
//            else
//            {
//                nodeActor->VisibilityOff();
//            }

            renderWindow->Render();
        }

Eigen::Matrix<int,3,1> NetworkLinkActor::vector2Clr(VectorDim clrVector) const
{
    float clrTol=100.0*FLT_EPSILON;
    if(clrVector(0)<-clrTol)
    {// first component not zero but begative, flip color
        clrVector*=-1.0;
    }
    else if(fabs(clrVector(0))<=clrTol)
    {// first component is zero, use second component
        if(clrVector(1)<-clrTol)
        {// second component not zero but begative, flip color
            clrVector*=-1.0;
        }
        else if(fabs(clrVector(1))<=clrTol)
        {// second component is zero, use third component
            if(clrVector(2)<-clrTol)
            {
                clrVector*=-1.0;
            }
        }
    }
    
    clrVector = (clrVector + VectorDim::Ones(dim) * clrVector.norm()).eval();
    clrVector.normalize();
    return (clrVector*255).cast<int>();
}

/*********************************************************************/
Eigen::Matrix<int,3,1> NetworkLinkActor::computeColor(const VectorDim& burgers, const VectorDim& chord, const VectorDim& planeNormal) const
{
    
    VectorDim clrVector(VectorDim::Zero());
    
    
    
    switch (linksColorBox->currentIndex())
    {
            
        default:
        {
            clrVector = burgers.normalized();
            //                    flipColor(colorVector);
            break;
        }

        case 1:
        {
            clrVector=planeNormal;
            break;
        }
            
        case 2:
        {
            const bool isSessile= planeNormal.squaredNorm()<FLT_EPSILON || fabs(burgers.dot(planeNormal))>FLT_EPSILON;
            clrVector(0)= isSessile? 1.0 : 0.1;
            clrVector(1)= isSessile? 0.5 : 0.4;
            clrVector(2)= isSessile? 0.0 : 0.9;
            break;
        }
            
        case 3:
        {
            const float u = acos(std::fabs(chord.normalized().dot(burgers.normalized())))*2.0/std::numbers::pi;
            //                            RGBcolor rgb(RGBmap::getColor(std::fabs(tubeTangents.col(k).normalized().dot(burgers.normalized())),0,1));
            //                            colorVector << rgb.r, rgb.g, rgb.b;
            clrVector=(VectorDim()<<1.0,0.647,0.0).finished()*u+VectorDim::UnitZ()*(1.0-u);
            break;
        }
            //                    break;
            
    }

//    switch (clr)
//    {
//
//
//        case colorSessile:
//        {
//            const bool isSessile= planeNormal.squaredNorm()<FLT_EPSILON || fabs(burgers.dot(planeNormal))>FLT_EPSILON;
//            clrVector(0)= isSessile? 1.0 : 0.1;
//            clrVector(1)= isSessile? 0.5 : 0.4;
//            clrVector(2)= isSessile? 0.0 : 0.9;
//            break;
//        }
//
//        case colorNormal:
//        {
//            clrVector=planeNormal;
//            break;
//        }
//
//        case colorEdgeScrew:
//        {
//            const float u = acos(std::fabs(chord.normalized().dot(burgers.normalized())))*2.0/std::numbers::pi;
//            //                            RGBcolor rgb(RGBmap::getColor(std::fabs(tubeTangents.col(k).normalized().dot(burgers.normalized())),0,1));
//            //                            colorVector << rgb.r, rgb.g, rgb.b;
//            clrVector=(VectorDim()<<1.0,0.647,0.0).finished()*u+VectorDim::UnitZ()*(1-u);
//            break;
//        }
//            //                    break;
//
//        default:
//            clrVector = burgers.normalized();
//            //                    flipColor(colorVector);
//            break;
//    }
    
    return vector2Clr(clrVector);
    
//    float clrTol=100.0*FLT_EPSILON;
//    if(clrVector(0)<-clrTol)
//    {// first component not zero but begative, flip color
//        clrVector*=-1.0;
//    }
//    else if(fabs(clrVector(0))<=clrTol)
//    {// first component is zero, use second component
//        if(clrVector(1)<-clrTol)
//        {// second component not zero but begative, flip color
//            clrVector*=-1.0;
//        }
//        else if(fabs(clrVector(1))<=clrTol)
//        {// second component is zero, use third component
//            if(clrVector(2)<-clrTol)
//            {
//                clrVector*=-1.0;
//            }
//        }
//    }
//
//    clrVector = (clrVector + VectorDim::Ones(dim) * clrVector.norm()).eval();
//    clrVector.normalize();
//    return (clrVector*255).cast<int>();
}
        
    
} // namespace model
#endif
