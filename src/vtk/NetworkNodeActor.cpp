/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_NetworkNodeActor_cpp_
#define model_NetworkNodeActor_cpp_

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
//#include <IDreader.h>
//#include <PlanarPolygon.h>
#include <DDconfigIO.h>
#include <MeshPlane.h>
#include <NetworkNodeActor.h>
#include <vtkRenderer.h>

namespace model
{
    
        
        /**********************************************************************/
        NetworkNodeActor::NetworkNodeActor(vtkGenericOpenGLRenderWindow* const renWin,vtkRenderer* const renderer) :
        /* init */ renderWindow(renWin)
        /* init */,mainLayout(new QGridLayout(this))
        /* init */,showNodes(new QCheckBox(this))
        /* init */,sliderNodeRadius(new QSlider(this))
        /* init */,showNodeLabels(new QCheckBox(this))
        /* init */,showSpecificNodeLabel(new QCheckBox(this))
        /* init */,showSpecificNodeLabelEdit(new QLineEdit("0"))
        /* init */,showVelocities(new QCheckBox(this))
        /* init */,velocityScaleEdit(new QLineEdit("1"))
        /* init */,nodePolyData(vtkSmartPointer<vtkPolyData>::New())
        /* init */,nodeGlyphs(vtkSmartPointer<vtkGlyph3D>::New())
        /* init */,nodeMapper(vtkSmartPointer<vtkPolyDataMapper>::New())
        /* init */,nodeActor(vtkSmartPointer<vtkActor>::New())
        /* init */,labelPolyData(vtkSmartPointer<vtkPolyData>::New())
        /* init */,labelMapper(vtkSmartPointer<vtkLabeledDataMapper>::New())
        /* init */,labelActor(vtkSmartPointer<vtkActor2D>::New())
        /* init */,velocityPolyData(vtkSmartPointer<vtkPolyData>::New())
        /* init */,velocityGlyphs(vtkSmartPointer<vtkGlyph3D>::New())
        /* init */,velocityMapper(vtkSmartPointer<vtkPolyDataMapper>::New())
        /* init */,velocityActor(vtkSmartPointer<vtkActor>::New())
        /* init */,specificNodeLabelPolyData(vtkSmartPointer<vtkPolyData>::New())
        /* init */,specificNodeLabelMapper(vtkSmartPointer<vtkLabeledDataMapper>::New())
        /* init */,specificNodeLabelActor(vtkSmartPointer<vtkActor2D>::New())
//        /* init */,singleNodeID(0)
        /* init */,nodeClr{{100,100,100},{0,255,255},{255,0,255},{1,1,1}}
        {
            
            const int initialRadius(12);
            
            showNodes->setText("nodes");
            showNodes->setChecked(false);
            nodeActor->SetVisibility(false);
            sliderNodeRadius->setEnabled(false);
            sliderNodeRadius->setMinimum(0);
            sliderNodeRadius->setMaximum(30);
            sliderNodeRadius->setValue(initialRadius);
            sliderNodeRadius->setOrientation(Qt::Horizontal);


            showSpecificNodeLabel->setText("selected nodes");
            showSpecificNodeLabel->setChecked(false);
            specificNodeLabelActor->SetVisibility(false);
            
            showNodeLabels->setText("node labels");
            showNodeLabels->setChecked(false);
            labelActor->SetVisibility(false);

            showVelocities->setText("velocities");
            showVelocities->setChecked(false);
            velocityActor->SetVisibility(false);
            velocityScaleEdit->setEnabled(false);


            mainLayout->addWidget(showNodes,0,0,1,1);
            mainLayout->addWidget(sliderNodeRadius,0,1,1,1);
            mainLayout->addWidget(showNodeLabels,1,0,1,1);
            mainLayout->addWidget(showSpecificNodeLabel,2,0,1,1);
            mainLayout->addWidget(showSpecificNodeLabelEdit,2,1,1,1);
            mainLayout->addWidget(showVelocities,3,0,1,1);
            mainLayout->addWidget(velocityScaleEdit,3,1,1,1);
            this->setLayout(mainLayout);

            connect(showNodes,SIGNAL(stateChanged(int)), this, SLOT(modify()));
            connect(showNodeLabels,SIGNAL(stateChanged(int)), this, SLOT(modify()));
            connect(showVelocities,SIGNAL(stateChanged(int)), this, SLOT(modify()));
            connect(velocityScaleEdit,SIGNAL(returnPressed()), this, SLOT(modify()));
            connect(showSpecificNodeLabel,SIGNAL(stateChanged(int)), this, SLOT(modify()));
            connect(showSpecificNodeLabelEdit,SIGNAL(returnPressed()), this, SLOT(modify()));
            connect(sliderNodeRadius,SIGNAL(valueChanged(int)), this, SLOT(modify()));

            
            nodeGlyphs->SetInputData(nodePolyData);
            nodeGlyphs->SetSourceConnection(vtkSmartPointer<vtkSphereSource>::New()->GetOutputPort());
            nodeGlyphs->ScalingOn();
            nodeGlyphs->SetScaleModeToScaleByVector();
            nodeGlyphs->SetScaleFactor(initialRadius);
            nodeGlyphs->SetColorModeToColorByScalar();
            nodeGlyphs->Update();
            nodeMapper->SetInputConnection(nodeGlyphs->GetOutputPort());
            nodeActor->SetMapper(nodeMapper);
            
            // Labels
            labelMapper->SetInputData(labelPolyData);
            labelMapper->SetLabelModeToLabelScalars();
            labelMapper->SetLabelFormat("%1.0f");
            labelMapper->GetLabelTextProperty()->SetFontSize(20);
            labelActor->SetMapper(labelMapper);
            labelActor->GetProperty()->SetColor(0.0, 0.0, 0.0); //(R,G,B)
            
            // Velocities
            velocityGlyphs->SetInputData(velocityPolyData);
            velocityGlyphs->SetSourceConnection(vtkSmartPointer<vtkArrowSource>::New()->GetOutputPort());
            velocityGlyphs->ScalingOn();
            velocityGlyphs->SetScaleModeToScaleByVector();
            velocityGlyphs->OrientOn();
            velocityGlyphs->ClampingOff();
            velocityGlyphs->SetVectorModeToUseVector();
            velocityGlyphs->SetIndexModeToOff();
            velocityMapper->SetInputConnection(velocityGlyphs->GetOutputPort());
            velocityMapper->ScalarVisibilityOff();
            velocityActor->SetMapper(velocityMapper);
            velocityActor->GetProperty()->SetColor(1.0, 0.0, 1.0); //(R,G,B)
            
            // Single node Label
           specificNodeLabelMapper->SetInputData(specificNodeLabelPolyData);
           specificNodeLabelMapper->SetLabelModeToLabelScalars();
           specificNodeLabelMapper->SetLabelFormat("%1.0f");
           specificNodeLabelMapper->GetLabelTextProperty()->SetFontSize(20);
           specificNodeLabelActor->SetMapper(specificNodeLabelMapper);
           specificNodeLabelActor->GetProperty()->SetColor(1.0, 0.0, 0.0); //(R,G,B)
//           specificNodeLabelActor->VisibilityOff();
            
            renderer->AddActor(nodeActor);
            renderer->AddActor(velocityActor);
            renderer->AddActor(labelActor);
            renderer->AddActor(specificNodeLabelActor);

        }
        

        
        /**********************************************************************/
        void NetworkNodeActor::updateConfiguration(const DDconfigIO<3>& configIO)
        {// https://stackoverflow.com/questions/6878263/remove-individual-points-from-vtkpoints
            std::cout<<"Updating nodes..."<<std::flush;
            const auto t0= std::chrono::system_clock::now();

            vtkSmartPointer<vtkPoints> nodePoints(vtkSmartPointer<vtkPoints>::New());
            vtkSmartPointer<vtkUnsignedCharArray> nodeColors(vtkSmartPointer<vtkUnsignedCharArray>::New());
            nodeColors->SetNumberOfComponents(3);

            vtkSmartPointer<vtkDoubleArray> nodeLabels(vtkSmartPointer<vtkDoubleArray>::New());
            nodeLabels->SetNumberOfComponents(1);
            nodeLabels->SetName("node IDs");
            
            vtkSmartPointer<vtkPoints> specificNodePoints(vtkSmartPointer<vtkPoints>::New());
            vtkSmartPointer<vtkDoubleArray> specificNodeLabels(vtkSmartPointer<vtkDoubleArray>::New());
            specificNodeLabels->SetNumberOfComponents(1);


            vtkSmartPointer<vtkDoubleArray> velocityVectors(vtkSmartPointer<vtkDoubleArray>::New());
            velocityVectors->SetNumberOfComponents(3);
            velocityVectors->SetName("nodeVelocity");
            
            std::set<size_t> singleNodes;
            if(showSpecificNodeLabel->isChecked())
            {
                try
                {
                    std::stringstream ss(showSpecificNodeLabelEdit->text() .toStdString());
                    size_t val;
                    while(ss>>val)
                    {
                        singleNodes.insert(val);
                    }
                }
                catch(std::exception& e)
                {
                    showSpecificNodeLabel->setText(QString::fromStdString(e.what()));
                }

            }
            
//            std::cout<<"Single nodes:"<<std::endl;
//            for(const auto& sn : singleNodes)
//            {
//                std::cout<<sn<<std::endl;
//            }

            for(const auto& node : configIO.nodes())
            {
                    if(singleNodes.find(node.sID)!=singleNodes.end())
                    {// current node is among selected
                        specificNodePoints->InsertNextPoint(node.P.data());
                        specificNodeLabels->InsertNextTuple1(node.sID);
                    }


                    
                nodePoints->InsertNextPoint(node.P.data());
                nodeLabels->InsertNextTuple1(node.sID);
                velocityVectors->InsertNextTuple(node.V.data()); // arrow vector
                nodeColors->InsertNextTypedTuple(node.meshLocation>2? this->nodeClr[3] : this->nodeClr[node.meshLocation]);

//                // Single node
//                if(node.sID==singleNodeID)
//                {
//                    singleNodePoint->InsertNextPoint(node.P.data());
//                    singleNodenodeLabels->InsertNextTuple1(node.sID);
//                }
            }

            nodePolyData->SetPoints(nodePoints);
            nodePolyData->GetPointData()->SetScalars(nodeColors);
            nodePolyData->Modified();

            labelPolyData->SetPoints(nodePoints);
            labelPolyData->GetPointData()->SetScalars(nodeLabels);
            labelPolyData->Modified();
            
            specificNodeLabelPolyData->SetPoints(specificNodePoints);
            specificNodeLabelPolyData->GetPointData()->SetScalars(specificNodeLabels);
            specificNodeLabelPolyData->Modified();
            
            velocityPolyData->SetPoints(nodePoints);
            velocityPolyData->GetPointData()->SetVectors(velocityVectors);
            velocityPolyData->Modified();
            std::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
        }
        
        /**********************************************************************/
        void NetworkNodeActor::modify()
        {
            
            nodeGlyphs->SetScaleFactor(sliderNodeRadius->value());
            sliderNodeRadius->setEnabled(showNodes->isChecked());
            nodeActor->SetVisibility(showNodes->isChecked());
            labelActor->SetVisibility(showNodeLabels->isChecked());
            specificNodeLabelActor->SetVisibility(showSpecificNodeLabel->isChecked());
            velocityActor->SetVisibility(showVelocities->isChecked());
            velocityScaleEdit->setEnabled(showVelocities->isChecked());
            const double vScaling(std::atof(velocityScaleEdit->text() .toStdString().c_str()));
            velocityGlyphs->SetScaleFactor(vScaling);

//                velocityActor->SetScale(vScaling);
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
//            if(this->showNodes)
//            {
//                nodeActor->VisibilityOn();
//            }
//            else
//            {
//                nodeActor->VisibilityOff();
//            }

            renderWindow->Render();
        }
        
    
} // namespace model
#endif
