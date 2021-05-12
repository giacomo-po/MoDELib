/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationNodeActor_h_
#define model_DislocationNodeActor_h_

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
#include <IDreader.h>
#include <PlanarPolygon.h>
#include <DDconfigIO.h>
#include <MeshPlane.h>
#include <DDconfigVtkBase.h>

namespace model
{
    struct DislocationNodeActor : public DDconfigVtkBase
    {
        
        vtkSmartPointer<vtkPolyData> nodePolyData;
        vtkSmartPointer<vtkGlyph3D> nodeGlyphs;
        vtkSmartPointer<vtkPolyDataMapper> nodeMapper;
        vtkSmartPointer<vtkActor> nodeActor;
        
        vtkSmartPointer<vtkPolyData> labelPolyData;
        vtkSmartPointer<vtkLabeledDataMapper> labelMapper;
        vtkSmartPointer<vtkActor2D> labelActor;

        vtkSmartPointer<vtkPolyData> velocityPolyData;
        vtkSmartPointer<vtkGlyph3D> velocityGlyphs;
        vtkSmartPointer<vtkPolyDataMapper> velocityMapper;
        vtkSmartPointer<vtkActor> velocityActor;
        
        vtkSmartPointer<vtkPolyData> singleNodeLabelPolyData;
        vtkSmartPointer<vtkLabeledDataMapper> singleNodeLabelMapper;
        vtkSmartPointer<vtkActor2D> singleNodeLabelActor;
        
        /**********************************************************************/
        DislocationNodeActor(vtkRenderer* const renderer) :
        /* init */ nodePolyData(vtkSmartPointer<vtkPolyData>::New()),
        /* init */ nodeGlyphs(vtkSmartPointer<vtkGlyph3D>::New()),
        /* init */ nodeMapper(vtkSmartPointer<vtkPolyDataMapper>::New()),
        /* init */ nodeActor(vtkSmartPointer<vtkActor>::New()),
        /* init */ labelPolyData(vtkSmartPointer<vtkPolyData>::New()),
        /* init */ labelMapper(vtkSmartPointer<vtkLabeledDataMapper>::New()),
        /* init */ labelActor(vtkSmartPointer<vtkActor2D>::New()),
        /* init */ velocityPolyData(vtkSmartPointer<vtkPolyData>::New()),
        /* init */ velocityGlyphs(vtkSmartPointer<vtkGlyph3D>::New()),
        /* init */ velocityMapper(vtkSmartPointer<vtkPolyDataMapper>::New()),
        /* init */ velocityActor(vtkSmartPointer<vtkActor>::New()),
        /* init */ singleNodeLabelPolyData(vtkSmartPointer<vtkPolyData>::New()),
        /* init */ singleNodeLabelMapper(vtkSmartPointer<vtkLabeledDataMapper>::New()),
        /* init */ singleNodeLabelActor(vtkSmartPointer<vtkActor2D>::New())
        {
            nodeGlyphs->SetInputData(nodePolyData);
            nodeGlyphs->SetSourceConnection(vtkSmartPointer<vtkSphereSource>::New()->GetOutputPort());
            nodeGlyphs->ScalingOn();
            nodeGlyphs->SetScaleModeToScaleByVector();
            nodeGlyphs->SetScaleFactor(2.0*tubeRadius*1.2);
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
            singleNodeLabelMapper->SetInputData(singleNodeLabelPolyData);
            singleNodeLabelMapper->SetLabelModeToLabelScalars();
            singleNodeLabelMapper->SetLabelFormat("%1.0f");
            singleNodeLabelMapper->GetLabelTextProperty()->SetFontSize(20);
            singleNodeLabelActor->SetMapper(singleNodeLabelMapper);
            singleNodeLabelActor->GetProperty()->SetColor(1.0, 0.0, 0.0); //(R,G,B)
            singleNodeLabelActor->VisibilityOff();
            
            renderer->AddActor(nodeActor);
            renderer->AddActor(velocityActor);
            renderer->AddActor(labelActor);
            renderer->AddActor(singleNodeLabelActor);

            modify();
        }
        

        
        /**********************************************************************/
        void updateConfiguration(const DDconfigIO<3>& configIO)
        {// https://stackoverflow.com/questions/6878263/remove-individual-points-from-vtkpoints
            std::cout<<"Updating nodes..."<<std::flush;
            const auto t0= std::chrono::system_clock::now();

            vtkSmartPointer<vtkPoints> nodePoints(vtkSmartPointer<vtkPoints>::New());
            vtkSmartPointer<vtkUnsignedCharArray> nodeColors(vtkSmartPointer<vtkUnsignedCharArray>::New());
            nodeColors->SetNumberOfComponents(3);

            vtkSmartPointer<vtkDoubleArray> nodeLabels(vtkSmartPointer<vtkDoubleArray>::New());
            nodeLabels->SetNumberOfComponents(1);
            nodeLabels->SetName("node IDs");
            
            vtkSmartPointer<vtkPoints> singleNodePoint(vtkSmartPointer<vtkPoints>::New());
            vtkSmartPointer<vtkDoubleArray> singleNodenodeLabels(vtkSmartPointer<vtkDoubleArray>::New());

            vtkSmartPointer<vtkDoubleArray> velocityVectors(vtkSmartPointer<vtkDoubleArray>::New());
            velocityVectors->SetNumberOfComponents(3);
            velocityVectors->SetName("nodeVelocity");

            for(const auto& node : configIO.nodeMap())
            {
                nodePoints->InsertNextPoint(node.second->P.data());
                nodeLabels->InsertNextTuple1(node.second->sID);
                velocityVectors->InsertNextTuple(node.second->V.data()); // arrow vactor
                nodeColors->InsertNextTypedTuple(node.second->meshLocation>2? this->nodeClr[3] : this->nodeClr[node.second->meshLocation]);

                // Single node
                if(node.second->sID==singleNodeID)
                {
                    singleNodePoint->InsertNextPoint(node.second->P.data());
                    singleNodenodeLabels->InsertNextTuple1(node.second->sID);
                }
            }

            nodePolyData->SetPoints(nodePoints);
            nodePolyData->GetPointData()->SetScalars(nodeColors);
            nodePolyData->Modified();

            labelPolyData->SetPoints(nodePoints);
            labelPolyData->GetPointData()->SetScalars(nodeLabels);
            labelPolyData->Modified();
            
            singleNodeLabelPolyData->SetPoints(singleNodePoint);
            singleNodeLabelPolyData->GetPointData()->SetScalars(singleNodenodeLabels);
            singleNodeLabelPolyData->Modified();
            
            velocityPolyData->SetPoints(nodePoints);
            velocityPolyData->GetPointData()->SetVectors(velocityVectors);
            velocityPolyData->Modified();
            model::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
        }
        
        /**********************************************************************/
        void modify()
        {
            
            nodeGlyphs->SetScaleFactor(2.0*this->tubeRadius*1.2);
            
            if(this->showVelocities)
            {
                velocityActor->VisibilityOn();
            }
            else
            {
                velocityActor->VisibilityOff();
            }
            
            if(this->showNodeIDs)
            {
                labelActor->VisibilityOn();
                
            }
            else
            {
                labelActor->VisibilityOff();
            }
            
            velocityGlyphs->SetScaleFactor(this->velocityFactor);
            
            
            if(this->showSingleNode)
            {
                // HERE WE SHOULD CHANGE THE NODE POSITION BASED ON NODE ID
                // OTHERWISE THE SELECTED NODE WILL BE VISIBLE ONLY UPON LOADING A NEW FRAME
                std::cout<<"RELOAD FRAME TO SHOW SELECTED NODE"<<std::endl;
                singleNodeLabelActor->VisibilityOn();
            }
            else
            {
                singleNodeLabelActor->VisibilityOff();
            }
            
            if(this->showNodes)
            {
                nodeActor->VisibilityOn();
            }
            else
            {
                nodeActor->VisibilityOff();
            }

            
        }
        
    };
    
} // namespace model
#endif
