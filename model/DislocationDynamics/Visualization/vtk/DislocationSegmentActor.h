/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationSegmentActor_H_
#define model_DislocationSegmentActor_H_

#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>

//#include <model/Mesh/SimplicialMesh.h>


namespace model
{
    struct DislocationSegmentActor
    {
        //    http://www.vtk.org/Wiki/VTK/Examples/Cxx/VisualizationAlgorithms/TubesWithVaryingRadiusAndColors
        
        vtkSmartPointer<vtkPoints> points;
        vtkSmartPointer<vtkCellArray> lines;
        vtkSmartPointer<vtkPolyData> polyData;
        vtkSmartPointer<vtkPolyDataMapper> lineMapper;
        vtkSmartPointer<vtkTubeFilter> tubeFilter;
        vtkSmartPointer<vtkPolyDataMapper> tubeMapper;
        
        
        DislocationSegmentActor() :
        /* init */ points(vtkSmartPointer<vtkPoints>::New()),
        /* init */ lines(vtkSmartPointer<vtkCellArray>::New()),
        /* init */ polyData(vtkSmartPointer<vtkPolyData>::New()),
        /* init */ lineMapper(vtkSmartPointer<vtkPolyDataMapper>::New()),
        /* init */ tubeFilter(vtkSmartPointer<vtkTubeFilter>::New()),
        /* init */ tubeMapper(vtkSmartPointer<vtkPolyDataMapper>::New())
        {
            double vX, vY, vZ;
            unsigned int nV = 256;      // No. of vertices
            unsigned int nCyc = 5;      // No. of spiral cycles
            double rT1 = 0.1, rT2 = 0.5;// Start/end tube radii
            double rS = 2000;              // Spiral radius
            double h = 10000;              // Height
            unsigned int nTv = 8;       // No. of surface elements for each tube vertex
            
            unsigned int i;
            
            // Create points and cells for the spiral
            for(i = 0; i < nV; i++)
            {
                // Spiral coordinates
                vX = rS * cos(2 * vtkMath::Pi() * nCyc * i / (nV - 1));
                vY = rS * sin(2 * vtkMath::Pi() * nCyc * i / (nV - 1));
                vZ = h * i / nV;
                points->InsertPoint(i, vX, vY, vZ);
            }
            
            lines->InsertNextCell(nV);
            for (i = 0; i < nV; i++)
            {
                lines->InsertCellPoint(i);
            }
            
            polyData->SetPoints(points);
            polyData->SetLines(lines);
            
            //        lineMapper->SetInputConnection(lines->GetOutputPort());
            lineMapper->SetInputData(polyData);
            
            tubeFilter->SetInputData(polyData);
            tubeFilter->SetRadius(500); //default is .5
            tubeFilter->SetNumberOfSides(50);
            tubeFilter->Update();
            
            tubeMapper->SetInputConnection(tubeFilter->GetOutputPort());
            tubeMapper->ScalarVisibilityOn();
            //        tubeMapper->SetScalarModeToUsePointFieldData();
            //        tubeMapper->SelectColorArray("Colors");
            
        }
        
        /**************************************************************************/
        vtkSmartPointer<vtkActor> lineActor() const
        {
            vtkSmartPointer<vtkActor> a =  vtkSmartPointer<vtkActor>::New();
            a->GetProperty()->SetColor(0.0,1.0,0.0); // Give some color to the tube
            a->SetMapper ( lineMapper );
            return a;
        }
        
        /**************************************************************************/
        vtkSmartPointer<vtkActor> tubeActor() const
        {
            vtkSmartPointer<vtkActor> a = vtkSmartPointer<vtkActor>::New();
            a->SetMapper(tubeMapper);
            a->GetProperty()->SetColor(1.0,0.0,0.0); // Give some color to the tube
            a->GetProperty()->SetOpacity(0.5); //Make the tube have some transparency.
            return a;
        }
        
    };
    
	
} // namespace model
#endif







