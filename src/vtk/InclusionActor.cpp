/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_InclusionActor_cpp_
#define model_InclusionActor_cpp_

#include <Eigen/Dense>

#include <QWidget>
#include <QGridLayout>
#include <QCheckBox>


#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkPolyData.h>
#include <vtkSphereSource.h>
#include <vtkGlyph3D.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkUnstructuredGrid.h>
#include <vtkLookupTable.h>

#include <InclusionActor.h>


// VTK documentation
// http://vtk.1045678.n5.nabble.com/VTK-slow-to-display-300-vtkSphereSource-in-real-time-td5740730.html
// http://www.vtk.org/Wiki/VTK/Examples/Cxx/Visualization/Scaleglyph3D
// http://hiraku0n.blogspot.com/2012/10/writing-scalar-and-vector-field-in-vtk.html
// http://vtk.1045678.n5.nabble.com/Glyphing-vtkImageData-scalars-3D-as-arrows-td3199837.html

namespace model
{
    
        void InclusionActor::updateConfiguration(const DDconfigIO<3>& configIO)
        {
            std::cout<<"Updating inclusions..."<<std::flush;
            const auto t0= std::chrono::system_clock::now();

            
            vtkSmartPointer<vtkPoints> points(vtkSmartPointer<vtkPoints>::New());
            vtkSmartPointer<vtkFloatArray> diameters(vtkSmartPointer<vtkFloatArray>::New());
            vtkSmartPointer<vtkFloatArray> colors(vtkSmartPointer<vtkFloatArray>::New());

            diameters->SetName("diameters");
            colors->SetName("colors");

            
            for(const auto& inclusion : configIO.eshelbyInclusions())
            {
                points->InsertNextPoint(inclusion.C(0), inclusion.C(1), inclusion.C(2));
                diameters->InsertNextValue(2.0*inclusion.a);  // origin of arrow
                colors->InsertNextValue(inclusion.phaseID);
            }

            grid->SetPoints(points);
            grid->GetPointData()->AddArray(diameters);
            grid->GetPointData()->SetActiveScalars("diameters"); // to set radius first
            grid->GetPointData()->AddArray(colors);

            std::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;

        }

        
        /**********************************************************************/
        InclusionActor::InclusionActor(vtkGenericOpenGLRenderWindow* const renWin,vtkRenderer* const renderer) :
        /* init */ renderWindow(renWin)
        /* init */,mainLayout(new QGridLayout(this))
        /* init */,showInclusions(new QCheckBox(this))
//        /* init */,points(vtkSmartPointer<vtkPoints>::New())
//        /* init */,diameters(vtkSmartPointer<vtkFloatArray>::New())
//        /* init */,colors(vtkSmartPointer<vtkFloatArray>::New())
        /* init */,grid(vtkSmartPointer<vtkUnstructuredGrid>::New())
        /* init */,sphereSource(vtkSmartPointer<vtkSphereSource>::New())
        /* init */,glyph3D(vtkSmartPointer<vtkGlyph3D>::New())
        /* init */,mapper(vtkSmartPointer<vtkPolyDataMapper>::New())
        /* init */,actor(vtkSmartPointer<vtkActor>::New())
        /* init */,lookUpColors(vtkSmartPointer<vtkLookupTable>::New())
        {
            showInclusions->setChecked(true);
            showInclusions->setText("show inclusions");
            mainLayout->addWidget(showInclusions,0,0,1,1);
            this->setLayout(mainLayout);
            connect(showInclusions,SIGNAL(stateChanged(int)), this, SLOT(modify()));

            
            
            
            
            sphereSource->SetPhiResolution(20);
            sphereSource->SetThetaResolution(20);
            
            // Set up glyph

            glyph3D->SetInputData(grid);
            glyph3D->SetSourceConnection(sphereSource->GetOutputPort());
            
            // Set up mapper
            mapper->SetInputConnection(glyph3D->GetOutputPort());
            //            mapper->ScalarVisibilityOff();
            
//            if(false)
//            {
//                mapper->SetScalarModeToUsePointFieldData(); // without, color are displayed regarding radius and not color label
//                mapper->SelectColorArray("colors"); // !!!to set color (nevertheless you will have nothing)
//                mapper->SetLookupTable(lookUpColors);
//
//            }
//            else
//            {
//                mapper->ScalarVisibilityOff();
//                actor->GetProperty()->SetColor(0.5,0.5,0.5); // Give some color to the mesh. (1,1,1) is white
//                actor->GetProperty()->SetOpacity(0.4); // Give some color to the mesh. (1,1,1) is white
//            }
            
            mapper->ScalarVisibilityOff();
            actor->GetProperty()->SetColor(0.5,0.5,0.5); // Give some color to the mesh. (1,1,1) is white
            actor->GetProperty()->SetOpacity(0.4); // Give some color to the mesh. (1,1,1) is white

            
            // Set up actor
            actor->SetMapper(mapper);
            
            
            // Add actor to renderer
            renderer->AddActor(actor);
            
            
        }
        
        
        /**********************************************************************/
        void InclusionActor::modify()
        {
            actor->SetVisibility(showInclusions->isChecked());

            renderWindow->Render();
        }
        
        
}
#endif
