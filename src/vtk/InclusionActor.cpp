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
#include <vtkPolyhedron.h>
#include <vtkIdList.h>
//#include <vtkXMLUnstructuredGridWriter.h>
#include <InclusionActor.h>


// VTK documentation
// http://vtk.1045678.n5.nabble.com/VTK-slow-to-display-300-vtkSphereSource-in-real-time-td5740730.html
// http://www.vtk.org/Wiki/VTK/Examples/Cxx/Visualization/Scaleglyph3D
// http://hiraku0n.blogspot.com/2012/10/writing-scalar-and-vector-field-in-vtk.html
// http://vtk.1045678.n5.nabble.com/Glyphing-vtkImageData-scalars-3D-as-arrows-td3199837.html

namespace model
{

    void InclusionActor::updateConfiguration()
    {
        std::cout<<"Updating inclusions..."<<std::flush;
        const auto t0= std::chrono::system_clock::now();
                
        vtkSmartPointer<vtkPoints> points(vtkSmartPointer<vtkPoints>::New());
        vtkSmartPointer<vtkFloatArray> diameters(vtkSmartPointer<vtkFloatArray>::New());
        vtkSmartPointer<vtkFloatArray> colors(vtkSmartPointer<vtkFloatArray>::New());
        
        diameters->SetName("diameters");
        colors->SetName("colors");
        
        for(const auto& inclusion : configFields.configIO.sphericalInclusions())
        {
            points->InsertNextPoint(inclusion.C(0), inclusion.C(1), inclusion.C(2));
            diameters->InsertNextValue(2.0*inclusion.a);  // origin of arrow
            colors->InsertNextValue(inclusion.phaseID);
        }
        
        grid->SetPoints(points);
        grid->GetPointData()->AddArray(diameters);
        grid->GetPointData()->SetActiveScalars("diameters"); // to set radius first
        grid->GetPointData()->AddArray(colors);
        
        vtkNew<vtkUnstructuredGrid> ugrid;
        vtkNew<vtkPoints> polyhedronPoints;
        vtkIdType polyhedronPointsIDs[configFields.configIO.polyhedronInclusionNodes().size()];
        long long k=0;
        for(const auto& node : configFields.configIO.polyhedronInclusionNodes())
        {
            polyhedronPoints->InsertNextPoint(node.P.data());
            polyhedronPointsIDs[k]=k;
            k++;
        }
        ugrid->SetPoints(polyhedronPoints);
        
        std::map<size_t,std::map<size_t,std::vector<size_t>>> faces;
        for(const auto& edge : configFields.configIO.polyhedronInclusionEdges())
        {
            const size_t& iID(edge.inclusionID);
            const size_t& fID(edge.faceID);
            const size_t& sourceID(edge.sourceID);
            faces[iID][fID].push_back(sourceID);
        }
        
        for(const auto& pair1 : faces)
        {
            vtkNew<vtkIdList> faces;
            for(const auto& pair2 : pair1.second)
            {
                faces->InsertNextId(pair2.second.size());
                for(const auto& nID : pair2.second)
                {
                    faces->InsertNextId(nID);
                }
            }
            ugrid->InsertNextCell(VTK_POLYHEDRON, configFields.configIO.polyhedronInclusionNodes().size(), polyhedronPointsIDs, pair1.second.size(), faces->GetPointer(0));
        }
        
        polyhedronMapper->SetInputData(ugrid);
        std::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
    }

    InclusionActor::InclusionActor(vtkGenericOpenGLRenderWindow* const renWin,vtkRenderer* const renderer,const DDconfigFields<3>& configFields_in) :
    /* init */ renderWindow(renWin)
    /* init */,mainLayout(new QGridLayout(this))
    /* init */,showInclusions(new QCheckBox(this))
    /* init */,sliderInclusionOpacity(new QSlider(this))
    /* init */,grid(vtkSmartPointer<vtkUnstructuredGrid>::New())
    /* init */,sphereSource(vtkSmartPointer<vtkSphereSource>::New())
    /* init */,glyph3D(vtkSmartPointer<vtkGlyph3D>::New())
    /* init */,mapper(vtkSmartPointer<vtkPolyDataMapper>::New())
    /* init */,actor(vtkSmartPointer<vtkActor>::New())
    /* init */,lookUpColors(vtkSmartPointer<vtkLookupTable>::New())
    /* init */,polyhedronMapper(vtkSmartPointer<vtkDataSetMapper>::New())
    /* init */,polyhedronActor(vtkSmartPointer<vtkActor>::New())
    /* init */,configFields(configFields_in)
    {
        showInclusions->setChecked(true);
        showInclusions->setText("show inclusions");
        mainLayout->addWidget(showInclusions,0,0,1,1);
        this->setLayout(mainLayout);
        connect(showInclusions,SIGNAL(stateChanged(int)), this, SLOT(modify()));
        
        sliderInclusionOpacity->setMinimum(0);
        sliderInclusionOpacity->setMaximum(10);
        sliderInclusionOpacity->setValue(5);
        sliderInclusionOpacity->setOrientation(Qt::Horizontal);
        connect(sliderInclusionOpacity,SIGNAL(valueChanged(int)), this, SLOT(modify()));
        mainLayout->addWidget(sliderInclusionOpacity,1,0,1,1);

        
        sphereSource->SetPhiResolution(20);
        sphereSource->SetThetaResolution(20);
        glyph3D->SetInputData(grid);
        glyph3D->SetSourceConnection(sphereSource->GetOutputPort());
        mapper->SetInputConnection(glyph3D->GetOutputPort());
        mapper->ScalarVisibilityOff();
        actor->GetProperty()->SetColor(0.5,0.5,0.5); // Give some color to the mesh. (1,1,1) is white
        actor->GetProperty()->SetOpacity(0.4); // Give some color to the mesh. (1,1,1) is white
        actor->SetMapper(mapper);
        
        polyhedronActor->SetMapper(polyhedronMapper);
        polyhedronActor->GetProperty()->SetColor(0.5,0.5,0.5); // Give some color to the mesh. (1,1,1) is white
        polyhedronActor->GetProperty()->SetOpacity(0.4); // Give some color to the mesh. (1,1,1) is white
                
        renderer->AddActor(actor);
        renderer->AddActor(polyhedronActor);
    }

    void InclusionActor::modify()
    {
        actor->SetVisibility(showInclusions->isChecked());
        actor->GetProperty()->SetOpacity(sliderInclusionOpacity->value()/10.0);

        renderWindow->Render();
    }

}
#endif
