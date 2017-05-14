// compile as cmake .
// make
// Examples:
// https://sites.google.com/site/siggraphasia2015democontest/vtk-tutorial

#include <string>

#include <vtkPolyDataMapper.h>
#include <vtkObjectFactory.h>
#include <vtkActor.h>
#include <vtkSmartPointer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
//#include <vtkPolyData.h>
#include <vtkSphereSource.h>

//#include <vtkProperty.h>
//#include <vtkPropPicker.h>
//#include <vtkTubeFilter.h>
//#include <vtkMath.h>
//#include <vtkObjectFactory.h> //vtkStandardNewMacro



#include <model/Mesh/SimplicialMesh.h>
#include <model/DislocationDynamics/Visualization/vtk/SimplicialMeshActor.h>
#include <model/DislocationDynamics/Visualization/vtk/DislocationSegmentActor.h>
//#include <model/DislocationDynamics/Visualization/vtk/DislocationActors.h>
#include <model/DislocationDynamics/Visualization/vtk/DDinteractionStyle.h>

#include <model/Utilities/EigenDataReader.h>







int main(int, char *[])
{
    
    // May be better to create: MeshRenderer, SegmentRenderer, and so on
    // so that one SegmentRenderer can be cleared of all actors while
    // MeshRenderer stays intact
    
    // Define the actors
    int meshID(0);
    model::EigenDataReader EDR;
    bool use_boundary=false;
    EDR.readScalarInFile("./DDinput.txt","use_boundary",use_boundary);
    if (use_boundary)
    {
        EDR.readScalarInFile("./DDinput.txt","meshID",meshID);
    }
//    model::SimplicialMeshActor meshActor();
//    meshActor.updated(meshID,renderer);
    
//    ddActors.read(frameID);
    
    // Define a renderer and add actors
    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    renderer->SetBackground(1,1,1); // Background color white

//    renderer->AddActor(meshActor.actor());
//    ddActors.clear();
//    ddActors.read(frameID);
//    ddActors.update(frameID,renderer);

//    for(auto& segmentActor : ddActors.segmentActors())
//    {
//        renderer->AddActor(segmentActor.tubeActor());
//    }
//    ddActors.addToRenderer(renderer);
    //renderer->AddActor(actor);

    // Define a render window and add the renderer
    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->AddRenderer(renderer);
    renderWindow->SetSize(1024,768); //(width, height)

    

//    renderWindow->Render();

    
    // Define the  window interactor
    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    renderWindowInteractor->SetRenderWindow(renderWindow);

    // Add interaction styles to the window
//    vtkSmartPointer<DDinteractionStyle> style = vtkSmartPointer<DDinteractionStyle>::New();
//    renderWindowInteractor->SetInteractorStyle(style);
//    style->SetCurrentRenderer(renderer);
    
    // Set the custom type to use for interaction.
    vtkSmartPointer<model::DDinteractionStyle> style =vtkSmartPointer<model::DDinteractionStyle>::New();
    renderWindowInteractor->SetInteractorStyle(style);
    style->SetDefaultRenderer(renderer);

    // Pupulate initial actors
    style->meshActor.update(meshID,renderer);
    style->ddActors.update(0,renderer);
    
//    model::SimplicialMeshActor meshActor;
//    meshActor.update(meshID,renderer);

    
    // Start
    renderWindow->LineSmoothingOn();
    renderWindow->PolygonSmoothingOn();
    renderWindow->PointSmoothingOn();
    renderWindow->SetMultiSamples(1);

    renderWindow->Render();
    renderWindowInteractor->Start();

    return EXIT_SUCCESS;
}

