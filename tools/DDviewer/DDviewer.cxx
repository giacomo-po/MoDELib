#include <vtkPolyDataMapper.h>
#include <vtkObjectFactory.h>
#include <vtkActor.h>
#include <vtkSmartPointer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPolyData.h>
#include <vtkSphereSource.h>
#include <vtkInteractorStyleTrackballCamera.h>

#include <vtkProperty.h>
#include <vtkTubeFilter.h>
#include <vtkMath.h>
//#include <vtkObjectFactory.h> //vtkStandardNewMacro

#include <vtkWindowToImageFilter.h>
#include <vtkPNGWriter.h>

#include <model/Mesh/SimplicialMesh.h>
#include <model/DislocationDynamics/Visualization/vtk/SimplicialMeshActor.h>
#include <model/DislocationDynamics/Visualization/vtk/DislocationSegmentActor.h>
#include <model/DislocationDynamics/Visualization/vtk/DislocationRenderer.h>


//To update the display once you get new data, you would just update the
//PolyData that is already attached to a mapper->actor->renderer (and
//                                                                call Modified() on it if necessary) and the renderer would
//automatically display the new points.

int frameID=0;

// Define interaction style
class KeyPressInteractorStyle : public vtkInteractorStyleTrackballCamera
{
public:
    static KeyPressInteractorStyle* New();
    vtkTypeMacro(KeyPressInteractorStyle, vtkInteractorStyleTrackballCamera);
    
    virtual void OnKeyPress()
    {
        // Get the keypress
        vtkRenderWindowInteractor *rwi = this->Interactor;
        std::string key = rwi->GetKeySym();
        
        // Output the key that was pressed
        std::cout << "Pressed " << key << std::endl;
        
        // Handle an arrow key
        if(key == "Up")
        {
            for(int n=0;n<299;n++){
            frameID+=10;
            std::cout << "The up arrow was pressed." << std::endl;
            ////            this->CurrentRenderer->ResetCamera();
            this->CurrentRenderer->Clear();
            this->CurrentRenderer->RemoveAllViewProps();
            
//            model::SimplicialMesh<3> mesh;
//            mesh.readMesh(1);
//            model::SimplicialMeshActor meshActor(mesh);
//            this->CurrentRenderer->AddActor(meshActor.actor());
            
            model::DislocationRenderer DDrenderer;
            DDrenderer.read(frameID);
            DDrenderer.addActors(this->CurrentRenderer);
            
            rwi->Render();
            }
        }
        
        if(key == "Down")
        {
            frameID-=10;
            std::cout << "The up arrow was pressed." << std::endl;
            ////            this->CurrentRenderer->ResetCamera();
            this->CurrentRenderer->Clear();
            this->CurrentRenderer->RemoveAllViewProps();
            
//            model::SimplicialMesh<3> mesh;
//            mesh.readMesh(0);
//            model::SimplicialMeshActor meshActor(mesh);
//            this->CurrentRenderer->AddActor(meshActor.actor());
//
            
            model::DislocationRenderer DDrenderer;
            DDrenderer.read(frameID);
            DDrenderer.addActors(this->CurrentRenderer);

            rwi->Render();
        }
        
        // Handle a "normal" key
        if(key == "a")
        {
            std::cout << "The a key was pressed." << std::endl;

            rwi->Render();

            vtkSmartPointer<vtkWindowToImageFilter> windowToImageFilter = vtkSmartPointer<vtkWindowToImageFilter>::New();
            windowToImageFilter->SetInput(rwi->GetRenderWindow()/*renderWindow*/);
            windowToImageFilter->SetMagnification(1); //set the resolution of the output image (3 times the current resolution of vtk render window)
            windowToImageFilter->SetInputBufferTypeToRGBA(); //also record the alpha (transparency) channel
            windowToImageFilter->ReadFrontBufferOff(); // read from the back buffer
            windowToImageFilter->Update();
            
            vtkSmartPointer<vtkPNGWriter> writer = vtkSmartPointer<vtkPNGWriter>::New();
            writer->SetFileName("screenshot2.png");
            writer->SetInputConnection(windowToImageFilter->GetOutputPort());
            writer->Write();
            rwi->Render();

        }
        
        // Forward events
        vtkInteractorStyleTrackballCamera::OnKeyPress();
    }
    
};
vtkStandardNewMacro(KeyPressInteractorStyle);

int main(int, char *[])
{
    
    // May be better to create: MeshRenderer, SegmentRenderer, and so on
    // so that one SegmentRenderer can be cleared of all actors while
    // MeshRenderer stays intact
    
//    model::SimplicialMesh<3> mesh;
//    mesh.readMesh(0);
    
    model::SimplicialMeshActor meshActor(model::SimplicialMesh<3>(0));
    
//    model::DislocationSegmentActor sa;
    
    model::DislocationRenderer DDrenderer;
    
    DDrenderer.read(0);
    
    // A renderer and render window
    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->AddRenderer(renderer);
    
    // An interactor
    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    renderWindowInteractor->SetRenderWindow(renderWindow);
    
    vtkSmartPointer<KeyPressInteractorStyle> style = vtkSmartPointer<KeyPressInteractorStyle>::New();
    renderWindowInteractor->SetInteractorStyle(style);
    style->SetCurrentRenderer(renderer);
    
    //renderer->AddActor(actor);
    renderer->AddActor(meshActor.actor());
    DDrenderer.addActors(renderer);

//    renderer->AddActor(sa.lineActor());
//    renderer->AddActor(sa.tubeActor());
    
    renderer->SetBackground(1,1,1); // Background color white
    
    renderWindow->Render();
    
    renderWindowInteractor->Start();
    
    return EXIT_SUCCESS;
}

