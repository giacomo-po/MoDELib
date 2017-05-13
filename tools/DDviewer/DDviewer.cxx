// compile as cmake .
// make
// Examples:
// https://sites.google.com/site/siggraphasia2015democontest/vtk-tutorial
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
#include <vtkBMPWriter.h>
#include <vtkJPEGWriter.h>

#include <model/Mesh/SimplicialMesh.h>
#include <model/DislocationDynamics/Visualization/vtk/SimplicialMeshActor.h>
#include <model/DislocationDynamics/Visualization/vtk/DislocationSegmentActor.h>
#include <model/DislocationDynamics/Visualization/vtk/DislocationActors.h>
#include <model/Utilities/EigenDataReader.h>


//To update the display once you get new data, you would just update the
//PolyData that is already attached to a mapper->actor->renderer (and
//                                                                call Modified() on it if necessary) and the renderer would
//automatically display the new points.

long int frameID=0;
model::DislocationActors ddActors;

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
            //for(int n=0;n<299;n++){
            frameID+=10;
//            std::cout << "The up arrow was pressed." << std::endl;
            ////            this->CurrentRenderer->ResetCamera();
            //this->CurrentRenderer->Clear();
//            this->CurrentRenderer->RemoveAllViewProps();
            
//            model::SimplicialMesh<3> mesh;
//            mesh.readMesh(1);
//            model::SimplicialMeshActor meshActor(mesh);
//            this->CurrentRenderer->AddActor(meshActor.actor());
            
            //model::DislocationActors ddActors;
//            ddActors.clear();
//            ddActors.read(frameID);
            ddActors.update(frameID,this->CurrentRenderer);
//            for(auto& segmentActor : ddActors.segmentActors())
//            {
//                this->CurrentRenderer->AddActor(segmentActor.tubeActor());
//            }
//            ddActors.addToRenderer(this->CurrentRenderer);
            
            rwi->Render();
            //}
        }
        
        if(key == "Down")
        {
            frameID-=10;
//            std::cout << "The up arrow was pressed." << std::endl;
            ////            this->CurrentRenderer->ResetCamera();
            //this->CurrentRenderer->Clear();
 //           this->CurrentRenderer->RemoveAllViewProps();
            
//            model::SimplicialMesh<3> mesh;
//            mesh.readMesh(0);
//            model::SimplicialMeshActor meshActor(mesh);
//            this->CurrentRenderer->AddActor(meshActor.actor());
//
            
//            model::DislocationActors ddActors;
//            ddActors.clear();
//            ddActors.read(frameID);
            ddActors.update(frameID,this->CurrentRenderer);

//            for(auto& segmentActor : ddActors.segmentActors())
//            {
//                this->CurrentRenderer->AddActor(segmentActor.tubeActor());
//            }
//            ddActors.addToRenderer(this->CurrentRenderer);

            rwi->Render();
        }
        
        if(key == "Right")
        {
            //for(int n=0;n<299;n++){
            while (std::cin.get()=='n')
            {
            frameID+=10;
            //            std::cout << "The up arrow was pressed." << std::endl;
            ////            this->CurrentRenderer->ResetCamera();
            //this->CurrentRenderer->Clear();
            //            this->CurrentRenderer->RemoveAllViewProps();
            
            //            model::SimplicialMesh<3> mesh;
            //            mesh.readMesh(1);
            //            model::SimplicialMeshActor meshActor(mesh);
            //            this->CurrentRenderer->AddActor(meshActor.actor());
            
            //model::DislocationActors ddActors;
            //            ddActors.clear();
            //            ddActors.read(frameID);
            ddActors.update(frameID,this->CurrentRenderer);
            //            for(auto& segmentActor : ddActors.segmentActors())
            //            {
            //                this->CurrentRenderer->AddActor(segmentActor.tubeActor());
            //            }
            //            ddActors.addToRenderer(this->CurrentRenderer);
            
            rwi->Render();
            }
            //}
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
            
            int imageType=2;
            switch (imageType)
            {
                case 0:
                {
                    vtkSmartPointer<vtkPNGWriter> writer = vtkSmartPointer<vtkPNGWriter>::New();
                    writer->SetFileName("screenshot2.png");
                    writer->SetInputConnection(windowToImageFilter->GetOutputPort());
                    writer->Write();
                    rwi->Render();
                    break;

                }
                case 1:
                {
                    vtkSmartPointer<vtkBMPWriter> writer = vtkSmartPointer<vtkBMPWriter>::New();
                    writer->SetFileName("screenshot2.bmp");
                    writer->SetInputConnection(windowToImageFilter->GetOutputPort());
                    writer->Write();
                    rwi->Render();
                    break;

                }
                case 2:
                {
                    vtkSmartPointer<vtkJPEGWriter> writer = vtkSmartPointer<vtkJPEGWriter>::New();
                    writer->SetFileName("screenshot2.jpg");
                    writer->SetInputConnection(windowToImageFilter->GetOutputPort());
                    writer->Write();
                    rwi->Render();
                    break;

                }
                    
                default:
                    break;
            }


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
    
    // Define the actors
    int meshID(0);
    model::EigenDataReader EDR;
    bool use_boundary=false;
    EDR.readScalarInFile("./DDinput.txt","use_boundary",use_boundary);
    if (use_boundary)
    {
        EDR.readScalarInFile("./DDinput.txt","meshID",meshID);
    }
    model::SimplicialMeshActor meshActor(meshID);
    
//    ddActors.read(frameID);
    
    // Define a renderer and add actors
    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    renderer->SetBackground(1,1,1); // Background color white

    renderer->AddActor(meshActor.actor());
//    ddActors.clear();
//    ddActors.read(frameID);
    ddActors.update(frameID,renderer);

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

    
    renderWindow->LineSmoothingOn();
    renderWindow->PolygonSmoothingOn();
    renderWindow->PointSmoothingOn();
    renderWindow->SetMultiSamples(1);

    renderWindow->Render();

    
    // Define the  window interactor
    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    renderWindowInteractor->SetRenderWindow(renderWindow);

    // Add interaction styles to the window
    vtkSmartPointer<KeyPressInteractorStyle> style = vtkSmartPointer<KeyPressInteractorStyle>::New();
    renderWindowInteractor->SetInteractorStyle(style);
    style->SetCurrentRenderer(renderer);

    // Start
    renderWindowInteractor->Start();

    return EXIT_SUCCESS;
}

