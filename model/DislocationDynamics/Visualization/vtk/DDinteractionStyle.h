/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DDinteractionStyle_H_
#define model_DDinteractionStyle_H_

#include <vtkRenderer.h>

#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkWindowToImageFilter.h>
#include <vtkPNGWriter.h>
#include <vtkBMPWriter.h>
#include <vtkJPEGWriter.h>
#include <vtkPropPicker.h>


#include <model/Network/Readers/VertexReader.h>
#include <model/Network/Readers/EdgeReader.h>
#include <model/DislocationDynamics/Visualization/vtk/DislocationSegmentActor.h>
#include <model/DislocationDynamics/Visualization/vtk/DislocationActors.h>

long int frameID=0;
model::DislocationActors ddActors;
bool saveImage=false;


namespace model
{
    
    
    class DDinteractionStyle : public vtkInteractorStyleTrackballCamera
    {
        
        long int currentFrameID;
        
        
        void loadFrame()
        {
            
            if(currentFrameID!=frameID)
            {
                if(ddActors.isGood(frameID,false) || ddActors.isGood(frameID,true))
                {
                    currentFrameID=frameID; //
                    std::cout<<"loading frame "<<currentFrameID<<std::endl;
                    
                    vtkRenderWindowInteractor *rwi = this->Interactor;
                    
                    if (this->LastPickedActor)
                    {// Before destroying all actors,restore properties of LastPickedActor
                        this->LastPickedActor->GetProperty()->DeepCopy(this->LastPickedProperty);
                        LastPickedActor = NULL; // LastPickedActor will be destroyed so we cannot further pick it
                    }
                    
                    // Update ddActors
                    ddActors.update(frameID,this->CurrentRenderer);
                    
                    // Update renderer
                    rwi->Render();
                    
                    if (saveImage)
                    {
                        vtkSmartPointer<vtkWindowToImageFilter> windowToImageFilter = vtkSmartPointer<vtkWindowToImageFilter>::New();
                        windowToImageFilter->SetInput(rwi->GetRenderWindow()/*renderWindow*/);
                        windowToImageFilter->SetMagnification(1); //set the resolution of the output image (3 times the current resolution of vtk render window)
                        windowToImageFilter->SetInputBufferTypeToRGBA(); //also record the alpha (transparency) channel
                        windowToImageFilter->ReadFrontBufferOff(); // read from the back buffer
                        windowToImageFilter->Update();
                        
                        int imageType=0;
                        switch (imageType)
                        {
                            case 0:
                            {
                                const std::string fileName="png/image_"+std::to_string(frameID)+".png";
                                vtkSmartPointer<vtkPNGWriter> writer = vtkSmartPointer<vtkPNGWriter>::New();
                                writer->SetFileName(fileName.c_str());
                                writer->SetInputConnection(windowToImageFilter->GetOutputPort());
                                writer->Write();
                                rwi->Render();
                                break;
                                
                            }
                            case 1:
                            {
                                const std::string fileName="bmp/image_"+std::to_string(frameID)+".bmp";
                                vtkSmartPointer<vtkBMPWriter> writer = vtkSmartPointer<vtkBMPWriter>::New();
                                writer->SetFileName(fileName.c_str());
                                writer->SetInputConnection(windowToImageFilter->GetOutputPort());
                                writer->Write();
                                rwi->Render();
                                break;
                                
                            }
                            case 2:
                            {
                                const std::string fileName="jpg/image_"+std::to_string(frameID)+".jpg";
                                vtkSmartPointer<vtkJPEGWriter> writer = vtkSmartPointer<vtkJPEGWriter>::New();
                                writer->SetFileName(fileName.c_str());
                                writer->SetInputConnection(windowToImageFilter->GetOutputPort());
                                writer->Write();
                                rwi->Render();
                                break;
                                
                            }
                                
                            default:
                                break;
                        }
                    }
                    
                }
                else
                {
                    std::cout<<"frame "<<frameID<<" not found. Reverting to "<<currentFrameID<<std::endl;
                    frameID=currentFrameID; // frameID is not a valid ID, return to last read
                }
            }
            
        }
        
        
    public:
        static DDinteractionStyle* New();
        vtkTypeMacro(DDinteractionStyle, vtkInteractorStyleTrackballCamera);
        
        DDinteractionStyle() :
        /* init list   */ currentFrameID(-1)
        {
            LastPickedActor = NULL;
            LastPickedProperty = vtkProperty::New();
        }
        
        virtual ~DDinteractionStyle()
        {
            LastPickedProperty->Delete();
        }
        
        virtual void OnLeftButtonDown()
        {
            int* clickPos = this->GetInteractor()->GetEventPosition();
            
            // Pick from this location.
            vtkSmartPointer<vtkPropPicker>  picker =
            vtkSmartPointer<vtkPropPicker>::New();
            picker->Pick(clickPos[0], clickPos[1], 0, this->GetDefaultRenderer());
            
            // If we picked something before, reset its property
            if (this->LastPickedActor)
            {
                this->LastPickedActor->GetProperty()->DeepCopy(this->LastPickedProperty);
            }
            this->LastPickedActor = picker->GetActor();
            if (this->LastPickedActor)
            {
                // Save the property of the picked actor so that we can
                // restore it next time
                this->LastPickedProperty->DeepCopy(this->LastPickedActor->GetProperty());
                // Highlight the picked actor by changing its properties
                this->LastPickedActor->GetProperty()->SetColor(1.0, 0.0, 0.0);
                this->LastPickedActor->GetProperty()->SetDiffuse(1.0);
                this->LastPickedActor->GetProperty()->SetSpecular(0.0);
            }
            
            // Forward events
            vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
        }
        
        virtual void OnKeyPress()
        {
            // Get the keypress
            vtkRenderWindowInteractor *rwi = this->Interactor;
            //        const std::string key = rwi->GetKeySym();
            std::string key = rwi->GetKeySym();
            
            // Output the key that was pressed
            //       std::cout << "Pressed " << key << std::endl;
            
            // Handle an arrow key
            if(key == "Up")
            {
                frameID-=10;
                loadFrame();
            }
            
            if(key == "Down")
            {
                frameID+=10;
                loadFrame();
            }
            
            if(key == "l")
            {
                std::cout<<"Enter frame# to load:"<<std::endl;
                std::cin>>frameID;
                loadFrame();
            }
            
            if(key == "s")
            {
                saveImage=!saveImage;
                std::cout<<"Saving image="<<saveImage<<std::endl;
            }
            
            // Forward events
            vtkInteractorStyleTrackballCamera::OnKeyPress();// use base class OnKeyPress()
        }
        
        
    private:
        vtkActor    *LastPickedActor;
        vtkProperty *LastPickedProperty;
    };
    vtkStandardNewMacro(DDinteractionStyle);
    
    
} // namespace model
#endif







