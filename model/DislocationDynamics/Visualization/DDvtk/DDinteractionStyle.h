/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DDinteractionStyle_H_
#define model_DDinteractionStyle_H_

#include <memory>
#include <stdlib.h>     //for using the function sleep


#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleTrackballCamera.h>
//#include <vtkInteractorStyleMultiTouchCamera.h>
#include <vtkWindowToImageFilter.h>
#include <vtkPNGWriter.h>
#include <vtkBMPWriter.h>
#include <vtkJPEGWriter.h>
#include <vtkPropPicker.h>
#include <vtkRendererCollection.h>

//#include <model/IO/VertexReader.h>
//#include <model/IO/EdgeReader.h>
#include <model/DislocationDynamics/Visualization/DDvtk/DislocationSegmentActor.h>
//#include <model/DislocationDynamics/Visualization/vtk/DislocationActors.h>
#include <model/DislocationDynamics/Visualization/DDvtk/PKActor.h>
#include <model/DislocationDynamics/Visualization/DDvtk/GlidePlaneActor.h>
#include <model/Utilities/TerminalColors.h>



namespace model
{
    
    //To update the display once you get new data, you would just update the
    //PolyData that is already attached to a mapper->actor->renderer (and
    //                                                                call Modified() on it if necessary) and the renderer would
    //automatically display the new points.
    
    class DDinteractionStyle :
    /* inherit */ public vtkInteractorStyleTrackballCamera
    //    public vtkInteractorStyleMultiTouchCamera
    {
        
        
        
        
        std::unique_ptr<DislocationSegmentActor> ddSegments;
        std::unique_ptr<PKActor> ddPK;
        std::unique_ptr<GlidePlaneActor> ddGP;
        
        
        
        int xCol;
        int yCol;
        
        double winFrac;
        
        std::string selectedKey;
        bool saveImage=false;
//        long int frameID;
        long int frameIncrement;
        long int currentFrameID;
        vtkActor    *LastPickedActor;
        vtkProperty *LastPickedProperty;
        
    public:

        static DDinteractionStyle* New();
        vtkTypeMacro(DDinteractionStyle, vtkInteractorStyleTrackballCamera);
        
        SimplicialMeshActor meshActor;
        vtkRenderer* ddRenderer;
        vtkRenderer* plotRenderer;

        
        /*************************************************************************/
        bool loadFrame(const long int& frameID)
        {
            
            bool frameLoaded=false;
            if(   currentFrameID!=frameID
               && (DislocationSegmentActor::VertexReaderType::isGood(frameID,false) || DislocationSegmentActor::VertexReaderType::isGood(frameID,true))
               && (DislocationSegmentActor::EdgeReaderType::isGood(frameID,false) || DislocationSegmentActor::EdgeReaderType::isGood(frameID,true))
               )
            {
                currentFrameID=frameID;
                std::cout<<greenBoldColor<<"Loading frame "<<currentFrameID<<defaultColor<<std::endl;
                
                vtkRenderWindowInteractor *rwi = this->Interactor;
                
                if (this->LastPickedActor)
                {// Before destroying all actors,restore properties of LastPickedActor
                    this->LastPickedActor->GetProperty()->DeepCopy(this->LastPickedProperty);
                    LastPickedActor = NULL; // LastPickedActor will be destroyed so we cannot further pick it
                }
                
                // Update ddActors
                meshActor.update(frameID);
                ddSegments.reset(new DislocationSegmentActor(frameID,ddRenderer));
                ddPK.reset(new PKActor(frameID,ddRenderer));
                ddGP.reset(new GlidePlaneActor(frameID,ddRenderer));
                
                // Update renderer
                rwi->Render();
                
                if (saveImage)
                {
                    vtkSmartPointer<vtkWindowToImageFilter> windowToImageFilter = vtkSmartPointer<vtkWindowToImageFilter>::New();
                    windowToImageFilter->SetInput(rwi->GetRenderWindow()/*renderWindow*/);
                    windowToImageFilter->SetMagnification(1); //set the resolution of the output image (3 times the current resolution of vtk render window)
                    windowToImageFilter->ReadFrontBufferOff(); // read from the back buffer
                    windowToImageFilter->Update();
                    
                    int imageType=0;
                    switch (imageType)
                    {
                        case 0:
                        {
                            const std::string fileName="png/image_"+std::to_string(frameID)+".png";
                            std::cout<<"saving image "<<fileName<<std::endl;
                            vtkSmartPointer<vtkPNGWriter> writer = vtkSmartPointer<vtkPNGWriter>::New();
                            writer->SetFileName(fileName.c_str());
                            writer->SetInputConnection(windowToImageFilter->GetOutputPort());
                            writer->Write();
                            break;
                            
                        }
                        case 1:
                        {
                            const std::string fileName="bmp/image_"+std::to_string(frameID)+".bmp";
                            std::cout<<"saving image "<<fileName<<std::endl;
                            vtkSmartPointer<vtkBMPWriter> writer = vtkSmartPointer<vtkBMPWriter>::New();
                            writer->SetFileName(fileName.c_str());
                            writer->SetInputConnection(windowToImageFilter->GetOutputPort());
                            writer->Write();
                            break;
                            
                        }
                        case 2:
                        {
                            const std::string fileName="jpg/image_"+std::to_string(frameID)+".jpg";
                            std::cout<<"saving image "<<fileName<<std::endl;
                            vtkSmartPointer<vtkJPEGWriter> writer = vtkSmartPointer<vtkJPEGWriter>::New();
                            writer->SetFileName(fileName.c_str());
                            writer->SetInputConnection(windowToImageFilter->GetOutputPort());
                            writer->Write();
                            break;
                            
                        }
                            
                        default:
                            break;
                    }
                    
                    rwi->Render();
                }
                frameLoaded=true;
            }
            else
            {
                std::cout<<"frame "<<frameID<<" not found. Reverting to "<<currentFrameID<<std::endl;
//                frameID=currentFrameID; // frameID is not a valid ID, return to last read
            }
            return frameLoaded;
        }
        
        
        
        
        /*************************************************************************/
        DDinteractionStyle() :
//        /* init list   */ ddRenderer(ddRen),
//        /* init list   */ plotRenderer(plotRen),
        /* init list   */ xCol(0),
        /* init list   */ yCol(0),
        /* init list   */ winFrac(0.5),
        /* init list   */ frameIncrement(1),
        /* init list   */ currentFrameID(-1)
        {
            LastPickedActor = NULL;
            LastPickedProperty = vtkProperty::New();
        }
        
        
        /*************************************************************************/
        virtual ~DDinteractionStyle()
        {
            LastPickedProperty->Delete();
        }
        
        /*************************************************************************/
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
        
        /*************************************************************************/
        virtual void OnChar()
        {/*! Overrides vtkInteractorStyleTrackballCamera::OnChar()
          * to avoid exiting the program on pressing "e"
          */
        }
        
        /*************************************************************************/
        virtual void OnKeyPress()
        {
            // Get the keypress
            vtkRenderWindowInteractor *rwi = this->Interactor;
            //        const std::string key = rwi->GetKeySym();
            std::string key = rwi->GetKeySym();
            
            if(key == "Escape")
            {
                std::cout << "Exiting DDvtk, goodbye!" << std::endl;
                exit(0);
            }
            
            // Handle an arrow key
            if(key == "Up")
            {
//                frameID-=frameIncrement;
                loadFrame(currentFrameID-frameIncrement);
            }
            
            if(key == "Down")
            {
//                frameID+=frameIncrement;
                loadFrame(currentFrameID+frameIncrement);
            }
            
//            if(key == "Right")
//            {
//                winFrac+=0.1;
//                if(winFrac>1.0)
//                {
//                    winFrac=1.0;
//                }
//                ddRenderer ->SetViewport(0.0,0,winFrac,1);
//                plotRenderer->SetViewport(winFrac,0,1,1);
//                
//                ddRenderer->Render();
//                plotRenderer->Render();
//                this->Interactor->Render();
//            }
            
//            if(key == "Left")
//            {                winFrac-=0.1;
//                if(winFrac<0.0)
//                {
//                    winFrac=0.0;
//                }
//                
//                ddRenderer->SetViewport(0.0,0,winFrac,1);
//                plotRenderer->SetViewport(winFrac,0,1,1);
//                
//                ddRenderer->Render();
//                plotRenderer->Render();
//                this->Interactor->Render();
//            }
            
            if(key == "Left")
            {
                if(loadFrame(currentFrameID-frameIncrement))
                {
//                    && !this->Interactor->GetShiftKey()
//                    std::cout<<this->Interactor->GetRepeatCount()<<std::endl;
                OnKeyPress();
//                    sleep(5000);         //make the programme waiting for 5 seconds

                }
//                while(loadFrame(currentFrameID-frameIncrement))
//                {
//                
//                }
            }
            
            if(key == "Right")
            {
                if(loadFrame(currentFrameID+frameIncrement))
                {
                OnKeyPress();
//                    sleep(5000);         //make the programme waiting for 5 seconds

                }

//                while(loadFrame(currentFrameID+frameIncrement))
//                {
//                    
//                }
            }
            
            if(key == "e")
            {
                selectedKey="e";
                std::cout<<"selecting objects: dislocation segments"<<std::endl;
                std::cout<<"    +/- to increase tube radius"<<std::endl;
                std::cout<<"      0 to show/hide zero-Burgers vector segments"<<std::endl;
                std::cout<<"      1 to show/hide node IDs"<<std::endl;
                std::cout<<"      2 to show/hide a specific node ID"<<std::endl;
            }
            

            
            if(key == "i")
            {
                std::cout<<"Enter frame increment (>0): "<<std::flush;
                long int temp;
                std::cin>>temp;
                if(temp>0)
                {
                    frameIncrement=temp;
                }
                else
                {
                    std::cout<<"frame increment must be >0. Reverting to frame increment="<<frameIncrement<<std::endl;
                }
            }
            
            
            if(key == "l")
            {
                std::cout<<"Enter frame# to load:"<<std::endl;
                long int frameID(0);
                std::cin>>frameID;
                loadFrame(frameID);
            }
            
            if(key == "m")
            {
                selectedKey="m";
                std::cout<<"selecting objects: mesh"<<std::endl;
                std::cout<<"    +/- to increase mesh displacement"<<std::endl;
                
            }
            
            if(key == "p")
            {
                if(selectedKey=="p")
                {
                    selectedKey=" ";
                    if(ddPK.get()!=nullptr)
                    {
                        PKActor::showPK=false;
                        ddPK->modify();
                        this->Interactor->Render();
                    }
                }
                else
                {
                    selectedKey="p";
                    std::cout<<"selecting objects: pk forces"<<std::endl;
                    std::cout<<"    +/- to increase vector size"<<std::endl;
                    if(ddPK.get()!=nullptr)
                    {
                        PKActor::showPK=true;
                        ddPK->modify();
                        this->Interactor->Render();
                    }
                    
                }
            }
            
            
            if(key == "s")
            {
                saveImage=!saveImage;
                std::cout<<"Saving image="<<saveImage<<std::endl;
            }
            
            
            if(key == "x")
            {
                std::cout<<"Enter column of x-axis data:"<<std::endl;
                int temp;
                std::cin>>temp;
                if(temp>0)
                {
                    xCol=temp;
                }
            }
            
            if(key == "y")
            {
                std::cout<<"Enter column of y-axis data:"<<std::endl;
                int temp;
                std::cin>>temp;
                if(temp>0)
                {
                    yCol=temp;
                }
            }
            
            
            if(key == "g")
            {
                if(selectedKey=="g")
                {
                    selectedKey=" ";
                    if(ddGP.get()!=nullptr)
                    {
                        GlidePlaneActor::showGlidePlane=false;
                        ddGP->modify();
                        this->Interactor->Render();
                    }
                }
                else
                {
                    selectedKey="g";
                    std::cout<<"selecting objects: glide planes"<<std::endl;
                    std::cout<<"    +/- to increase opacity"<<std::endl;
                    if(ddGP.get()!=nullptr)
                    {
                        GlidePlaneActor::showGlidePlane=true;
                        ddGP->modify();
                        this->Interactor->Render();
                    }
                }
            }

            if(key == "v")
            {
                if(selectedKey=="v")
                {
                    selectedKey=" ";
                    if(ddSegments.get()!=nullptr)
                    {
                        DislocationSegmentActor::showNodes=false;
                        ddSegments->modify();
                        this->Interactor->Render();
                    }
                }
                else
                {
                    selectedKey="v";
                    std::cout<<"selecting objects: Dislocation Nodes"<<std::endl;
                    std::cout<<"    +/- to increase vector size"<<std::endl;
                    if(ddSegments.get()!=nullptr)
                    {
                        DislocationSegmentActor::showNodes=true;
                        ddSegments->modify();
                        this->Interactor->Render();
                    }
                }
            }
            
            if(key == "w")
            {
                if(selectedKey=="w")
                {
                    selectedKey=" ";
                    if(ddSegments.get()!=nullptr)
                    {
                        DislocationSegmentActor::showVelocities=false;
                        ddSegments->modify();
                        this->Interactor->Render();
                    }
                }
                else
                {
                    selectedKey="w";
                    std::cout<<"selecting objects: nodal velocities"<<std::endl;
                    std::cout<<"    +/- to increase vector size"<<std::endl;
                    if(ddSegments.get()!=nullptr)
                    {
                        DislocationSegmentActor::showVelocities=true;
                        ddSegments->modify();
                        this->Interactor->Render();
                    }
                    
                    
                    
                    
                }
            }
            
            if(selectedKey=="e")
            {
                
                if(key == "equal")
                {
                    DislocationSegmentActor::tubeRadius*=2.0;
                    ddSegments->modify();
                    std::cout<<"tube radius="<<DislocationSegmentActor::tubeRadius<<std::endl;
                    this->Interactor->Render();
                }
                if(key == "minus")
                {
                    DislocationSegmentActor::tubeRadius*=0.5;
                    ddSegments->modify();
                    std::cout<<"tube radius="<<DislocationSegmentActor::tubeRadius<<std::endl;
                    this->Interactor->Render();
                }
                if(key == "0")
                {
                    DislocationSegmentActor::showZeroBuergers=!DislocationSegmentActor::showZeroBuergers;
                    ddSegments->modify();
                    this->Interactor->Render();
                }
                
                if(key == "1")
                {
                    DislocationSegmentActor::showNodeIDs=!DislocationSegmentActor::showNodeIDs;
                    ddSegments->modify();
                    this->Interactor->Render();
                }
                
                if(key == "2")
                {
                    DislocationSegmentActor::showSingleNode=!DislocationSegmentActor::showSingleNode;
                    if(DislocationSegmentActor::showSingleNode)
                    {
                        std::cout << "Enter node ID "<<std::flush;
                        std::cin >> DislocationSegmentActor::singleNodeID;
                        std::cout <<std::endl;
                    }
                    ddSegments->modify();
                    this->Interactor->Render();
                }
                
            }
            
            if(selectedKey=="m")
            {
                if(key == "equal")
                {
                    SimplicialMeshActor::dispCorr*=2.0;
                    std::cout<<"displacement amplification="<<SimplicialMeshActor::dispCorr<<std::endl;
                    meshActor.modifyPts();
                    this->Interactor->Render();
                }
                if(key == "minus")
                {
                    SimplicialMeshActor::dispCorr*=0.5;
                    std::cout<<"displacement amplification="<<SimplicialMeshActor::dispCorr<<std::endl;
                    meshActor.modifyPts();
                    this->Interactor->Render();
                }
                
            }
            
            if(selectedKey=="p")
            {
                if(key == "equal" && ddPK.get()!=nullptr)
                {
                    PKActor::pkFactor*=2.0;
                    ddPK->modify();
                    std::cout<<"force scaling="<<PKActor::pkFactor<<std::endl;
                    this->Interactor->Render();
                }
                if(key == "minus" && ddPK.get()!=nullptr)
                {
                    PKActor::pkFactor*=0.5;
                    ddPK->modify();
                    std::cout<<"force scaling="<<PKActor::pkFactor<<std::endl;
                    this->Interactor->Render();
                }
                
            }
            
            if(selectedKey=="w")
            {
                if(key == "equal" && ddSegments.get()!=nullptr)
                {
                    DislocationSegmentActor::velocityFactor*=2.0;
                    ddSegments->modify();
                    std::cout<<"velocity scaling="<<DislocationSegmentActor::velocityFactor<<std::endl;
                    this->Interactor->Render();
                }
                if(key == "minus" && ddSegments.get()!=nullptr)
                {
                    DislocationSegmentActor::velocityFactor*=0.5;
                    ddSegments->modify();
                    std::cout<<"velocity scaling="<<DislocationSegmentActor::velocityFactor<<std::endl;
                    this->Interactor->Render();
                }
                
            }
            
            if(selectedKey=="g")
            {
                if(key == "equal" && ddGP.get()!=nullptr)
                {
                    GlidePlaneActor::opacity*=2.0;
                    ddGP->modify();
                    this->Interactor->Render();
                }
                if(key == "minus" && ddGP.get()!=nullptr)
                {
                    GlidePlaneActor::opacity*=0.5;
                    ddGP->modify();
                    this->Interactor->Render();
                }
                
            }
            
            
            
            // Forward events
            vtkInteractorStyleTrackballCamera::OnKeyPress();// use base class OnKeyPress()
            /*
             http://www.vtk.org/doc/nightly/html/classvtkInteractorStyle.htmlf
             Keypress j / Keypress t: toggle between joystick (position sensitive) and trackball (motion sensitive) styles. In joystick style, motion occurs continuously as long as a mouse button is pressed. In trackball style, motion occurs when the mouse button is pressed and the mouse pointer moves.
             Keypress c / Keypress a: toggle between camera and actor modes. In camera mode, mouse events affect the camera position and focal point. In actor mode, mouse events affect the actor that is under the mouse pointer.
             Button 1: rotate the camera around its focal point (if camera mode) or rotate the actor around its origin (if actor mode). The rotation is in the direction defined from the center of the renderer's viewport towards the mouse position. In joystick mode, the magnitude of the rotation is determined by the distance the mouse is from the center of the render window.
             Button 2: pan the camera (if camera mode) or translate the actor (if actor mode). In joystick mode, the direction of pan or translation is from the center of the viewport towards the mouse position. In trackball mode, the direction of motion is the direction the mouse moves. (Note: with 2-button mice, pan is defined as <Shift>-Button 1.)
             Button 3: zoom the camera (if camera mode) or scale the actor (if actor mode). Zoom in/increase scale if the mouse position is in the top half of the viewport; zoom out/decrease scale if the mouse position is in the bottom half. In joystick mode, the amount of zoom is controlled by the distance of the mouse pointer from the horizontal centerline of the window.
             Keypress 3: toggle the render window into and out of stereo mode. By default, red-blue stereo pairs are created. Some systems support Crystal Eyes LCD stereo glasses; you have to invoke SetStereoTypeToCrystalEyes() on the rendering window.
             Keypress e: exit the application.
             Keypress f: fly to the picked point
             Keypress p: perform a pick operation. The render window interactor has an internal instance of vtkCellPicker that it uses to pick.
             Keypress r: reset the camera view along the current view direction. Centers the actors and moves the camera so that all actors are visible.
             Keypress s: modify the representation of all actors so that they are surfaces.
             Keypress u: invoke the user-defined function. Typically, this keypress will bring up an interactor that you can type commands in. Typing u calls UserCallBack() on the vtkRenderWindowInteractor, which invokes a vtkCommand::UserEvent. In other words, to define a user-defined callback, just add an observer to the vtkCommand::UserEvent on the vtkRenderWindowInteractor object.
             Keypress w: modify the representation of all actors so that they are wireframe.
             */
        }
        
    };
    vtkStandardNewMacro(DDinteractionStyle);
    
    
} // namespace model
#endif







