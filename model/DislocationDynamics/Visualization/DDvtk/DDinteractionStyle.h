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
#include <fstream>

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
#include <vtkObjectFactory.h>

#include <model/DislocationDynamics/Visualization/DDvtk/DislocationSegmentActor.h>
#include <model/DislocationDynamics/Visualization/DDvtk/PKActor.h>
#include <model/DislocationDynamics/Visualization/DDvtk/GlidePlaneActor.h>
#include <model/DislocationDynamics/Visualization/DDvtk/InclusionActor.h>
#include <model/Utilities/TerminalColors.h>
//#include <model/IO/EigenDataReader.h>
#include <model/IO/TextFileParser.h>
#include <model/DislocationDynamics/Visualization/DDvtk/PlotActor.h>


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
        std::unique_ptr<InclusionActor> inclusions;
        std::unique_ptr<PlotActor> plot;
        vtkSmartPointer<vtkAxesActor> axes;
        vtkSmartPointer<vtkOrientationMarkerWidget> widget;

        
        int xCol;
        int yCol;
        double winFrac;
        std::string selectedKey;
        bool saveImage;
        int imageType;
        int imageMagnification;
        bool imageTransparentBackground;
        long int frameIncrement;
        long int currentFrameID;
        long int lastFrameID;
        vtkActor    *LastPickedActor;
        vtkProperty *LastPickedProperty;
        
        std::map<int,std::string> FlabelsMap;
        bool autoSpin;
        double degPerStep;
        Eigen::Matrix<double,3,1> spinAxis;
        
        bool axisWidgetEnabled;
    public:
        
        static DDinteractionStyle* New();
        vtkTypeMacro(DDinteractionStyle, vtkInteractorStyleTrackballCamera);
        
        SimplicialMeshActor meshActor;
        vtkRenderer* ddRenderer;
        vtkRenderer* plotRenderer;
        
        /**********************************************************************/
        DDinteractionStyle() :
        axes(vtkSmartPointer<vtkAxesActor>::New()),
        widget(vtkSmartPointer<vtkOrientationMarkerWidget>::New()),
        /* init list   */ xCol(0),
        /* init list   */ yCol(2),
        /* init list   */ winFrac(0.5),
        /* init list   */ saveImage(false),
        /* init list   */ imageType(1),
        /* init list   */ imageMagnification(1),
        /* init list   */ imageTransparentBackground(false),
        /* init list   */ frameIncrement(TextFileParser("./inputFiles/DD.txt").readScalar<int>("outputFrequency",false)),
        /* init list   */ currentFrameID(-1),
        /* init list   */ lastFrameID(currentFrameID),
        /* init list   */ autoSpin(false),
        /* init list   */ degPerStep(0.0),
        /* init list   */ spinAxis(Eigen::Matrix<double,3,1>::Zero()),
        /* init list   */ axisWidgetEnabled(true)
        {
            LastPickedActor = NULL;
            LastPickedProperty = vtkProperty::New();
            
//            model::EigenDataReader EDR;
//            EDR.readScalarInFile("./DDinput.txt","outputFrequency",frameIncrement);
            
            //                        PlotActor pa(plotRenderer);
            
            std::string filename("F/F_labels.txt");
            std::ifstream ifs ( filename.c_str() , std::ifstream::in );
            if (ifs.is_open())
            {
                std::string line;
                while (std::getline(ifs, line))
                {
                    std::stringstream ss(line);
                    
                    int colID;
                    ss>>colID;
                    std::string label;
                    std::string temp;
                    
                    while (ss >> temp)
                    {
                        label+=" "+temp;
                    }
                    
                    FlabelsMap.emplace(colID,label);
                    
                }
                
            }
            else
            {
                std::cout<<"CANNOT READ "<<filename<<std::endl;
            }
            
            
        }
        
        
        void init(vtkRenderer* _ddRenderer,vtkRenderer* _plotRenderer,const int& meshID)
        {
        
            ddRenderer=_ddRenderer;
            plotRenderer=_plotRenderer;
            
            meshActor.init(meshID,ddRenderer);
            
            loadFrame(0);
            
            plotRenderer->SetBackground(1,1,1);
            plotRenderer->SetViewport(0.5,0,1.0,1);
            
            widget->SetOutlineColor( 0.9300, 0.5700, 0.1300 );
            widget->SetOrientationMarker( axes );
            widget->SetInteractor( this->Interactor );
            widget->SetViewport( 0.0, 0.0, 0.4, 0.4 );
            widget->SetEnabled( axisWidgetEnabled );
            widget->InteractiveOn();

            
            plotRenderer->SetBackground(1,1,1);
            plotRenderer->SetViewport(0.5,0,1.0,1);
            
//            ddRenderer->ResetCamera();
//            renderWindow->Render();
        }
        
        
        /**********************************************************************/
        bool loadFrame(const long int& frameID)
        {
            
            bool frameLoaded=false;
            if(   currentFrameID!=frameID
               && (EVLio<3>::isBinGood(frameID) || EVLio<3>::isTxtGood(frameID))
//               (DislocationSegmentActor::VertexReaderType().isGood(frameID,false) || DislocationSegmentActor::VertexReaderType().isGood(frameID,true))
//               && (DislocationSegmentActor::EdgeReaderType().isGood(frameID,false) || DislocationSegmentActor::EdgeReaderType().isGood(frameID,true))
               )
            {
                lastFrameID=currentFrameID;
                currentFrameID=frameID;
                std::cout<<greenBoldColor<<"Loading frame "<<currentFrameID<<defaultColor<<std::endl;
                
                vtkRenderWindowInteractor *rwi = this->Interactor;
                
                if (this->LastPickedActor)
                {// Before destroying all actors,restore properties of LastPickedActor
                    this->LastPickedActor->GetProperty()->DeepCopy(this->LastPickedProperty);
                    LastPickedActor = NULL; // LastPickedActor will be destroyed so we cannot further pick it
                }
                
                // Update ddActors
                meshActor.update(frameID/*,lastFrameID,degPerStep,spinAxis*/);
                ddSegments.reset(new DislocationSegmentActor(frameID/*,lastFrameID,degPerStep,spinAxis*/,ddRenderer));
                ddPK.reset(new PKActor(frameID,ddRenderer));
                ddGP.reset(new GlidePlaneActor(frameID,ddRenderer));
                inclusions.reset(new InclusionActor(0,ddRenderer));
                plot.reset(new PlotActor(plotRenderer,xCol,yCol,currentFrameID,FlabelsMap));
                

                const double spinAxisNorm(spinAxis.norm());
                if(autoSpin && degPerStep && spinAxisNorm)
                {
                    spinAxis/=spinAxisNorm;
                    
                    Eigen::AngleAxisd R(degPerStep*M_PI/180.0,spinAxis);
                    double viewUp[3];
                    this->GetDefaultRenderer( )->GetActiveCamera( )->GetViewUp( viewUp[0], viewUp[1], viewUp[2] );
                    
                    double focalPoint[3];
                    this->GetDefaultRenderer( )->GetActiveCamera( )->GetFocalPoint( focalPoint[0], focalPoint[1], focalPoint[2] );
                    
                    double cameraPosition[3];
                    this->GetDefaultRenderer( )->GetActiveCamera( )->GetPosition( cameraPosition[0], cameraPosition[1], cameraPosition[2] );

                    Eigen::Map<Eigen::Vector3d> ViewUp(viewUp);
                    Eigen::Map<Eigen::Vector3d> FocalPoint(focalPoint);
                    Eigen::Map<Eigen::Vector3d> CameraPosition(cameraPosition);
                    
                    const Eigen::Vector3d newPosition=FocalPoint+R*(CameraPosition-FocalPoint);
                    const Eigen::Vector3d newViewUp=R*ViewUp;
                    
                    this->GetDefaultRenderer()->GetActiveCamera()->SetViewUp( newViewUp(0), newViewUp(1), newViewUp(2) );
//                    this->GetDefaultRenderer()->GetActiveCamera()->SetFocalPoint( focalPoint[0], focalPoint[1], focalPoint[2] );
                    this->GetDefaultRenderer()->GetActiveCamera()->SetPosition( newPosition[0], newPosition[1], newPosition[2] );

                }
                
                // Update renderer
                rwi->Render();
                
                if (saveImage)
                {
                    vtkSmartPointer<vtkWindowToImageFilter> windowToImageFilter = vtkSmartPointer<vtkWindowToImageFilter>::New();
                    windowToImageFilter->SetInput(rwi->GetRenderWindow()/*renderWindow*/);
                    windowToImageFilter->SetMagnification(imageMagnification); //set the resolution of the output image (3 times the current resolution of vtk render window)
                    if(imageTransparentBackground)
                    {
                        windowToImageFilter->SetInputBufferTypeToRGBA(); //also record the alpha (transparency) channel
                    }
//                    windowToImageFilter->ReadFrontBufferOff(); // read from the back buffer
                    windowToImageFilter->Update();
                    
                    //                    int imageType=0;
                    switch (imageType)
                    {
                        case 1:
                        {
                            const std::string fileName="png/image_"+std::to_string(frameID)+".png";
                            std::cout<<"saving image "<<fileName<<std::endl;
                            vtkSmartPointer<vtkPNGWriter> writer = vtkSmartPointer<vtkPNGWriter>::New();
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
                            
                        case 3:
                        {
                            const std::string fileName="bmp/image_"+std::to_string(frameID)+".bmp";
                            std::cout<<"saving image "<<fileName<<std::endl;
                            vtkSmartPointer<vtkBMPWriter> writer = vtkSmartPointer<vtkBMPWriter>::New();
                            writer->SetFileName(fileName.c_str());
                            writer->SetInputConnection(windowToImageFilter->GetOutputPort());
                            writer->Write();
                            break;
                            
                        }
                            
                            
                        default:
                            break;
                    }
                    
                    //                    rwi->ResetCamera();
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
        
        
        
        

        
        
        /**********************************************************************/
        virtual ~DDinteractionStyle()
        {
            LastPickedProperty->Delete();
        }
        
        /**********************************************************************/
        virtual void OnRightButtonDown()
        {
            
            double viewUp[3];
            this->GetDefaultRenderer( )->GetActiveCamera( )->GetViewUp( viewUp[0], viewUp[1], viewUp[2] );
            
            double focalPoint[3];
            this->GetDefaultRenderer( )->GetActiveCamera( )->GetFocalPoint( focalPoint[0], focalPoint[1], focalPoint[2] );
            
            double cameraPosition[3];
            this->GetDefaultRenderer( )->GetActiveCamera( )->GetPosition( cameraPosition[0], cameraPosition[1], cameraPosition[2] );
            
            double viewAngle= this->GetDefaultRenderer( )->GetActiveCamera( )->GetViewAngle(  );
            
            double parallelScale=this->GetDefaultRenderer( )->GetActiveCamera( )->GetParallelScale( );
            
            
            std::cout<<greenBoldColor<<"writing to cameraState.txt"<<defaultColor<<std::endl;
            std::cout<<"viewUp="<<viewUp[0]<<" "<<viewUp[1]<<" "<<viewUp[2]<<std::endl;
            std::cout<<"focalPoint="<<focalPoint[0]<<" "<<focalPoint[1]<<" "<<focalPoint[2]<<std::endl;
            std::cout<<"cameraPosition="<<cameraPosition[0]<<" "<<cameraPosition[1]<<" "<<cameraPosition[2]<<std::endl;
            std::cout<<"viewAngle="<<viewAngle<<std::endl;
            std::cout<<"parallelScale="<<parallelScale<<std::endl;
            
            std::ofstream myfile;
            myfile.open ("cameraState.txt");
            myfile<<viewUp[0]<<" "<<viewUp[1]<<" "<<viewUp[2]<<std::endl;
            myfile<<focalPoint[0]<<" "<<focalPoint[1]<<" "<<focalPoint[2]<<std::endl;
            myfile<<cameraPosition[0]<<" "<<cameraPosition[1]<<" "<<cameraPosition[2]<<std::endl;
            myfile<<viewAngle<<std::endl;
            myfile<<parallelScale<<std::endl;
            myfile.close();
            
            // Forward events
            vtkInteractorStyleTrackballCamera::OnRightButtonDown();
        }
        
        /*************************************************************************/
        virtual void OnLeftButtonDown()
        {
//            autoSpin=false;
            
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
            std::string key = rwi->GetKeySym();
            
            if(key == "Escape")
            {
                std::cout << "Exiting DDvtk, goodbye!" << std::endl;
                exit(0);
            }
            
            if(key == "c")
            {
                autoSpin=false;
                
                std::cout<<greenBoldColor<<"reading cameraState.txt"<<defaultColor<<std::endl;
                std::ifstream ifs ( "cameraState.txt" , std::ifstream::in );
                if (ifs.is_open())
                {
                    std::string line;
                    
                    std::getline(ifs, line);
                    double viewUp[3];
                    std::stringstream viewUpss(line);
                    viewUpss>>viewUp[0]>>viewUp[1]>>viewUp[2];
                    
                    std::getline(ifs, line);
                    double focalPoint[3];
                    std::stringstream focalPointss(line);
                    focalPointss>>focalPoint[0]>>focalPoint[1]>>focalPoint[2];
                    
                    std::getline(ifs, line);
                    double cameraPosition[3];
                    std::stringstream cameraPositionss(line);
                    cameraPositionss>>cameraPosition[0]>>cameraPosition[1]>>cameraPosition[2];
                    
                    std::getline(ifs, line);
                    double viewAngle;
                    std::stringstream viewAngless(line);
                    viewAngless>>viewAngle;
                    
                    std::getline(ifs, line);
                    double parallelScale;
                    std::stringstream parallelScaless(line);
                    parallelScaless>>parallelScale;
                    
                    
                    std::cout<<"viewUp="<<viewUp[0]<<" "<<viewUp[1]<<" "<<viewUp[2]<<std::endl;
                    std::cout<<"focalPoint="<<focalPoint[0]<<" "<<focalPoint[1]<<" "<<focalPoint[2]<<std::endl;
                    std::cout<<"cameraPosition="<<cameraPosition[0]<<" "<<cameraPosition[1]<<" "<<cameraPosition[2]<<std::endl;
                    std::cout<<"viewAngle="<<viewAngle<<std::endl;
                    std::cout<<"parallelScale="<<parallelScale<<std::endl;
                    
                    this->GetDefaultRenderer()->GetActiveCamera()->SetViewUp( viewUp[0], viewUp[1], viewUp[2] );
                    this->GetDefaultRenderer()->GetActiveCamera()->SetFocalPoint( focalPoint[0], focalPoint[1], focalPoint[2] );
                    this->GetDefaultRenderer()->GetActiveCamera()->SetPosition( cameraPosition[0], cameraPosition[1], cameraPosition[2] );
                    this->GetDefaultRenderer( )->GetActiveCamera( )->SetViewAngle( viewAngle );
                    this->GetDefaultRenderer( )->GetActiveCamera( )->SetParallelScale( parallelScale );
                    
                    this->Interactor->Render();
                    
                }
                else
                {
                    std::cout<<"file not found."<<std::endl;
                }
                ifs.close();
                
                
                
            }
            
            // Handle an arrow key
            if(key == "Up")
            {
                //                frameID-=frameIncrement;
                const bool success=loadFrame(currentFrameID-frameIncrement);
                if(this->Interactor->GetShiftKey() && success)
                {
                    //                    std::cout << "Shift held. ";
                    std::cout<<"Loading all frames"<<std::endl;
                    OnKeyPress();
                }
                
            }
            
            if(key == "Down")
            {
                //                frameID+=frameIncrement;
                const bool success=loadFrame(currentFrameID+frameIncrement);
                
                if(this->Interactor->GetShiftKey() && success)
                {
                    //                    std::cout << "Shift held. ";
                    std::cout<<"Loading all frames"<<std::endl;
                    OnKeyPress();
                }
            }
            
            if(key == "Right")
            {
                winFrac+=0.1;
                if(winFrac>1.0)
                {
                    winFrac=1.0;
                }
                ddRenderer ->SetViewport(0.0,0,winFrac,1);
                plotRenderer->SetViewport(winFrac,0,1,1);
                
                ddRenderer->Render();
                plotRenderer->Render();
                this->Interactor->Render();
            }
            
            if(key == "Left")
            {                winFrac-=0.1;
                if(winFrac<0.0)
                {
                    winFrac=0.0;
                }
                
                ddRenderer->SetViewport(0.0,0,winFrac,1);
                plotRenderer->SetViewport(winFrac,0,1,1);
                
                ddRenderer->Render();
                plotRenderer->Render();
                this->Interactor->Render();
            }
            


            if(key == "a")
            {
                selectedKey="a";
                std::cout<<"selecting objects: axis"<<std::endl;
                axisWidgetEnabled=!axisWidgetEnabled;
                widget->SetEnabled( axisWidgetEnabled );
                ddRenderer->Render();
                this->Interactor->Render();
            }
            
            if(key == "e")
            {
                selectedKey="e";
                std::cout<<"selecting objects: dislocation segments"<<std::endl;
                std::cout<<"    +/- to increase tube radius"<<std::endl;
                std::cout<<"      0 to show/hide zero-Burgers vector segments"<<std::endl;
                std::cout<<"      1 to show/hide boundary segments"<<std::endl;
                std::cout<<"      2 to turn on/off segment radii scaled by Burgers norm"<<std::endl;
                std::cout<<"      3 to color grain-boundary segments in black "<<std::endl;
                std::cout<<"      4 to color segments by Burgers vector "<<std::endl;
                std::cout<<"      5 to color segments by glissile/sessile "<<std::endl;

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
            {// std::cin return the current element of the input buffer, if the buffer is empty it waits for user input
                
                autoSpin=false;
                
                std::cout<<"Enter frame# to load:"<<std::endl;
                long int frameID;
                std::cin>>frameID;
                std::cout<<"entered "<<frameID<<std::endl;
                
//                cin.seekg( std::ios_base::seekdir::end);
//                long int frameID(0);
//                while(!(std::cin>>frameID))
//                {
////                    cin.seekg( std::ios_base::seekdir::end);
//                    std::cin.clear();
//                    std::cin.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
////                    std::cout <<"   invalid input, try again: "<<std::flush;
//                    //std::cin>>frameID;
//                    // code
//                }
//                std::cout<<std::endl;
//                std::cout<<"frameID="<<frameID<<std::endl;
//                std::cin>>frameID;
//                std::cin.clear();
//                std::cin.ignore(1000000,'\n');

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
            
            if(key == "r")
            {
                autoSpin=true;
                selectedKey="a";
                std::cout<<"Enter spin axis x y z: "<<std::flush;
                float x,y,z;
                std::cin>>x>>y>>z;
                spinAxis<<x,y,z;
                std::cout<<"Enter spin angle per step (deg/step): "<<std::flush;
                std::cin>>degPerStep;
            }
            
            if(key == "s")
            {
                
                selectedKey="s";
                //saveImage=true;
                std::cout<<"Save images menu:"<<std::endl;
                std::cout<<"    press space to start/stop saving images"<<std::endl;
                std::cout<<"    press 1 to save in png format"<<std::endl;
                std::cout<<"    press 2 to save in jgp format"<<std::endl;
                std::cout<<"    press 3 to save in bmp format"<<std::endl;
                std::cout<<"    press 0 to enable/disable transparent background"<<std::endl;
                std::cout<<"    press +/- to increase/decrease image resolution"<<std::endl;
                
            }
            
            
            if(key == "x")
            {
                std::cout<<"Enter column of x-axis data:"<<std::endl;
                for(const auto& pair : FlabelsMap)
                {
                    std::cout<<pair.first<<"    "<<pair.second<<"\n";
                }
                std::cin>>xCol;
                if(xCol<0)
                {
                    std::cout<<"wrong column number "<<xCol<<std::flush;
                    xCol=0;
                    std::cout<<". reverting to "<<xCol<<std::endl;
                }
                plot.reset(new PlotActor(plotRenderer,xCol,yCol,currentFrameID,FlabelsMap));
                this->Interactor->Render();
            }
            
            if(key == "y")
            {
                std::cout<<"Enter column of y-axis data:"<<std::endl;
                for(const auto& pair : FlabelsMap)
                {
                    std::cout<<pair.first<<"    "<<pair.second<<"\n";
                }
                std::cin>>yCol;
                if(yCol<0)
                {
                    std::cout<<"wrong column number "<<yCol<<std::flush;
                    yCol=1;
                    std::cout<<". reverting to "<<yCol<<std::endl;
                }
                plot.reset(new PlotActor(plotRenderer,xCol,yCol,currentFrameID,FlabelsMap));
                this->Interactor->Render();
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
            
            
            
            if(key == "t")
            {
                if(selectedKey=="t")
                {
                    selectedKey=" ";
                    if(ddSegments.get()!=nullptr)
                    {
                        DislocationSegmentActor::showSlippedArea=false;
                        ddSegments->modify();
                        this->Interactor->Render();
                    }
                }
                else
                {
                    selectedKey="t";
                    std::cout<<"selecting objects: Slipped Areas"<<std::endl;
                    std::cout<<"    +/- to increase/decrease opacity"<<std::endl;
                    
                    if(ddSegments.get()!=nullptr)
                    {
                        DislocationSegmentActor::showSlippedArea=true;
                        ddSegments->modify();
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
                    std::cout<<"      1 to show/hide node IDs"<<std::endl;
                    std::cout<<"      2 to show/hide a specific node ID"<<std::endl;
                    
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
                    DislocationSegmentActor::showBoundarySegments=!DislocationSegmentActor::showBoundarySegments;
                    ddSegments->modify();
                    this->Interactor->Render();
                }
                if(key == "2")
                {
                    DislocationSegmentActor::scaleRadiusByBurgers=!DislocationSegmentActor::scaleRadiusByBurgers;
                    std::cout<<"scaleRadiusByBurgers="<<DislocationSegmentActor::scaleRadiusByBurgers<<std::endl;
                    ddSegments->modify();
                    this->Interactor->Render();
                }
                if(key == "3")
                {
                    DislocationSegmentActor::blackGrainBoundarySegments=!DislocationSegmentActor::blackGrainBoundarySegments;
                    std::cout<<"blackGrainBoundarySegments="<<DislocationSegmentActor::blackGrainBoundarySegments<<std::endl;
                    ddSegments->modify();
                    this->Interactor->Render();
                }
                if(key == "4")
                {
                    DislocationSegmentActor::clr=DislocationSegmentActor::colorBurgers;
                    std::cout<<"DislocationSegment color scheme = Burgers. Reload frame to update colors."<<std::endl;
//                    ddSegments->modify();
//                    this->Interactor->Render();
                }
                if(key == "5")
                {
                    DislocationSegmentActor::clr=DislocationSegmentActor::colorSessile;
                    std::cout<<"DislocationSegment color scheme = Glissile/Sessile. Reload frame to update colors."<<std::endl;
                    //                    ddSegments->modify();
                    //                    this->Interactor->Render();
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
            
            
            if(selectedKey=="s")
            {
                if(key == "space")
                {
                    saveImage=!saveImage;
                    if(saveImage)
                    {
                        std::cout<<"saving images ON"<<std::endl;
                    }
                    else
                    {
                        std::cout<<"saving images OFF"<<std::endl;
                    }
                }
                
                if(key == "1" && saveImage)
                {
                    imageType=1;
                    std::cout<<"selecting png output"<<std::endl;
                }
                
                if(key == "2" && saveImage)
                {
                    imageType=2;
                    std::cout<<"selecting jpg output"<<std::endl;
                }
                
                if(key == "3" && saveImage)
                {
                    imageType=3;
                    std::cout<<"selecting bmp output"<<std::endl;
                }
                
                if(key == "0" && saveImage)
                {
                    imageTransparentBackground=!imageTransparentBackground;
                    if(imageTransparentBackground)
                    {
                        std::cout<<"transparent background ON"<<std::endl;
                    }
                    else
                    {
                        std::cout<<"transparent background OFF"<<std::endl;
                    }
                }
                
                if(key == "equal" && saveImage)
                {
                    imageMagnification++;
                    std::cout<<"image magnification "<<imageMagnification<<"X"<<std::endl;
                    
                }
                
                if(key == "minus" && saveImage)
                {
                    imageMagnification--;
                    if(imageMagnification<1)
                    {
                        imageMagnification=1;
                    }
                    std::cout<<"image magnification "<<imageMagnification<<"X"<<std::endl;
                }
                
            }
            
            if(selectedKey=="v")
            {
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
            
            if(selectedKey=="t")
            {
                if(key == "equal" && ddSegments.get()!=nullptr)
                {
                    DislocationSegmentActor::slippedAreaOpacity*=1.25;
                    ddSegments->modify();
                    this->Interactor->Render();
                }
                if(key == "minus" && ddSegments.get()!=nullptr)
                {
                    DislocationSegmentActor::slippedAreaOpacity/=1.25;
                    ddSegments->modify();
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







