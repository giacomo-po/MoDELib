/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2017 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DDvtk_H_
#define model_DDvtk_H_

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

#include <vtkJPEGReader.h>
#include <vtkImageMapper.h>
#include <vtkActor2D.h>

#include <vtkContextView.h>
#include <vtkContextScene.h>
#include <vtkPen.h>
#include <vtkXYPlotActor.h>
#include <vtkAxisActor2D.h>
#include <vtkPolyDataMapper2D.h>
#include <vtkProperty2D.h>
#include <vtkPoints2D.h>
#include <vtkDelimitedTextReader.h>
#include <vtkTextProperty.h>
#include <vtkAxesActor.h>
#include <vtkOrientationMarkerWidget.h>
#include <vtkContext2D.h>
#include <vtkBrush.h>
#include <vtkSliderWidget.h>
#include <vtkSliderRepresentation2D.h>

//#include <vtkProperty.h>
//#include <vtkPropPicker.h>
//#include <vtkTubeFilter.h>
//#include <vtkMath.h>
//#include <vtkObjectFactory.h> //vtkStandardNewMacro



#include <SimplicialMeshActor.h>
#include <DislocationSegmentActor.h>
//#include <DislocationActors.h>
#include <PlotActor.h>
#include <DDinteractionStyle.h>

//#include <EigenDataReader.h>
#include <TextFileParser.h>
//#include <FieldActor.h>

//#include <VertexReader.h>


namespace model
{
    
    class vtkSliderCallback : public vtkCommand
    {
    public:
        static vtkSliderCallback *New()
        {
            return new vtkSliderCallback;
        }
        virtual void Execute(vtkObject *caller, unsigned long, void*)
        {
            vtkSliderWidget *sliderWidget =            reinterpret_cast<vtkSliderWidget*>(caller);
//            this->SphereSource->SetPhiResolution(static_cast<vtkSliderRepresentation *>(sliderWidget->GetRepresentation())->GetValue());
//            this->SphereSource->SetThetaResolution(static_cast<vtkSliderRepresentation *>(sliderWidget->GetRepresentation())->GetValue());
        }
//        vtkSliderCallback():SphereSource(0) {}
//        vtkSphereSource *SphereSource;
    };
    
    
    /**************************************************************************/
    /**************************************************************************/
    struct DDvtk
    {
        
        vtkSmartPointer<vtkRenderWindow> renderWindow;
        vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor;
//        vtkSmartPointer<vtkRenderer> ddRenderer;
        vtkSmartPointer<vtkRenderer> plotRenderer;
        vtkSmartPointer<model::DDinteractionStyle> style;

//        vtkSmartPointer<vtkAxesActor> axes;
        
        
        /**********************************************************************/
        DDvtk() :
        /* init */ renderWindow(vtkSmartPointer<vtkRenderWindow>::New()),
        /* init */ renderWindowInteractor(vtkSmartPointer<vtkRenderWindowInteractor>::New()),
//        /* init */ ddRenderer(vtkSmartPointer<vtkRenderer>::New()),
        /* init */ style(vtkSmartPointer<model::DDinteractionStyle>::New()),
        /* init */ plotRenderer(vtkSmartPointer<vtkRenderer>::New())
//        axes(vtkSmartPointer<vtkAxesActor>::New())
        {
            

//            vtkSmartPointer<FieldActor> myCallback = vtkSmartPointer<FieldActor>::New();
////            myCallback->mesh=&style->meshActor.mesh;
////            myCallback->plane = style->meshActor.clipPlane;
////            myCallback->ddSegments = &style->ddConfig.segments;
////            myCallback->nodeMap=&style->ddConfig.nodeMap();
//            vtkSmartPointer<vtkImplicitPlaneRepresentation> rep = vtkSmartPointer<vtkImplicitPlaneRepresentation>::New();
//            rep->SetPlaceFactor(0.15); // This must be set prior to placing the widget
//            rep->PlaceWidget(style->meshActor.clipActor->GetBounds());
//            rep->SetNormal(style->meshActor.clipPlane->GetNormal());
//            vtkSmartPointer<vtkImplicitPlaneWidget2> planeWidget = vtkSmartPointer<vtkImplicitPlaneWidget2>::New();
//            planeWidget->SetInteractor(renderWindowInteractor);
//            planeWidget->SetRepresentation(rep);
//            planeWidget->AddObserver(vtkCommand::InteractionEvent,myCallback);
//            style->ddRenderer->AddActor(myCallback->meshActor);
//            planeWidget->On();

            
            
//            planeWidget->OutsideBoundsOff();
//
//            vtkSmartPointer<vtkSliderRepresentation2D> sliderRep =  vtkSmartPointer<vtkSliderRepresentation2D>::New();
//            sliderRep->SetMinimumValue(1.0);
//            sliderRep->SetMaximumValue(100.0);
//            sliderRep->SetValue(10.0);
//            sliderRep->SetTitleText("Element Size");
//
//            sliderRep->GetSliderProperty()->SetColor(1,0,0);//red
//
//            // Change the color of the text indicating what the slider controls
//            sliderRep->GetTitleProperty()->SetColor(1,0,0);//red
//
//            // Change the color of the text displaying the value
//            sliderRep->GetLabelProperty()->SetColor(1,0,0);//red
//
//            // Change the color of the knob when the mouse is held on it
//            sliderRep->GetSelectedProperty()->SetColor(0,1,0);//green
//
//            // Change the color of the bar
//            sliderRep->GetTubeProperty()->SetColor(1,1,0);//yellow
//
//            // Change the color of the ends of the bar
//            sliderRep->GetCapProperty()->SetColor(1,1,0);//yellow
//            sliderRep->GetPoint1Coordinate()->SetCoordinateSystemToDisplay();
//            sliderRep->GetPoint1Coordinate()->SetValue(40 ,40);
//            sliderRep->GetPoint2Coordinate()->SetCoordinateSystemToDisplay();
//            sliderRep->GetPoint2Coordinate()->SetValue(100, 40);
//
//            vtkSmartPointer<vtkSliderWidget> sliderWidget =
//            vtkSmartPointer<vtkSliderWidget>::New();
//            sliderWidget->SetInteractor(renderWindowInteractor);
//            sliderWidget->SetRepresentation(sliderRep);
//            sliderWidget->SetAnimationModeToAnimate();
//            sliderWidget->EnabledOn();
//            vtkSmartPointer<vtkSliderCallback> callback =
//            vtkSmartPointer<vtkSliderCallback>::New();
////            callback->SphereSource = sphereSource;
//
//            sliderWidget->AddObserver(vtkCommand::InteractionEvent,callback);
//            renderWindowInteractor->Initialize();

            


            
//            int meshID=TextFileParser("inputFiles/DD.txt").readScalar<int>("meshID",false);
//            model::EigenDataReader EDR;
//            bool use_boundary=false;
//            EDR.readScalarInFile("./DDinput.txt","use_boundary",use_boundary);
//            if (use_boundary)
//            {
//                EDR.readScalarInFile("./DDinput.txt","meshID",meshID);
//            }

            // https://en.wikipedia.org/wiki/Computer_display_standard
//            renderWindow->SetSize(1024,768); // XGA (width, height)
            renderWindow->SetSize(1920,1080); // HD (width, height)

            renderWindowInteractor->SetRenderWindow(renderWindow);
//            ddRenderer->SetBackground(1,1,1); // Background color white
//            ddRenderer->SetViewport(0.0,0,0.5,1);
            
//            renderWindow->SetAlphaBitPlanes(true);
//            // 2. Force to not pick a framebuffer with a multisample buffer
//            // (as initial value is 8):
//            renderWindow->SetMultiSamples(0);
//            
//            // 3. Choose to use depth peeling (if supported) (initial value is 0 (false)):
//            ddRenderer->SetUseDepthPeeling(true);
//            
//            // 4. Set depth peeling parameters
//            // - Set the maximum number of rendering passes (initial value is 4):
//            ddRenderer->SetMaximumNumberOfPeels(20);
//            // - Set the occlusion ratio (initial value is 0.0, exact image):
//            ddRenderer->SetOcclusionRatio(0.5);


            
            renderWindow->AddRenderer(style->ddRenderer);

            renderWindow->AddRenderer(plotRenderer);

            renderWindowInteractor->SetInteractorStyle(style);

//            style->SetDefaultRenderer(ddRenderer);
            
            style->init(plotRenderer);


            
            // Start
            style->ddRenderer->ResetCamera();
            renderWindow->LineSmoothingOn();
            renderWindow->PolygonSmoothingOn();
            renderWindow->PointSmoothingOn();
            renderWindow->SetMultiSamples(1);
            renderWindow->Render();
            

//            planeWidget->On();

            renderWindowInteractor->Start();
            

        }
        
    };
    
    
}
#endif

