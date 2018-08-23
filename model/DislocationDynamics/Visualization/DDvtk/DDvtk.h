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

//#include <vtkProperty.h>
//#include <vtkPropPicker.h>
//#include <vtkTubeFilter.h>
//#include <vtkMath.h>
//#include <vtkObjectFactory.h> //vtkStandardNewMacro



#include <model/DislocationDynamics/Visualization/DDvtk/SimplicialMeshActor.h>
#include <model/DislocationDynamics/Visualization/DDvtk/DislocationSegmentActor.h>
//#include <model/DislocationDynamics/Visualization/vtk/DislocationActors.h>
#include <model/DislocationDynamics/Visualization/DDvtk/PlotActor.h>
#include <model/DislocationDynamics/Visualization/DDvtk/DDinteractionStyle.h>

//#include <model/IO/EigenDataReader.h>
#include <model/IO/TextFileParser.h>

//#include <model/IO/VertexReader.h>
#include <model/IO/IDreader.h>


namespace model
{
    
    /**************************************************************************/
    /**************************************************************************/
    struct DDvtk
    {
        
        vtkSmartPointer<vtkRenderWindow> renderWindow;
        vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor;
        vtkSmartPointer<vtkRenderer> ddRenderer;
        vtkSmartPointer<vtkRenderer> plotRenderer;
        vtkSmartPointer<model::DDinteractionStyle> style;

//        vtkSmartPointer<vtkAxesActor> axes;
        
        
        /**********************************************************************/
        DDvtk() :
        /* init */ renderWindow(vtkSmartPointer<vtkRenderWindow>::New()),
        /* init */ renderWindowInteractor(vtkSmartPointer<vtkRenderWindowInteractor>::New()),
        /* init */ ddRenderer(vtkSmartPointer<vtkRenderer>::New()),
        /* init */ style(vtkSmartPointer<model::DDinteractionStyle>::New()),
        /* init */ plotRenderer(vtkSmartPointer<vtkRenderer>::New())
//        axes(vtkSmartPointer<vtkAxesActor>::New())
        {
            
            int meshID=TextFileParser("inputFiles/DD.txt").readScalar<int>("meshID",false);
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
            ddRenderer->SetBackground(1,1,1); // Background color white
            ddRenderer->SetViewport(0.0,0,0.5,1);
            
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


            renderWindow->AddRenderer(ddRenderer);
            renderWindow->AddRenderer(plotRenderer);

            renderWindowInteractor->SetInteractorStyle(style);
            style->SetDefaultRenderer(ddRenderer);
            
            style->init(ddRenderer,plotRenderer,meshID);
//            style->ddRenderer=ddRenderer;
//            style->plotRenderer=plotRenderer;
//
//
//            style->meshActor.init(meshID,ddRenderer);
//            style->loadFrame(0);
//            
//            plotRenderer->SetBackground(1,1,1);
//            plotRenderer->SetViewport(0.5,0,1.0,1);
            
//            PlotActor pa(plotRenderer);
            
            // The axes are positioned with a user transform
//            axes->SetUserTransform(transform);

//            axes->SetTotalLength(1000,1000,1000);
//            ddRenderer->AddActor(axes);

//            vtkSmartPointer<vtkOrientationMarkerWidget> widget =
//            vtkSmartPointer<vtkOrientationMarkerWidget>::New();
//            widget->SetOutlineColor( 0.9300, 0.5700, 0.1300 );
//            widget->SetOrientationMarker( axes );
//            widget->SetInteractor( renderWindowInteractor );
//            widget->SetViewport( 0.0, 0.0, 0.4, 0.4 );
//            widget->SetEnabled( 1 );
//            widget->InteractiveOn();
            
//            renderWindow->Render();
            
            // Start
            ddRenderer->ResetCamera();
            renderWindow->LineSmoothingOn();
            renderWindow->PolygonSmoothingOn();
            renderWindow->PointSmoothingOn();
            renderWindow->SetMultiSamples(1);
            renderWindow->Render();
            renderWindowInteractor->Start();
            

        }
        
    };
    
    
}
#endif

