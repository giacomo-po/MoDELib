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

#include <model/IO/EigenDataReader.h>
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


        
        /**********************************************************************/
        DDvtk() :
        /* init */ renderWindow(vtkSmartPointer<vtkRenderWindow>::New()),
        /* init */ renderWindowInteractor(vtkSmartPointer<vtkRenderWindowInteractor>::New()),
        /* init */ ddRenderer(vtkSmartPointer<vtkRenderer>::New()),
        /* init */ style(vtkSmartPointer<model::DDinteractionStyle>::New()),
        /* init */ plotRenderer(vtkSmartPointer<vtkRenderer>::New())
        {
            
            int meshID(0);
            model::EigenDataReader EDR;
            bool use_boundary=false;
            EDR.readScalarInFile("./DDinput.txt","use_boundary",use_boundary);
            if (use_boundary)
            {
                EDR.readScalarInFile("./DDinput.txt","meshID",meshID);
            }

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
            renderWindowInteractor->SetInteractorStyle(style);
            style->SetDefaultRenderer(ddRenderer);
//            style->SetCurrentRenderer(ddRenderer);
            style->ddRenderer=ddRenderer;
            style->plotRenderer=plotRenderer;


            style->meshActor.init(meshID,ddRenderer);
            style->loadFrame(0);
            
            plotRenderer->SetBackground(1,1,1);
            plotRenderer->SetViewport(0.5,0,1.0,1);
            
            PlotActor pa(plotRenderer);
            
            renderWindow->AddRenderer(plotRenderer);
            

            
            
            // Start
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




//    line = chart->AddPlot(vtkChart::LINE);
//    line->SetInputData(table, 0, 2);
//    line->SetColor(255, 0, 0, 255);
//    line->SetWidth(5.0);

// For dotted line, the line type can be from 2 to 5 for different dash/dot
// patterns (see enum in vtkPen containing DASH_LINE, value 2):
//#ifndef WIN32
//    line->GetPen()->SetLineType(vtkPen::DASH_LINE);
//#endif
// (ifdef-ed out on Windows because DASH_LINE does not work on Windows
//  machines with built-in Intel HD graphics card...)

//view->GetRenderWindow()->SetMultiSamples(0);
//
//
//
//            // Define the actors
//            int meshID(0);
//            model::EigenDataReader EDR;
//            bool use_boundary=false;
//            EDR.readScalarInFile("./DDinput.txt","use_boundary",use_boundary);
//            if (use_boundary)
//            {
//                EDR.readScalarInFile("./DDinput.txt","meshID",meshID);
//            }
//
//            // Define first renderer
//            vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
//            renderer->SetBackground(1,1,1); // Background color white
//            renderer->SetViewport(0.0,0,0.5,1);
//            view->GetRenderWindow()->AddRenderer(renderer);
//
//            vtkSmartPointer<model::DDinteractionStyle> style =vtkSmartPointer<model::DDinteractionStyle>::New();
//            view->GetInteractor()->SetInteractorStyle(style);
//            style->SetDefaultRenderer(renderer);
//            style->meshActor.init(meshID,renderer);
////            style->loadFrame(0);
//
//
//            vtkSmartPointer<vtkAxesActor> axes = vtkSmartPointer<vtkAxesActor>::New();
//
//            vtkSmartPointer<vtkOrientationMarkerWidget> widget = vtkSmartPointer<vtkOrientationMarkerWidget>::New();
//            widget->SetOutlineColor( 0.9300, 0.5700, 0.1300 );
//            widget->SetOrientationMarker( axes );
//            widget->SetInteractor( view->GetInteractor() );
//            widget->SetViewport( 0.0, 0.0, 0.4, 0.4 );
//            widget->SetEnabled( 1 );
//            widget->InteractiveOn();
//
//            renderer->ResetCamera();
//            view->GetRenderWindow()->Render();
//
//
//            // Start interactor
//            view->GetInteractor()->Initialize();
//            view->GetInteractor()->Start();

// Plot section
//            vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();
//
//            vtkSmartPointer<vtkFloatArray> arrX = vtkSmartPointer<vtkFloatArray>::New();
//            arrX->SetName("disp");
//            table->AddColumn(arrX);
//
//            vtkSmartPointer<vtkFloatArray> arrC = vtkSmartPointer<vtkFloatArray>::New();
//            arrC->SetName("stress");
//            table->AddColumn(arrC);
//
//            //    vtkSmartPointer<vtkFloatArray> arrS = vtkSmartPointer<vtkFloatArray>::New();
//            //    arrS->SetName("Sine");
//            //    table->AddColumn(arrS);
//
//
//            model::IDreader<'F',1,200,double> vReader;
//            Eigen::Matrix<double,1,200> temp(Eigen::Matrix<double,1,200>::Zero());
//
//
//            if (vReader.isGood(0,true))
//            {
//                vReader.read(0,true);
//
//                table->SetNumberOfRows(vReader.size());
//                int i=0;
//                int xCcol=0;
//                int yCol=1;
//
//                for (const auto& row : vReader)
//                {
//                    table->SetValue(i, 0, row.second[xCcol]);
//                    table->SetValue(i, 1, row.second[yCol]);
//
//                    i++;
//                }
//            }
//            else
//            {
//                model::cout<<"could not read runID from F/F_0.txt"<<std::endl;
//                //        runID=0;
//            }
//
//
//
//            // Add multiple line plots, setting the colors etc
//            vtkSmartPointer<vtkChartXY> chart = vtkSmartPointer<vtkChartXY>::New();
//
//            vtkPlot *line = chart->AddPlot(vtkChart::LINE);
//            line->SetInputData(table, 0, 1);
//            line->SetColor(0, 0, 255, 255);
//            line->SetWidth(1.0);
//
//            VTK_CREATE(vtkContextActor, actor);
//            VTK_CREATE(APIDiagram, diagram);
//            actor->GetScene()->AddItem(chart);
//            plotRenderer->AddActor(actor);



