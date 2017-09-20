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

#include <vtkChartXY.h>
#include <vtkTable.h>
#include <vtkPlot.h>
#include <vtkFloatArray.h>
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

//#include <vtkProperty.h>
//#include <vtkPropPicker.h>
//#include <vtkTubeFilter.h>
//#include <vtkMath.h>
//#include <vtkObjectFactory.h> //vtkStandardNewMacro



#include <model/Mesh/SimplicialMesh.h>
#include <model/DislocationDynamics/Visualization/DDvtk/SimplicialMeshActor.h>
#include <model/DislocationDynamics/Visualization/DDvtk/DislocationSegmentActor.h>
//#include <model/DislocationDynamics/Visualization/vtk/DislocationActors.h>
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
        
        /**********************************************************************/
        DDvtk()
        {
            vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();
            
            vtkSmartPointer<vtkFloatArray> arrX = vtkSmartPointer<vtkFloatArray>::New();
            arrX->SetName("disp");
            table->AddColumn(arrX);
            
            vtkSmartPointer<vtkFloatArray> arrC = vtkSmartPointer<vtkFloatArray>::New();
            arrC->SetName("stress");
            table->AddColumn(arrC);
            
            //    vtkSmartPointer<vtkFloatArray> arrS = vtkSmartPointer<vtkFloatArray>::New();
            //    arrS->SetName("Sine");
            //    table->AddColumn(arrS);
            
            
            model::IDreader<'F',1,200,double> vReader;
            Eigen::Matrix<double,1,200> temp(Eigen::Matrix<double,1,200>::Zero());
            
            
            if (vReader.isGood(0,true))
            {
                vReader.read(0,true);
                
                table->SetNumberOfRows(vReader.size());
                int i=0;
                int xCcol=0;
                int yCol=1;
                
                for (const auto& row : vReader)
                {
                    table->SetValue(i, 0, row.second[xCcol]);
                    table->SetValue(i, 1, row.second[yCol]);
                    
                    i++;
                }
                
                //        vReader.read(0,true);
                //
                //        if(runID<0)
                //        {
                //            if(vReader.size())
                //            {
                //                runID=vReader.rbegin()->first;
                //                temp=vReader.rbegin()->second;
                //            }
                //            else
                //            {
                //                runID=0;
                //            }
                //        }
                //        else
                //        {
                //            const auto iter=vReader.find(runID);
                //            if(iter!=vReader.end())
                //            {// runID has been found
                //                temp=iter->second;
                //            }
                //            else
                //            {
                //                assert(0 && "runID NOT FOUND IN F/F_0.txt");
                //            }
                //        }
            }
            else
            {
                model::cout<<"could not read runID from F/F_0.txt"<<std::endl;
                //        runID=0;
            }
            
            // Create a table with some points in it
            
            // Fill in the table with some example values
            //    int numPoints = 69;
            //    float inc = 7.5 / (numPoints-1);
            //    table->SetNumberOfRows(numPoints);
            //    for (int i = 0; i < numPoints; ++i)
            //    {
            ////        table->SetValue(i, 0, i * inc);
            ////        table->SetValue(i, 1, cos(i * inc));
            ////        table->SetValue(i, 2, sin(i * inc));
            //    }
            
            // Set up the view
            vtkSmartPointer<vtkContextView> view = vtkSmartPointer<vtkContextView>::New();
            view->GetRenderWindow()->SetSize(1024,768);
            view->GetRenderer()->SetBackground(1.0, 1.0, 1.0);
            view->GetRenderer()->SetViewport(0.5,0,1.0,1);
            
            // Add multiple line plots, setting the colors etc
            vtkSmartPointer<vtkChartXY> chart = vtkSmartPointer<vtkChartXY>::New();
            view->GetScene()->AddItem(chart);
            
            vtkPlot *line = chart->AddPlot(vtkChart::LINE);
            line->SetInputData(table, 0, 1);
            line->SetColor(0, 0, 255, 255);
            line->SetWidth(1.0);
            
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
            
            
            
            // Define the actors
            int meshID(0);
            model::EigenDataReader EDR;
            bool use_boundary=false;
            EDR.readScalarInFile("./DDinput.txt","use_boundary",use_boundary);
            if (use_boundary)
            {
                EDR.readScalarInFile("./DDinput.txt","meshID",meshID);
            }
            
            // Define first renderer
            vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
            renderer->SetBackground(1,1,1); // Background color white
            renderer->SetViewport(0.0,0,0.5,1);
            view->GetRenderWindow()->AddRenderer(renderer);
            
            vtkSmartPointer<model::DDinteractionStyle> style =vtkSmartPointer<model::DDinteractionStyle>::New();
            view->GetInteractor()->SetInteractorStyle(style);
            style->SetDefaultRenderer(renderer);
            style->meshActor.init(meshID,renderer);
//            style->loadFrame(0);
            
            
            vtkSmartPointer<vtkAxesActor> axes = vtkSmartPointer<vtkAxesActor>::New();
            
            vtkSmartPointer<vtkOrientationMarkerWidget> widget = vtkSmartPointer<vtkOrientationMarkerWidget>::New();
            widget->SetOutlineColor( 0.9300, 0.5700, 0.1300 );
            widget->SetOrientationMarker( axes );
            widget->SetInteractor( view->GetInteractor() );
            widget->SetViewport( 0.0, 0.0, 0.4, 0.4 );
            widget->SetEnabled( 1 );
            widget->InteractiveOn();
            
            renderer->ResetCamera();
            view->GetRenderWindow()->Render();
            
            
            // Start interactor
            view->GetInteractor()->Initialize();
            view->GetInteractor()->Start();

        }
        
    };
    
    
}
#endif







