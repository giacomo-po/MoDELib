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
    

    
    // Define the actors
    int meshID(0);
    model::EigenDataReader EDR;
    bool use_boundary=false;
    EDR.readScalarInFile("./DDinput.txt","use_boundary",use_boundary);
    if (use_boundary)
    {
        EDR.readScalarInFile("./DDinput.txt","meshID",meshID);
    }
    
    // Define a render window
    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->SetSize(1024,768); //(width, height)

    // Define the  window interactor
    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    renderWindowInteractor->SetRenderWindow(renderWindow);

    // Define first renderer
    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    renderer->SetBackground(1,1,1); // Background color white
    renderWindow->AddRenderer(renderer);

    // Set the custom type to use for interaction.
    vtkSmartPointer<model::DDinteractionStyle> style =vtkSmartPointer<model::DDinteractionStyle>::New();
    renderWindowInteractor->SetInteractorStyle(style);
    style->SetDefaultRenderer(renderer);
    
    // Pupulate initial actors
    style->meshActor.init(meshID,renderer);
    style->ddActors.update(0,renderer);

    

    
    if(true)
    {
        // Define second renderer
        vtkSmartPointer<vtkRenderer> renderer2 = vtkSmartPointer<vtkRenderer>::New();
        //renderer2->SetBackground(0.9,1,1); // Background color white
        renderWindow->AddRenderer(renderer2);
        renderer2->SetViewport(0.5,0,1,1);
        renderer2->SetBackground(1, 1, 1);
        
        renderer->SetViewport(0.0,0,0.5,1);
        

        
        std::string inputFilename="F/F_0.txt";
        vtkSmartPointer<vtkDelimitedTextReader> reader =  vtkSmartPointer<vtkDelimitedTextReader>::New();
        reader->SetFileName(inputFilename.c_str());
        reader->DetectNumericColumnsOn();
        reader->SetFieldDelimiterCharacters(" ");
        reader->Update();
        
        vtkTable* table = reader->GetOutput();
        std::cout << "Table has " << table->GetNumberOfRows() << " rows." << std::endl;
        std::cout << "Table has " << table->GetNumberOfColumns() << " columns." << std::endl;
        
        
        
                          int DIM = table->GetNumberOfRows();
                          vtkDataArray *dataArray1 = vtkDataArray::CreateDataArray(VTK_FLOAT);
                          dataArray1->SetNumberOfTuples(DIM);
        
                          vtkDataArray *dataArray2 =
                       vtkDataArray::CreateDataArray(VTK_FLOAT);
                          dataArray2->SetNumberOfTuples(DIM);
        
                          int t;
        int c1=0;
        int c2=2;
                          for (t = 0; t < DIM; t++)
                           {
                                    float x = (table->GetValue(t,c1)).ToDouble();
                                    float y = (table->GetValue(t,c2)).ToDouble();
                                    dataArray1->SetTuple(t, &x);
                                    dataArray2->SetTuple(t, &y);
                                }
        
                          vtkFieldData *fieldData = vtkFieldData::New();
                          fieldData->AllocateArrays(2);
                          fieldData->AddArray(dataArray1);
                          fieldData->AddArray(dataArray2);
        
                          vtkDataObject *dataObject = vtkDataObject::New();
                          dataObject->SetFieldData(fieldData);
        
                          vtkXYPlotActor *plot = vtkXYPlotActor::New();
                          plot->AddDataObjectInput(dataObject);
//                          plot->SetTitle("Plot");
                          plot->SetXTitle("X-Axis");
                          plot->SetYTitle("Y-Axis");
//                          plot->SetXValuesToValue(); // not sure what this does
                          plot->SetWidth(0.9);
                          plot->SetHeight(0.9);
                          plot->SetPosition(0.05, 0.05);
                          plot->LegendOn();
                          plot->PickableOff();
                          plot->PlotCurvePointsOn();
                          plot->PlotCurveLinesOff();
        
                          plot->SetDataObjectXComponent(0, 0);
                          plot->SetDataObjectYComponent(0, 1);
                          plot->SetPlotColor(0, 0.0, 0.0, 1.0);
                          plot->SetPlotLabel(0, "My Label");
                          //plot->GetProperty()->SetColor(0.0, 0.0, 0.0);
        
        plot->GetXAxisActor2D()->GetLabelTextProperty()->SetColor(0.0, 0.0, 0.0);
        plot->GetYAxisActor2D()->GetLabelTextProperty()->SetOrientation(90.0);
        
                          renderer2->AddActor2D(plot);
//        renderer2->Render();
    }

    if(false) // this works to add an image
        {
            vtkSmartPointer<vtkJPEGReader> pngReader = vtkSmartPointer<vtkJPEGReader>::New();
            pngReader->SetFileName("IMG_0013.jpg");
            pngReader->Update();
            
            vtkSmartPointer<vtkImageMapper> imageMapper = vtkSmartPointer<vtkImageMapper>::New();
            imageMapper->SetInputConnection(pngReader->GetOutputPort());
            imageMapper->SetColorWindow(255);
            imageMapper->SetColorLevel(127.5);
            
            vtkSmartPointer<vtkActor2D> imageActor = vtkSmartPointer<vtkActor2D>::New();
            imageActor->SetMapper(imageMapper);
            imageActor->SetPosition2(0.5,0.5);
            
            //vtkRenderer renderer = renderWindowControl1.RenderWindow.GetRenderers().GetFirstRenderer();
            renderer->AddActor2D(imageActor);
        }

    if(false) // try to add stress-strain plot
    {
        
        //
        //    vtkSmartPointer<vtkAxisActor2D> axisActor=vtkSmartPointer<vtkAxisActor2D>::New();
        //    axisActor->SetRange(1000,1000);
        //    axisActor->AxisVisibilityOn();
        //    renderer->AddActor2D(axisActor);
        
        
        //    vtkSmartPointer<vtkPoints> points=vtkSmartPointer<vtkPoints>::New();
        //    points->InsertPoint(0, 0, 0,0);
        //    points->InsertPoint(1, 1000, 1000,0);
        
        //    vtkSmartPointer<vtkCellArray> lines=vtkSmartPointer<vtkCellArray>::New();
        //    lines->InsertCellPoint(0);
        //    lines->InsertCellPoint(1);
        
        // http://www.vtk.org/Wiki/VTK/Examples/Cxx/InfoVis/DelimitedTextReader
        std::string inputFilename="F/F_0.txt";
        vtkSmartPointer<vtkDelimitedTextReader> reader =  vtkSmartPointer<vtkDelimitedTextReader>::New();
        reader->SetFileName(inputFilename.c_str());
        reader->DetectNumericColumnsOn();
        reader->SetFieldDelimiterCharacters(" ");
        reader->Update();
        
        vtkTable* table = reader->GetOutput();
        std::cout << "Table has " << table->GetNumberOfRows() << " rows." << std::endl;
        std::cout << "Table has " << table->GetNumberOfColumns() << " columns." << std::endl;
        
        // http://public.kitware.com/pipermail/vtkusers/2009-February/050753.html
        vtkSmartPointer<vtkXYPlotActor> actor = vtkSmartPointer<vtkXYPlotActor>::New();
        //actor->AddDataSetInput(table, 0, 3);
        
        //vtkSmartPointer<vtkChartXY> chart = vtkSmartPointer<vtkChartXY>::New();
        //view->GetScene()->AddItem(chart);
        //chart->SetShowLegend(true);
        
        //chart->AddPlot(vtkChart::POINTS)
        
        //        vtkPlot *points = chart->AddPlot(vtkChart::POINTS);
        //        points->SetInputData(table, 0, 1);
        //        points = chart->AddPlot(vtkChart::POINTS);
        //        points->SetColor(0, 0, 0, 255);
        //        points->SetWidth(1.0);
        //        vtkPlotPoints::SafeDownCast(points)->SetMarkerStyle(vtkPlotPoints::CROSS);
        //
        //        points = chart->AddPlot(vtkChart::POINTS);
        //        points->SetInputData(table, 0, 2);
        //        points->SetColor(0, 0, 0, 255);
        //        points->SetWidth(1.0);
        //        vtkPlotPoints::SafeDownCast(points)->SetMarkerStyle(vtkPlotPoints::PLUS);
        //
        //        points = chart->AddPlot(vtkChart::POINTS);
        //        points->SetInputData(table, 0, 3);
        //        points->SetColor(0, 0, 255, 255);
        //        points->SetWidth(1.0);
        //        vtkPlotPoints::SafeDownCast(points)->SetMarkerStyle(vtkPlotPoints::CIRCLE);
        //
        //
        //        vtkSmartPointer<vtkPoints> points=vtkSmartPointer<vtkPoints>::New();
        //
        //
        //        for(vtkIdType i = 0; i < table->GetNumberOfRows(); i++)
        //        {
        //            std::cout << "x: " << (table->GetValue(i,0)).ToDouble()
        //            << " y: " << (table->GetValue(i,1)).ToDouble()
        //            << " z: " << (table->GetValue(i,2)).ToDouble();
        //
        //            points->InsertNextPoint((table->GetValue(i,0)).ToDouble(),
        //                                    (table->GetValue(i,1)).ToDouble(),
        //                                    (table->GetValue(i,2)).ToDouble());
        //
        ////            double n[3];
        ////            n[0] = (table->GetValue(i,3)).ToDouble();
        ////            n[1] = (table->GetValue(i,4)).ToDouble();
        ////            n[2] = (table->GetValue(i,5)).ToDouble();
        ////
        ////            std::cout << " n: " << n[0] << " " << n[1] << " " << n[2] << std::endl;
        ////            normals->InsertNextTuple(n);
        //        }
        //
        //        std::cout << "There are " << points->GetNumberOfPoints()
        //        << " points." << std::endl;
        //
        //    // http://www.vtk.org/Wiki/VTK/Examples/Cxx/Visualization/Visualize2DPoints
        //    vtkSmartPointer<vtkPolyData> polyData= vtkSmartPointer<vtkPolyData>::New();
        //
        //    polyData->SetPoints(points);
        ////    polyData->SetLines(lines);
        //
        //    // populate polydata
        //    vtkSmartPointer<vtkPolyDataMapper2D> mapper = vtkSmartPointer<vtkPolyDataMapper2D>::New();
        //    mapper->SetInputData( polyData );
        //    mapper->ScalarVisibilityOff();
        //
        //    
        //    double color[3] = {1,0,0};
        
        //    actor->SetMapper( mapper );
        //    actor->GetProperty()->SetColor( color );
        //    actor->GetProperty()->SetPointSize( 200.0 );
        
        renderer->AddActor2D(actor);
    }
    

    // Start
    renderWindow->LineSmoothingOn();
    renderWindow->PolygonSmoothingOn();
    renderWindow->PointSmoothingOn();
    renderWindow->SetMultiSamples(1);
    
    
//    vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();
//    
//    vtkSmartPointer<vtkFloatArray> arrX = vtkSmartPointer<vtkFloatArray>::New();
//    arrX->SetName("X Axis");
//    table->AddColumn(arrX);
//    
//    vtkSmartPointer<vtkFloatArray> arrC = vtkSmartPointer<vtkFloatArray>::New();
//    arrC->SetName("Cosine");
//    table->AddColumn(arrC);
//    
//    vtkSmartPointer<vtkFloatArray> arrS = vtkSmartPointer<vtkFloatArray>::New();
//    arrS->SetName("Sine");
//    table->AddColumn(arrS);
//    
//    // Fill in the table with some example values
//    int numPoints = 69;
//    float inc = 7.5 / (numPoints-1);
//    table->SetNumberOfRows(numPoints);
//    for (int i = 0; i < numPoints; ++i)
//    {
//        table->SetValue(i, 0, i * inc);
//        table->SetValue(i, 1, cos(i * inc));
//        table->SetValue(i, 2, sin(i * inc));
//    }
//    
//    // Set up the view
//    vtkSmartPointer<vtkContextView> view = vtkSmartPointer<vtkContextView>::New();
//    view->SetRenderWindow(renderWindow);
////    view->GetRenderer()->SetBackground(1.0, 1.0, 1.0);
//    
//    // Add multiple line plots, setting the colors etc
//    vtkSmartPointer<vtkChartXY> chart =
//    vtkSmartPointer<vtkChartXY>::New();
//    view->GetScene()->AddItem(chart);
//    vtkPlot *line = chart->AddPlot(vtkChart::LINE);
//#if VTK_MAJOR_VERSION <= 5
//    line->SetInput(table, 0, 1);
//#else
//    line->SetInputData(table, 0, 1);
//#endif
//    line->SetColor(0, 255, 0, 255);
//    line->SetWidth(1.0);
//    line = chart->AddPlot(vtkChart::LINE);
//#if VTK_MAJOR_VERSION <= 5
//    line->SetInput(table, 0, 2);
//#else
//    line->SetInputData(table, 0, 2);
//#endif
//    line->SetColor(255, 0, 0, 255);
//    line->SetWidth(5.0);
//    
//    // For dotted line, the line type can be from 2 to 5 for different dash/dot
//    // patterns (see enum in vtkPen containing DASH_LINE, value 2):
//#ifndef WIN32
//    line->GetPen()->SetLineType(vtkPen::DASH_LINE);
//#endif
//    // (ifdef-ed out on Windows because DASH_LINE does not work on Windows
//    //  machines with built-in Intel HD graphics card...)
//    
//    //view->GetRenderWindow()->SetMultiSamples(0);
//    
//    // Start interactor
//    view->GetInteractor()->Initialize();

    

    renderWindow->Render();
    renderWindowInteractor->Start();

    return EXIT_SUCCESS;
}

