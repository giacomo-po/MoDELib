/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SimplicialMeshActor_cpp_
#define model_SimplicialMeshActor_cpp_


#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkTriangle.h>
#include <vtkCellData.h>
#include <vtkPlane.h>
#include <vtkClipPolyData.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkPolyDataMapper.h>
#include <vtkDataSetMapper.h>
#include <vtkPolyData.h>
#include <vtkStripper.h>
#include <vtkFeatureEdges.h>
#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkPolyDataMapper.h>
#include <vtkClipPolyData.h>
#include <vtkPlane.h>
#include <vtkSphereSource.h>
#include <vtkPlane.h>
#include <vtkCommand.h>
#include <vtkImplicitPlaneWidget2.h>
#include <vtkImplicitPlaneRepresentation.h>
#include <vtkNamedColors.h>
#include <vtkLookupTable.h>
#include <vtkUnsignedCharArray.h>
#include <vtkPointData.h>
#include <vtkNamedColors.h>
#include <vtkTextProperty.h>
#include <random>
#include <algorithm>

#include <TextFileParser.h>

#include <SimplicialMeshActor.h>


namespace model
{

        SimplicialMeshActor::SimplicialMeshActor(vtkSmartPointer<vtkGenericOpenGLRenderWindow> renderWindow_in,vtkRenderer* const renderer_in, const SimplicialMesh<3>& mesh_in) :
        /* init */ mainLayout(new QGridLayout(this))
        /* init */,showMesh(new QCheckBox(this))
//        /* init */,showMeshLabel(new QLabel("show mesh"))
        /* init */,showFaceBoundaries(new QCheckBox(this))
//        /* init */,showFaceBoundariesLabel(new QLabel("show face edges"))
//        /* init */,showGrainColors(new QCheckBox(this))
//        /* init */,showGrainColorsLabel(new QLabel("show grain colors"))
        /* init */,showRegionBoundaries(new QCheckBox(this))
        /* init */,sliderRegionBoundaries(new QSlider(this))
        /* init */,showClipPlane(new QCheckBox(this))
        /* init */,sliderClipPlane(new QSlider(this))
        /* init */,showAxes(new QCheckBox(this))
        /* init */,showPeriodicityVectors(new QCheckBox(this))
        ///* init */,sliderAxes(new QSlider(this))
//        /* init */,showRegionBoundariesLabel(new QLabel("show grain boundaries"))
        /* init */,renderWindow(renderWindow_in)
        /* init */,renderer(renderer_in)
        /* init */,mesh(mesh_in)
        /* init */,dispCorr(1.0)
        /* init */,pts(vtkSmartPointer<vtkPoints>::New())
        /* init */,polydata(vtkSmartPointer<vtkPolyData>::New())
        /* init */,mapper(vtkSmartPointer<vtkPolyDataMapper>::New())
        /* init */,actor(vtkSmartPointer<vtkActor>::New())
        /* init */,facePts(vtkSmartPointer<vtkPoints>::New())
        /* init */,facePolydata(vtkSmartPointer<vtkPolyData>::New())
        /* init */,faceMapper(vtkSmartPointer<vtkPolyDataMapper>::New())
        /* init */,faceActor(vtkSmartPointer<vtkActor>::New())
        /* init */,gbColors(vtkSmartPointer<vtkUnsignedCharArray>::New())
        /* init */,gbPoints(vtkSmartPointer<vtkPoints>::New())
        /* init */,gbTriangles(vtkSmartPointer<vtkCellArray>::New())
        /* init */,gbTrianglePolyData(vtkSmartPointer<vtkPolyData>::New())
        /* init */,gbMapper(vtkSmartPointer<vtkPolyDataMapper>::New())
        /* init */,gbActor(vtkSmartPointer<vtkActor>::New())
        /* init */,clipPlane(vtkSmartPointer<vtkPlane>::New())
        /* init */,clipper(vtkSmartPointer<vtkClipPolyData>::New())
        /* init */,clippedPolyData(vtkSmartPointer<vtkPolyData>::New())
        /* init */,clipMapper(vtkSmartPointer<vtkDataSetMapper>::New())
        /* init */,clipActor(vtkSmartPointer<vtkActor>::New())
        /* init */,cubeAxesActor(vtkSmartPointer<vtkCubeAxesActor>::New())
        /* init */,periodicityPolydata(vtkSmartPointer<vtkPolyData>::New())
        /* init */,periodicityGlyphs(vtkSmartPointer<vtkGlyph3D>::New())
        /* init */,periodicityMapper(vtkSmartPointer<vtkPolyDataMapper>::New())
        /* init */,periodicityActor(vtkSmartPointer<vtkActor>::New())

//        /* init */,axes(vtkSmartPointer<vtkAxesActor>::New())
//        /* init */,axesWidget(vtkSmartPointer<vtkOrientationMarkerWidget>::New())
//        /* init */,myCallback(vtkSmartPointer<FieldActorCallback>::New())
//        /* init */,rep(vtkSmartPointer<vtkImplicitPlaneRepresentation>::New())
//        /* init */,planeWidget(vtkSmartPointer<vtkImplicitPlaneWidget2>::New())
        {

            showFaceBoundaries->setChecked(true);
            showFaceBoundaries->setText("Face Edges");
            faceActor->SetVisibility(true);
            showMesh->setChecked(false);
            showMesh->setText("Mesh");
            actor->SetVisibility(false);
//            showGrainColors->setChecked(false);
//            showGrainColors->setText("grain colors");
//            gbMapper->SetScalarVisibility(false);
            showRegionBoundaries->setChecked(false);
            showRegionBoundaries->setText("Grain Boundaries");
            sliderRegionBoundaries->setEnabled(false);
            gbActor->SetVisibility(false);

            showClipPlane->setText("Clip Plane");
            showClipPlane->setChecked(false);
            sliderClipPlane->setEnabled(false);
            clipActor->SetVisibility(false);
            
            showAxes->setText("Axes");
            showAxes->setChecked(false);
            showAxes->setEnabled(true);
            cubeAxesActor->SetVisibility(false);
            
            showPeriodicityVectors->setText("Periodicity Vectors");
            periodicityGlyphs->SetInputData(periodicityPolydata);
            periodicityGlyphs->SetSourceConnection(vtkSmartPointer<vtkArrowSource>::New()->GetOutputPort());
            periodicityGlyphs->ScalingOn();
            periodicityGlyphs->SetScaleModeToScaleByVector();
            periodicityGlyphs->OrientOn();
            periodicityGlyphs->ClampingOff();
            periodicityGlyphs->SetVectorModeToUseVector();
            periodicityGlyphs->SetIndexModeToOff();
            periodicityMapper->SetInputConnection(periodicityGlyphs->GetOutputPort());
            periodicityMapper->ScalarVisibilityOff();
            periodicityActor->SetMapper(periodicityMapper);
            periodicityActor->GetProperty()->SetColor(1.0, 1.0, 0.0); //(R,G,B)
            periodicityActor->SetVisibility(false);
            renderer->AddActor(periodicityActor);
            
            vtkSmartPointer<vtkPoints> periodicityPoints = vtkSmartPointer<vtkPoints>::New();
            vtkSmartPointer<vtkDoubleArray> periodicityVectors(vtkSmartPointer<vtkDoubleArray>::New());
            periodicityVectors->SetNumberOfComponents(3);
            for(const auto& pShift : mesh.periodicBasis())
            {
                periodicityPoints->InsertNextPoint(mesh.xMin().data());
                periodicityVectors->InsertNextTuple(pShift.data());
            }
            periodicityPolydata->SetPoints(periodicityPoints);
            periodicityPolydata->GetPointData()->SetVectors(periodicityVectors);
            periodicityPolydata->Modified();

            
//            axesWidget->SetOrientationMarker( axes );
////            axesWidget->SetOutlineColor( 0.9300, 0.5700, 0.1300 );
//            std::cout<<"Interactor"<<std::endl;
//            renderWindow->SetInteractor(renderWindow->MakeRenderWindowInteractor());
////            std::cout<(renderWindow->GetInteractor())<<std::endl;
//            std::cout<<bool(renderWindow->GetInteractor())<<std::endl;
//            axesWidget->SetInteractor( renderWindow->GetInteractor());
//            renderer->AddActor(axes);
//
////            axesWidget->SetViewport( 0.0, 0.0, 0.4, 0.4 );
//            axesWidget->SetEnabled( true );
//            axesWidget->InteractiveOn();


            mainLayout->addWidget(showMesh,0,0,1,1);
//            mainLayout->addWidget(showMeshLabel,0,1,1,1);
            mainLayout->addWidget(showFaceBoundaries,1,0,1,1);
//            mainLayout->addWidget(showFaceBoundariesLabel,1,1,1,1);
//            mainLayout->addWidget(showGrainColors,2,0,1,1);
  //          mainLayout->addWidget(showGrainColorsLabel,2,1,1,1);
            mainLayout->addWidget(showRegionBoundaries,2,0,1,1);
            mainLayout->addWidget(sliderRegionBoundaries,2,1,1,2);
            mainLayout->addWidget(showClipPlane,3,0,1,1);
            mainLayout->addWidget(sliderClipPlane,3,1,1,2);
            mainLayout->addWidget(showAxes,4,0,1,1);
            mainLayout->addWidget(showPeriodicityVectors,5,0,1,1);
            sliderRegionBoundaries->setOrientation(Qt::Horizontal);
            sliderRegionBoundaries->setMinimum(0);
            sliderRegionBoundaries->setMaximum(10);
            sliderRegionBoundaries->setValue(5);
            sliderClipPlane->setOrientation(Qt::Horizontal);
            sliderClipPlane->setMinimum(0);
            sliderClipPlane->setMaximum(10);
            sliderClipPlane->setValue(5);


   //         mainLayout->addWidget(showRegionBoundariesLabel,3,1,1,1);
            mainLayout->setColumnStretch(0, 1);
            mainLayout->setColumnStretch(1, 12);
            this->setLayout(mainLayout);

            connect(showMesh,SIGNAL(stateChanged(int)), this, SLOT(modify()));
            connect(showFaceBoundaries,SIGNAL(stateChanged(int)), this, SLOT(modify()));
//            connect(showGrainColors,SIGNAL(stateChanged(int)), this, SLOT(modify()));
            connect(showRegionBoundaries,SIGNAL(stateChanged(int)), this, SLOT(modify()));
            connect(showClipPlane,SIGNAL(stateChanged(int)), this, SLOT(modify()));
            connect(sliderRegionBoundaries,SIGNAL(valueChanged(int)), this, SLOT(modify()));
            connect(sliderClipPlane,SIGNAL(valueChanged(int)), this, SLOT(modify()));
            connect(showAxes,SIGNAL(stateChanged(int)), this, SLOT(modify()));
            connect(showPeriodicityVectors,SIGNAL(stateChanged(int)), this, SLOT(modify()));

            
             
            
            

            const Eigen::Matrix<double,3,1> c(0.5*(mesh.xMax()+mesh.xMin()));
            clipPlane->SetOrigin(c(0),c(1),c(2));
            clipPlane->SetNormal(1.0, 0.0, 0.0);

            
            polydata->Allocate();
            
            size_t pointID(0);
            for (const auto& node : mesh.observer<0>())
            {
                if(node.second->isBoundarySimplex())
                {
                    pts->InsertNextPoint(node.second->P0(0),
                                         node.second->P0(1),
                                         node.second->P0(2));
                    
                    //                    pts->InsertNextPoint(edge.second->child(1).P0(0),
                    //                                         edge.second->child(1).P0(1),
                    //                                         edge.second->child(1).P0(2));
                    sIDtoVtkPointsMap.emplace(node.second->xID(0),std::make_pair(pointID,Eigen::Matrix<double,3,1>::Zero()));
                    pointID++;
                }
            }
            
            //            std::cout<<"done inserting points"<<std::endl;
            
            
            
            //            size_t connectivityID=0;
            for (const auto& edge : mesh.observer<1>())
            {
                if(edge.second->isBoundarySimplex())
                {
                    const auto iter0(sIDtoVtkPointsMap.find(edge.second->child(0).xID(0)));
                    if(iter0!=sIDtoVtkPointsMap.end())
                    {
                        const auto iter1(sIDtoVtkPointsMap.find(edge.second->child(1).xID(0)));
                        if(iter1!=sIDtoVtkPointsMap.end())
                        {
                            vtkIdType connectivity[2];
                            connectivity[0] = iter0->second.first;
                            connectivity[1] = iter1->second.first;
                            polydata->InsertNextCell(VTK_LINE,2,connectivity); //Connects the first and fourth point we inserted into a line

                        }
                        else
                        {
                            throw std::runtime_error("child1 not found in sIDtoVtkPointsMap");
                        }
                    }
                    else
                    {
                        throw std::runtime_error("child0 not found in sIDtoVtkPointsMap");
                    }
                  }
            }
            
            
            
            
            polydata->SetPoints(pts);
            
            mapper->SetInputData(polydata);
            
            actor->SetMapper ( mapper );
            //            actor->GetProperty()->SetLineWidth(0.1);
            actor->GetProperty()->SetLineWidth(0.5);
            
            actor->GetProperty()->SetColor(0.5,0.5,0.5); // Give some color to the mesh. (1,1,1) is white
            //            actor->GetProperty()->SetColor(0.0,0.0,0.0); // Give some color to the mesh. (1,1,1) is white
            
            //            actor->GetProperty()->SetOpacity(0.15); //Make the mesh have some transparency.
            actor->GetProperty()->SetOpacity(0.5); //Make the mesh have some transparency.
            
            facePolydata->Allocate();
            size_t faceConnectivityID=0;
            for(auto region : mesh.regions())
            {// Sum number of external faces for final check
                for(auto& face : region.second->faces())
                {
                    for(size_t k=0;k<face.second->convexHull().size();++k)
                    {
                        //                        const int k1(k<convexHull().size()-1? k+1 : 0);
                        
                        facePts->InsertNextPoint(face.second->convexHull()[k]->P0(0),
                                                 face.second->convexHull()[k]->P0(1),
                                                 face.second->convexHull()[k]->P0(2));
                        
                        vtkIdType connectivity[2];
                        connectivity[0] = faceConnectivityID;
                        connectivity[1] = k<face.second->convexHull().size()-1? faceConnectivityID+1 : faceConnectivityID+1-face.second->convexHull().size();//connectivityID+1;
                        
                        facePolydata->InsertNextCell(VTK_LINE,2,connectivity); //Connects the first and fourth point we inserted into a line
                        
                        faceConnectivityID+=1;
                    }
                }
            }
            facePolydata->SetPoints(facePts);
            faceMapper->SetInputData(facePolydata);
            faceActor->SetMapper ( faceMapper );
            faceActor->GetProperty()->SetLineWidth(5.0);
            faceActor->GetProperty()->SetColor(0.0,0.0,0.0); // Give some color to the mesh. (1,1,1) is white
            faceActor->GetProperty()->SetOpacity(0.5); //Make the mesh have some transparency.
            
            
            
            gbColors->SetNumberOfComponents(3);
            gbColors->SetName("Colors");
            
            std::mt19937 generator;
            std::uniform_int_distribution<> distribution(0,255);
            std::map<int,Eigen::Matrix<int,1,3>> grainColors;
            for(const auto& region : mesh.regions())
            {
                grainColors.emplace(region.first,Eigen::Matrix<int,1,3>::Random()*255);
                //                grainColors.emplace(region.first,Eigen::Matrix<int,1,3>::Constant(clr)*255);
            }
            
            
//            if(showGrainColors)
//            {
                const double faceOffset(1.0e-1);
                size_t triPtID=0;

                for(const auto& meshTriangle : mesh.observer<2>())
                {// Color triangles the grain they belong to. TO DO: triangles can be replaced by mesh faces
                    if(meshTriangle.second->isBoundarySimplex())
                    {
                        const Eigen::Matrix<double,3,3> vM(meshTriangle.second->vertexPositionMatrix());
                        const int regionID=meshTriangle.second->parents().begin()->second->region->regionID;
                        const Eigen::Matrix<int,1,3>& regionClr(grainColors.at(regionID));
                        
                        gbPoints->InsertNextPoint ( vM(0,0), vM(1,0), vM(2,0) );
                        gbPoints->InsertNextPoint ( vM(0,1), vM(1,1), vM(2,1) );
                        gbPoints->InsertNextPoint ( vM(0,2), vM(1,2), vM(2,2) );
                        vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
                        triangle->GetPointIds()->SetId (0,triPtID+0);
                        triangle->GetPointIds()->SetId (1,triPtID+1);
                        triangle->GetPointIds()->SetId (2,triPtID+2);
                        gbColors->InsertNextTuple3(regionClr(0),regionClr(1),regionClr(2)); // use this to assig color to each vertex
                        gbTriangles->InsertNextCell ( triangle );
                        triPtID+=3;
                    }
                    else if(meshTriangle.second->isRegionBoundarySimplex())
                    {
                        const Eigen::Matrix<double,3,3> vM(meshTriangle.second->vertexPositionMatrix());
                        
                        const int regionID1=meshTriangle.second->parents().begin()->second->region->regionID;
                        const int regionID2=meshTriangle.second->parents().rbegin()->second->region->regionID;
                        
                        const Eigen::Matrix<double,3,1> n1=-meshTriangle.second->outNormal(regionID1).normalized();
                        const Eigen::Matrix<double,3,1> n2=-meshTriangle.second->outNormal(regionID2).normalized();
                        
                        const Eigen::Matrix<int,1,3>& regionClr1(grainColors.at(regionID1));
                        const Eigen::Matrix<int,1,3>& regionClr2(grainColors.at(regionID2));
                        
                        // Triangle on first grain
                        gbPoints->InsertNextPoint ( vM(0,0)+faceOffset*n1(0), vM(1,0)+faceOffset*n1(1), vM(2,0)+faceOffset*n1(2) );
                        gbPoints->InsertNextPoint ( vM(0,1)+faceOffset*n1(0), vM(1,1)+faceOffset*n1(1), vM(2,1)+faceOffset*n1(2) );
                        gbPoints->InsertNextPoint ( vM(0,2)+faceOffset*n1(0), vM(1,2)+faceOffset*n1(1), vM(2,2)+faceOffset*n1(2) );
                        vtkSmartPointer<vtkTriangle> triangle1 = vtkSmartPointer<vtkTriangle>::New();
                        triangle1->GetPointIds()->SetId (0,triPtID+0);
                        triangle1->GetPointIds()->SetId (1,triPtID+1);
                        triangle1->GetPointIds()->SetId (2,triPtID+2);
                        gbTriangles->InsertNextCell ( triangle1 );
                        gbColors->InsertNextTuple3(regionClr1(0),regionClr1(1),regionClr1(2)); // use this to assig color to each vertex
                        triPtID+=3;
                        
                        // Triangle on second grain
                        gbPoints->InsertNextPoint ( vM(0,0)+faceOffset*n2(0), vM(1,0)+faceOffset*n2(1), vM(2,0)+faceOffset*n2(2) );
                        gbPoints->InsertNextPoint ( vM(0,1)+faceOffset*n2(0), vM(1,1)+faceOffset*n2(1), vM(2,1)+faceOffset*n2(2) );
                        gbPoints->InsertNextPoint ( vM(0,2)+faceOffset*n2(0), vM(1,2)+faceOffset*n2(1), vM(2,2)+faceOffset*n2(2) );
                        vtkSmartPointer<vtkTriangle> triangle2 = vtkSmartPointer<vtkTriangle>::New();
                        triangle2->GetPointIds()->SetId (0,triPtID+0);
                        triangle2->GetPointIds()->SetId (1,triPtID+1);
                        triangle2->GetPointIds()->SetId (2,triPtID+2);
                        gbTriangles->InsertNextCell ( triangle2 );
                        gbColors->InsertNextTuple3(regionClr2(0),regionClr2(1),regionClr2(2)); // use this to assig color to each vertex
                        triPtID+=3;
                    }
                }
//            }
            
            gbTrianglePolyData->SetPoints ( gbPoints );
            gbTrianglePolyData->SetPolys ( gbTriangles );
            gbTrianglePolyData->GetCellData()->SetScalars(gbColors);
            gbMapper->SetInputData(gbTrianglePolyData);
            gbActor->SetMapper(gbMapper);
            gbActor->GetProperty()->SetOpacity(0.2);
            gbActor->GetProperty()->SetColor(0.0,0.5,0.5);
            
            clipper->SetInputData(gbTrianglePolyData);
            clipper->SetClipFunction(clipPlane);
            clipper->SetValue(0);
            clipper->Update();
            clipMapper->SetInputConnection(clipper->GetOutputPort());
            clipActor->SetMapper(clipMapper);
            
            vtkNew<vtkNamedColors> colors;
            vtkColor3d axisColor=colors->GetColor3d("LightSkyBlue");
            const int fontSize=12;
            
            
//            cubeAxesActor->GetTextActor()->SetTextScaleModeToNone();
//              cubeAxesActor->SetUseTextActor3D(1);
//            cubeAxesActor->SetUse2DMode(true);
                //cubeAxesActor->SetUseTextActor3D(false);
              cubeAxesActor->SetBounds(polydata->GetBounds());
              cubeAxesActor->SetCamera(renderer->GetActiveCamera());
//            cubeAxesActor->SetTitleOffset(20+fontSize);
//            cubeAxesActor->GetAxisLabelTextProperty()->SetFontSize(fontsize);
//            cubeAxesActor->GetAxisTitleTextProperty()->SetFontSize(fontsize);
            
              cubeAxesActor->GetTitleTextProperty(0)->SetColor(axisColor.GetData());
              cubeAxesActor->GetTitleTextProperty(0)->SetFontSize(fontSize);
              cubeAxesActor->GetLabelTextProperty(0)->SetColor(axisColor.GetData());
            cubeAxesActor->GetLabelTextProperty(0)->SetFontSize(fontSize);
            cubeAxesActor->XAxisTickVisibilityOn();
            cubeAxesActor->GetXAxesLinesProperty()->SetColor(axisColor.GetData());
            cubeAxesActor->GetXAxesLinesProperty()->SetLineWidth(4);
            cubeAxesActor->DrawXGridlinesOn();
            cubeAxesActor->GetXAxesGridlinesProperty()->SetColor(axisColor.GetData());
            cubeAxesActor->SetXTitle("X");

              cubeAxesActor->GetTitleTextProperty(1)->SetColor(axisColor.GetData());
              cubeAxesActor->GetTitleTextProperty(1)->SetFontSize(fontSize);
              cubeAxesActor->GetLabelTextProperty(1)->SetColor(axisColor.GetData());
            cubeAxesActor->GetLabelTextProperty(1)->SetFontSize(fontSize);
            cubeAxesActor->YAxisTickVisibilityOn();
            cubeAxesActor->GetYAxesLinesProperty()->SetColor(axisColor.GetData());
            cubeAxesActor->GetYAxesLinesProperty()->SetLineWidth(4);
            cubeAxesActor->DrawYGridlinesOn();
            cubeAxesActor->GetYAxesGridlinesProperty()->SetColor(axisColor.GetData());
            cubeAxesActor->SetYTitle("Y");

            
              cubeAxesActor->GetTitleTextProperty(2)->SetColor(axisColor.GetData());
              cubeAxesActor->GetTitleTextProperty(2)->SetFontSize(fontSize);
              cubeAxesActor->GetLabelTextProperty(2)->SetColor(axisColor.GetData());
            cubeAxesActor->GetLabelTextProperty(2)->SetFontSize(fontSize);
            cubeAxesActor->ZAxisTickVisibilityOn();
            cubeAxesActor->GetZAxesLinesProperty()->SetColor(axisColor.GetData());
            cubeAxesActor->GetZAxesLinesProperty()->SetLineWidth(4);
            cubeAxesActor->DrawZGridlinesOn();
            cubeAxesActor->GetZAxesGridlinesProperty()->SetColor(axisColor.GetData());
            cubeAxesActor->SetZTitle("Z");
            
            
            //            cubeAxesActor->SetGridLineLocation(cubeAxesActor->VTK_GRID_LINES_FURTHEST);

//            cubeAxesActor->XAxisMinorTickVisibilityOff();
//              cubeAxesActor->YAxisMinorTickVisibilityOff();
//              cubeAxesActor->ZAxisMinorTickVisibilityOff();
//            cubeAxesActor->SetColor(0.0,0.0,0.0,0.0);

              cubeAxesActor->SetFlyModeToStaticTriad();
            
            
            renderer->AddActor(gbActor);
            renderer->AddActor(clipActor); // enable to clip grain boundaries
            renderer->AddActor(cubeAxesActor);

            
            // Update
            //            update(0/*,0,0.0,Eigen::Matrix<float,3,1>::Zero()*/);
            modifyPts();

            // Render
            renderer->AddActor(actor);
            renderer->AddActor(faceActor);
            renderer->ResetCamera();


            
//            myCallback->plane=clipPlane;
//            rep->SetPlaceFactor(1.0); // This must be set prior to placing the widget
//            rep->PlaceWidget(clipActor->GetBounds());
//            rep->SetNormal(clipPlane->GetNormal());
//            //            planeWidget->SetInteractor(renderWindow->GetInteractor());
//            planeWidget->AddObserver(vtkCommand::InteractionEvent,myCallback);
//            planeWidget->SetRepresentation(rep);
//            planeWidget->setEnabled(true );

//            modify();
        }
        

        
 
        
        void SimplicialMeshActor::modify()
        {
            actor->SetVisibility(showMesh->isChecked());
            faceActor->SetVisibility(showFaceBoundaries->isChecked());
            gbActor->SetVisibility(showRegionBoundaries->isChecked());
            sliderRegionBoundaries->setEnabled(showRegionBoundaries->isChecked());
            gbActor->GetProperty()->SetOpacity(sliderRegionBoundaries->value()/10.0);
            clipActor->SetVisibility(showClipPlane->isChecked());
            sliderClipPlane->setEnabled(showClipPlane->isChecked());
            clipActor->GetProperty()->SetOpacity(sliderClipPlane->value()/10.0);

            cubeAxesActor->SetVisibility(showAxes->isChecked());

            periodicityActor->SetVisibility(showPeriodicityVectors->isChecked());
            
            renderWindow->Render();
        }


        /**************************************************************************/
        void SimplicialMeshActor::modifyPts()
        {
            for (const auto& node : mesh.observer<0>())
            {
                if(node.second->isBoundarySimplex())
                {
                    
                    const auto nodeIter(sIDtoVtkPointsMap.find(node.second->xID(0)));
                    if(nodeIter!=sIDtoVtkPointsMap.end())
                    {
                        Eigen::Matrix<double,3,1> newP(node.second->P0+dispCorr*nodeIter->second.second);
                        pts->SetPoint(nodeIter->second.first,newP(0),newP(1),newP(2));

                    }
                    else
                    {
                        throw std::runtime_error("node not found in sIDtoVtkPointsMap");
                    }
                    
                }
            }
            pts->Modified();
        }

} // namespace model

 // namespace model
#endif







