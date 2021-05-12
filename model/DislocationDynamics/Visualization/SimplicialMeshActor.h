/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SimplicialMeshActor_H_
#define model_SimplicialMeshActor_H_

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
#include <random>
#include <algorithm>

#include <TextFileParser.h>

#include <SimplicialMesh.h>
#include <DDauxIO.h>

//#include <VertexReader.h>

// See this for plane interactor
// https://www.vtk.org/Wiki/VTK/Examples/Cxx/Widgets/ImplicitPlaneWidget2

// camera
// https://www.vtk.org/pipermail/vtkusers/2014-January/082864.html

namespace model
{
    

    
    
    //    struct SimplicialMeshActor : public VertexReader<'D',4,float>//: public vtkSmartPointer<vtkActor>
    struct SimplicialMeshActor  //public IDreader<'D',1,3,float>//: public vtkSmartPointer<vtkActor>
    {
        
        //        typedef VertexReader<'D',4,float> DispContainerType;
//        typedef IDreader<'D',1,3,float> DispContainerType;
        
        static double dispCorr;
        
        std::map<size_t,std::pair<size_t,Eigen::Matrix<double,3,1>>> sIDtoVtkPointsMap;
        
        vtkSmartPointer<vtkPoints> pts;
        vtkSmartPointer<vtkPolyData> polydata;
        vtkSmartPointer<vtkPolyDataMapper> mapper;
        vtkSmartPointer<vtkActor> actor;

        vtkSmartPointer<vtkPoints> facePts;
        vtkSmartPointer<vtkPolyData> facePolydata;
        vtkSmartPointer<vtkPolyDataMapper> faceMapper;
        vtkSmartPointer<vtkActor> faceActor;

        
        vtkSmartPointer<vtkUnsignedCharArray> gbColors;
        vtkSmartPointer<vtkPoints> gbPoints;
        vtkSmartPointer<vtkCellArray> gbTriangles;
        vtkSmartPointer<vtkPolyData> gbTrianglePolyData;
        vtkSmartPointer<vtkPolyDataMapper> gbMapper;
        vtkSmartPointer<vtkActor> gbActor;
        
        vtkSmartPointer<vtkPlane> clipPlane;
        vtkSmartPointer<vtkClipPolyData> clipper;
        vtkSmartPointer<vtkPolyData> clippedPolyData;
        vtkSmartPointer<vtkDataSetMapper> clipMapper;
        vtkSmartPointer<vtkActor> clipActor;
        
//        vtkSmartPointer<vtkSphereSource> fieldSphereSource;
//        vtkSmartPointer<vtkPlane> fieldPlane;
//        vtkSmartPointer<vtkClipPolyData> fieldClipper;
//        vtkSmartPointer<vtkPolyDataMapper> fieldMapper;
//        vtkSmartPointer<vtkActor> fieldActor;
//        vtkSmartPointer<vtkIPWCallback> fieldCallback;
//        vtkSmartPointer<vtkImplicitPlaneRepresentation> fieldRep;
//        vtkSmartPointer<vtkImplicitPlaneWidget2> fieldPlaneWidget;
        
        vtkRenderer* const renderer;
        const SimplicialMesh<3>& mesh;
        static bool showGrainColors;
        static bool showRegionBoundaries;
        static bool showFaceBoundaries;
        static bool showMesh;

        /**************************************************************************/
        SimplicialMeshActor(vtkRenderer* const renderer_in, const SimplicialMesh<3>& mesh_in) :
        /* init */ renderer(renderer_in)
        /* init */,mesh(mesh_in)
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
        {
            
            Eigen::Matrix<double,3,1> c(0.5*(mesh.xMax()+mesh.xMin()));
            clipPlane->SetOrigin(c(0),c(1),c(2));
            clipPlane->SetNormal(1.0, 0.0, 0.0);
            
//            mesh.readMesh(TextFileParser("./inputFiles/polycrystal.txt").readString("meshFile",true));
            
            
            
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
                    assert(iter0!=sIDtoVtkPointsMap.end() && "child0 not found in sIDtoVtkPointsMap");
                    
                    const auto iter1(sIDtoVtkPointsMap.find(edge.second->child(1).xID(0)));
                    assert(iter1!=sIDtoVtkPointsMap.end() && "child0 not found in sIDtoVtkPointsMap");
                    
                    vtkIdType connectivity[2];
                    connectivity[0] = iter0->second.first;
                    connectivity[1] = iter1->second.first;
                    polydata->InsertNextCell(VTK_LINE,2,connectivity); //Connects the first and fourth point we inserted into a line
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
                    for(int k=0;k<face.second->convexHull().size();++k)
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
            faceActor->GetProperty()->SetLineWidth(2.0);
            faceActor->GetProperty()->SetColor(0.0,0.0,0.0); // Give some color to the mesh. (1,1,1) is white
            faceActor->GetProperty()->SetOpacity(0.5); //Make the mesh have some transparency.
            
            
            
            gbColors->SetNumberOfComponents(3);
            gbColors->SetName("Colors");
            size_t triPtID=0;
            
            std::mt19937 generator;
            std::uniform_int_distribution<> distribution(0,255);
            
            
            std::map<int,Eigen::Matrix<int,1,3>> grainColors;
            for(const auto& region : mesh.regions())
            {
                const int clr=distribution(generator); // a random SlipSystem ID
                grainColors.emplace(region.first,Eigen::Matrix<int,1,3>::Random()*255);
                //                grainColors.emplace(region.first,Eigen::Matrix<int,1,3>::Constant(clr)*255);
            }
            
            
            if(showGrainColors)
            {// Color triangles the grain they belong to. TO DO: triangles can be replaced by mesh faces
                const double faceOffset(1.0e-1);
                
                for(const auto& meshTriangle : mesh.observer<2>())
                {
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
            }
            
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
            
            
            renderer->AddActor(gbActor);
            //            renderer->AddActor(clipActor); // enable to clip grain boundaries
            
            
            // Update
            //            update(0/*,0,0.0,Eigen::Matrix<float,3,1>::Zero()*/);
            modifyPts();
            
            // Render
            renderer->AddActor(actor);
            renderer->AddActor(faceActor);
            
        }
        
//        /**************************************************************************/
//        void init(vtkRenderer* renderer)
//        {
//
//            mesh.readMesh(TextFileParser("./inputFiles/polycrystal.txt").readString("meshFile",true));
//
//
//
//            polydata->Allocate();
//
//            size_t pointID(0);
//            for (const auto& node : mesh.observer<0>())
//            {
//                if(node.second->isBoundarySimplex())
//                {
//                    pts->InsertNextPoint(node.second->P0(0),
//                                         node.second->P0(1),
//                                         node.second->P0(2));
//
////                    pts->InsertNextPoint(edge.second->child(1).P0(0),
////                                         edge.second->child(1).P0(1),
////                                         edge.second->child(1).P0(2));
//                    sIDtoVtkPointsMap.emplace(node.second->xID(0),std::make_pair(pointID,Eigen::Matrix<double,3,1>::Zero()));
//                    pointID++;
//                }
//            }
//
////            std::cout<<"done inserting points"<<std::endl;
//
//
//
////            size_t connectivityID=0;
//            for (const auto& edge : mesh.observer<1>())
//            {
//                if(edge.second->isBoundarySimplex())
//                {
//                    const auto iter0(sIDtoVtkPointsMap.find(edge.second->child(0).xID(0)));
//                    assert(iter0!=sIDtoVtkPointsMap.end() && "child0 not found in sIDtoVtkPointsMap");
//
//                    const auto iter1(sIDtoVtkPointsMap.find(edge.second->child(1).xID(0)));
//                    assert(iter1!=sIDtoVtkPointsMap.end() && "child0 not found in sIDtoVtkPointsMap");
//
//                    vtkIdType connectivity[2];
//                    connectivity[0] = iter0->second.first;
//                    connectivity[1] = iter1->second.first;
//                    polydata->InsertNextCell(VTK_LINE,2,connectivity); //Connects the first and fourth point we inserted into a line
//                }
//            }
//
//
//
//
//            polydata->SetPoints(pts);
//
//            mapper->SetInputData(polydata);
//
//            actor->SetMapper ( mapper );
////            actor->GetProperty()->SetLineWidth(0.1);
//            actor->GetProperty()->SetLineWidth(0.5);
//
//            actor->GetProperty()->SetColor(0.5,0.5,0.5); // Give some color to the mesh. (1,1,1) is white
////            actor->GetProperty()->SetColor(0.0,0.0,0.0); // Give some color to the mesh. (1,1,1) is white
//
//            //            actor->GetProperty()->SetOpacity(0.15); //Make the mesh have some transparency.
//            actor->GetProperty()->SetOpacity(0.5); //Make the mesh have some transparency.
//
//            facePolydata->Allocate();
//            size_t faceConnectivityID=0;
//            for(auto region : mesh.regions())
//            {// Sum number of external faces for final check
//                for(auto& face : region.second->faces())
//                {
//                    for(int k=0;k<face.second->convexHull().size();++k)
//                    {
////                        const int k1(k<convexHull().size()-1? k+1 : 0);
//
//                        facePts->InsertNextPoint(face.second->convexHull()[k]->P0(0),
//                                                 face.second->convexHull()[k]->P0(1),
//                                                 face.second->convexHull()[k]->P0(2));
//
//                        vtkIdType connectivity[2];
//                        connectivity[0] = faceConnectivityID;
//                        connectivity[1] = k<face.second->convexHull().size()-1? faceConnectivityID+1 : faceConnectivityID+1-face.second->convexHull().size();//connectivityID+1;
//
//                        facePolydata->InsertNextCell(VTK_LINE,2,connectivity); //Connects the first and fourth point we inserted into a line
//
//                        faceConnectivityID+=1;
//                    }
//                }
//            }
//            facePolydata->SetPoints(facePts);
//            faceMapper->SetInputData(facePolydata);
//            faceActor->SetMapper ( faceMapper );
//            faceActor->GetProperty()->SetLineWidth(2.0);
//            faceActor->GetProperty()->SetColor(0.0,0.0,0.0); // Give some color to the mesh. (1,1,1) is white
//            faceActor->GetProperty()->SetOpacity(0.5); //Make the mesh have some transparency.
//
//
//
//            gbColors->SetNumberOfComponents(3);
//            gbColors->SetName("Colors");
//            size_t triPtID=0;
//
//            std::mt19937 generator;
//            std::uniform_int_distribution<> distribution(0,255);
//
//
//            std::map<int,Eigen::Matrix<int,1,3>> grainColors;
//            for(const auto& region : mesh.regions())
//            {
//                const int clr=distribution(generator); // a random SlipSystem ID
//                grainColors.emplace(region.first,Eigen::Matrix<int,1,3>::Random()*255);
////                grainColors.emplace(region.first,Eigen::Matrix<int,1,3>::Constant(clr)*255);
//            }
//
//
//            if(showGrainColors)
//            {// Color triangles the grain they belong to. TO DO: triangles can be replaced by mesh faces
//                const double faceOffset(1.0e-1);
//
//                for(const auto& meshTriangle : mesh.observer<2>())
//                {
//                    if(meshTriangle.second->isBoundarySimplex())
//                    {
//                        const Eigen::Matrix<double,3,3> vM(meshTriangle.second->vertexPositionMatrix());
//                        const int regionID=meshTriangle.second->parents().begin()->second->region->regionID;
//                        const Eigen::Matrix<int,1,3>& regionClr(grainColors.at(regionID));
//
//                        gbPoints->InsertNextPoint ( vM(0,0), vM(1,0), vM(2,0) );
//                        gbPoints->InsertNextPoint ( vM(0,1), vM(1,1), vM(2,1) );
//                        gbPoints->InsertNextPoint ( vM(0,2), vM(1,2), vM(2,2) );
//                        vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
//                        triangle->GetPointIds()->SetId (0,triPtID+0);
//                        triangle->GetPointIds()->SetId (1,triPtID+1);
//                        triangle->GetPointIds()->SetId (2,triPtID+2);
//                        gbColors->InsertNextTuple3(regionClr(0),regionClr(1),regionClr(2)); // use this to assig color to each vertex
//                        gbTriangles->InsertNextCell ( triangle );
//                        triPtID+=3;
//                    }
//                    else if(meshTriangle.second->isRegionBoundarySimplex())
//                    {
//                        const Eigen::Matrix<double,3,3> vM(meshTriangle.second->vertexPositionMatrix());
//
//                        const int regionID1=meshTriangle.second->parents().begin()->second->region->regionID;
//                        const int regionID2=meshTriangle.second->parents().rbegin()->second->region->regionID;
//
//                        const Eigen::Matrix<double,3,1> n1=-meshTriangle.second->outNormal(regionID1).normalized();
//                        const Eigen::Matrix<double,3,1> n2=-meshTriangle.second->outNormal(regionID2).normalized();
//
//                        const Eigen::Matrix<int,1,3>& regionClr1(grainColors.at(regionID1));
//                        const Eigen::Matrix<int,1,3>& regionClr2(grainColors.at(regionID2));
//
//                        // Triangle on first grain
//                        gbPoints->InsertNextPoint ( vM(0,0)+faceOffset*n1(0), vM(1,0)+faceOffset*n1(1), vM(2,0)+faceOffset*n1(2) );
//                        gbPoints->InsertNextPoint ( vM(0,1)+faceOffset*n1(0), vM(1,1)+faceOffset*n1(1), vM(2,1)+faceOffset*n1(2) );
//                        gbPoints->InsertNextPoint ( vM(0,2)+faceOffset*n1(0), vM(1,2)+faceOffset*n1(1), vM(2,2)+faceOffset*n1(2) );
//                        vtkSmartPointer<vtkTriangle> triangle1 = vtkSmartPointer<vtkTriangle>::New();
//                        triangle1->GetPointIds()->SetId (0,triPtID+0);
//                        triangle1->GetPointIds()->SetId (1,triPtID+1);
//                        triangle1->GetPointIds()->SetId (2,triPtID+2);
//                        gbTriangles->InsertNextCell ( triangle1 );
//                        gbColors->InsertNextTuple3(regionClr1(0),regionClr1(1),regionClr1(2)); // use this to assig color to each vertex
//                        triPtID+=3;
//
//                        // Triangle on second grain
//                        gbPoints->InsertNextPoint ( vM(0,0)+faceOffset*n2(0), vM(1,0)+faceOffset*n2(1), vM(2,0)+faceOffset*n2(2) );
//                        gbPoints->InsertNextPoint ( vM(0,1)+faceOffset*n2(0), vM(1,1)+faceOffset*n2(1), vM(2,1)+faceOffset*n2(2) );
//                        gbPoints->InsertNextPoint ( vM(0,2)+faceOffset*n2(0), vM(1,2)+faceOffset*n2(1), vM(2,2)+faceOffset*n2(2) );
//                        vtkSmartPointer<vtkTriangle> triangle2 = vtkSmartPointer<vtkTriangle>::New();
//                        triangle2->GetPointIds()->SetId (0,triPtID+0);
//                        triangle2->GetPointIds()->SetId (1,triPtID+1);
//                        triangle2->GetPointIds()->SetId (2,triPtID+2);
//                        gbTriangles->InsertNextCell ( triangle2 );
//                        gbColors->InsertNextTuple3(regionClr2(0),regionClr2(1),regionClr2(2)); // use this to assig color to each vertex
//                        triPtID+=3;
//                    }
//                }
//            }
//
//            gbTrianglePolyData->SetPoints ( gbPoints );
//            gbTrianglePolyData->SetPolys ( gbTriangles );
//            gbTrianglePolyData->GetCellData()->SetScalars(gbColors);
//            gbMapper->SetInputData(gbTrianglePolyData);
//            gbActor->SetMapper(gbMapper);
//            gbActor->GetProperty()->SetOpacity(0.2);
//            gbActor->GetProperty()->SetColor(0.0,0.5,0.5);
//
//            clipper->SetInputData(gbTrianglePolyData);
//            clipper->SetClipFunction(clipPlane);
//            clipper->SetValue(0);
//            clipper->Update();
//            clipMapper->SetInputConnection(clipper->GetOutputPort());
//            clipActor->SetMapper(clipMapper);
//
//
//            renderer->AddActor(gbActor);
////            renderer->AddActor(clipActor); // enable to clip grain boundaries
//
//
//            // Update
////            update(0/*,0,0.0,Eigen::Matrix<float,3,1>::Zero()*/);
//            modifyPts();
//
//            // Render
//            renderer->AddActor(actor);
//            renderer->AddActor(faceActor);
////            renderer->AddActor(fieldActor);
//
//        }
        
        /**************************************************************************/
        void update(const DDauxIO<3>& ddAux)
        {
            
            
            

            for(const auto& node : ddAux.meshNodes())
            {
                auto iter1(sIDtoVtkPointsMap.find(node.nodeID));
                if(iter1!=sIDtoVtkPointsMap.end())
                {
                    iter1->second.second=node.displacement;
                }
//                assert(iter1!=sIDtoVtkPointsMap.end() && "child0 not found in sIDtoVtkPointsMap");
            }
            

            
            modifyPts();

//            if(DispContainerType::isGood(frameID,true))
//            {
//                dispFileIsGood=true;
//                DispContainerType::read(frameID,true);
//
//                modifyPts();
//
//
//            }
//            else
//            {
//                dispFileIsGood=false;
//            }
            
//            const double axisNorm(spinAxis.norm());
//            if(anglePerStep && axisNorm>0.0)
//            {
//                spinAxis/=axisNorm;
////                std::cout<<"rotating"<< frameID*anglePerStep<<std::endl;
//                const double splinAngle((frameID-lastFrameID)*anglePerStep);
//                std::cout<<"mesh "<<splinAngle<<std::endl;
//
//                actor->RotateWXYZ(splinAngle,spinAxis(0),spinAxis(1),spinAxis(2));
//                gbActor->RotateWXYZ(splinAngle,spinAxis(0),spinAxis(1),spinAxis(2));
//                clipActor->RotateWXYZ(splinAngle,spinAxis(0),spinAxis(1),spinAxis(2));
//
//            }
            
            
        }
        
        /**************************************************************************/
        void modifyPts()
        {
//            if(dispFileIsGood)
//            {
//                size_t k=0;
//                for (const auto& edge : mesh.observer<1>())
//                {
//                    if(edge.second->isBoundarySimplex())
//                    {
//                        DispContainerType::const_iterator iterD1(DispContainerType::find(edge.second->child(0).xID(0)));
//                        DispContainerType::const_iterator iterD2(DispContainerType::find(edge.second->child(1).xID(0)));
//                        //                        VertexReader<'D',4,float>::const_iterator iterD3(DispContainerType::find(triangleID.second[2]));
//
//                        assert(iterD1!=DispContainerType::end() && "MESH NODE NOT FOUND IN D FILE");
//                        assert(iterD2!=DispContainerType::end() && "MESH NODE NOT FOUND IN D FILE");
//                        //                        assert(iterD3!=DispContainerType::end() && "MESH NODE NOT FOUND IN D FILE");
//
//
//                        const float x1=edge.second->child(0).P0(0)+dispCorr*iterD1->second[0];
//                        const float y1=edge.second->child(0).P0(1)+dispCorr*iterD1->second[1];
//                        const float z1=edge.second->child(0).P0(2)+dispCorr*iterD1->second[2];
//                        pts->SetPoint(k,x1,y1,z1);
//
//                        const float x2=edge.second->child(1).P0(0)+dispCorr*iterD2->second[0];
//                        const float y2=edge.second->child(1).P0(1)+dispCorr*iterD2->second[1];
//                        const float z2=edge.second->child(1).P0(2)+dispCorr*iterD2->second[2];
//                        pts->SetPoint(k+1,x2,y2,z2);
//
//                        k=k+2;
//                    }
//                }
//
//                pts->Modified();
//            }
            
//            for(const auto& pair : sIDtoVtkPointsMap)
//            {
//                const auto nodeIter(mesh.observer<1>())
//                pts->SetPoint(pair.second.first,x2,y2,z2)
//            }
            
            for (const auto& node : mesh.observer<0>())
            {
                if(node.second->isBoundarySimplex())
                {
                    
                    const auto nodeIter(sIDtoVtkPointsMap.find(node.second->xID(0)));
                    assert(nodeIter!=sIDtoVtkPointsMap.end() && "node not found in sIDtoVtkPointsMap");
                    
                    Eigen::Matrix<double,3,1> newP(node.second->P0+dispCorr*nodeIter->second.second);
                    pts->SetPoint(nodeIter->second.first,newP(0),newP(1),newP(2));
                    
                }
            }
            pts->Modified();
            
//            size_t t=0;
//            for(const auto& meshTriangle : mesh.observer<2>())
//            {
//                if(meshTriangle.second->isBoundarySimplex())
//                {
//                    gbColors->SetTuple3(t,regionClr(1),regionClr(2)); // use this to assig color to each vertex
//                    t++;
//                }
//            }
            
            gbMapper->SetScalarVisibility(showGrainColors);
            
            if(showRegionBoundaries)
            {
                gbActor->VisibilityOn();
            }
            else
            {
                gbActor->VisibilityOff();
            }
//            gbActor->GetProperty()->SetColor(0.0,0.5,0.5);
            
            if(showFaceBoundaries)
            {
                faceActor->VisibilityOn();
            }
            else
            {
                faceActor->VisibilityOff();

            }
            
            if(showMesh)
            {
                actor->VisibilityOn();
            }
            else
            {
                actor->VisibilityOff();
                
            }
            
        }
        
    };
    
    double SimplicialMeshActor::dispCorr=1.0;
    bool SimplicialMeshActor::showGrainColors=true;
    bool SimplicialMeshActor::showRegionBoundaries=false;
    bool SimplicialMeshActor::showFaceBoundaries=true;
    bool SimplicialMeshActor::showMesh=true;

} // namespace model
#endif







