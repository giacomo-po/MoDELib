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
#include <random>


#include <model/Mesh/SimplicialMesh.h>
//#include <model/IO/VertexReader.h>
#include <model/IO/IDreader.h>

// See this for plane interactor
// https://www.vtk.org/Wiki/VTK/Examples/Cxx/Widgets/ImplicitPlaneWidget2

// camera
// https://www.vtk.org/pipermail/vtkusers/2014-January/082864.html

namespace model
{
    
    //    struct SimplicialMeshActor : public VertexReader<'D',4,float>//: public vtkSmartPointer<vtkActor>
    struct SimplicialMeshActor : public IDreader<'D',1,3,float>//: public vtkSmartPointer<vtkActor>
    {
        
        //        typedef VertexReader<'D',4,float> DispContainerType;
        typedef IDreader<'D',1,3,float> DispContainerType;
        
        static double dispCorr;
        
        vtkSmartPointer<vtkPoints> pts;
        vtkSmartPointer<vtkPolyData> polydata;
        vtkSmartPointer<vtkPolyDataMapper> mapper;
        vtkSmartPointer<vtkActor> actor;
        
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
        
        SimplicialMesh<3> mesh;
        bool dispFileIsGood;
        static bool showGrainColors;
        static bool showRegionBoundaries;
        
        
        /**************************************************************************/
        SimplicialMeshActor() :
        /* init */ pts(vtkSmartPointer<vtkPoints>::New()),
        /* init */ polydata(vtkSmartPointer<vtkPolyData>::New()),
        /* init */ mapper(vtkSmartPointer<vtkPolyDataMapper>::New()),
        /* init */ actor(vtkSmartPointer<vtkActor>::New()),
        gbColors(vtkSmartPointer<vtkUnsignedCharArray>::New()),
        gbPoints(vtkSmartPointer<vtkPoints>::New()),
        gbTriangles(vtkSmartPointer<vtkCellArray>::New()),
        gbTrianglePolyData(vtkSmartPointer<vtkPolyData>::New()),
        gbMapper(vtkSmartPointer<vtkPolyDataMapper>::New()),
        gbActor(vtkSmartPointer<vtkActor>::New()),
        clipPlane(vtkSmartPointer<vtkPlane>::New()),
        clipper(vtkSmartPointer<vtkClipPolyData>::New()),
        clippedPolyData(vtkSmartPointer<vtkPolyData>::New()),
        clipMapper(vtkSmartPointer<vtkDataSetMapper>::New()),
        clipActor(vtkSmartPointer<vtkActor>::New()),
        /* init */ dispFileIsGood(false)
        {
            
            Eigen::Matrix<double,3,1> c(0.5*(mesh.xMax()+mesh.xMin()));
            clipPlane->SetOrigin(c(0),c(1),c(2));
            clipPlane->SetNormal(1.0, 0.0, 0.0);
            
        }
        
        /**************************************************************************/
        void init(const int& meshID,vtkRenderer* renderer)
        {
            mesh.readMesh(meshID);
            
            polydata->Allocate();
            
            size_t connectivityID=0;
            for (const auto& edge : mesh.observer<1>())
            {
                if(edge.second->isBoundarySimplex())
                {
                    pts->InsertNextPoint(edge.second->child(0).P0(0),
                                         edge.second->child(0).P0(1),
                                         edge.second->child(0).P0(2));
                    
                    pts->InsertNextPoint(edge.second->child(1).P0(0),
                                         edge.second->child(1).P0(1),
                                         edge.second->child(1).P0(2));
                    
                    vtkIdType connectivity[2];
                    connectivity[0] = connectivityID;
                    connectivity[1] = connectivityID+1;
                    polydata->InsertNextCell(VTK_LINE,2,connectivity); //Connects the first and fourth point we inserted into a line
                    
                    connectivityID+=2;
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
            {
                
                for(const auto& meshTriangle : mesh.observer<2>())
                {
                    if(meshTriangle.second->isBoundarySimplex())
                    {
                        const int regionID=(*meshTriangle.second->parents().begin())->region->regionID;
                        
                        const Eigen::Matrix<double,3,3> vM(meshTriangle.second->vertexPositionMatrix());
                        
                        gbPoints->InsertNextPoint ( vM(0,0), vM(1,0), vM(2,0) );
                        gbPoints->InsertNextPoint ( vM(0,1), vM(1,1), vM(2,1) );
                        gbPoints->InsertNextPoint ( vM(0,2), vM(1,2), vM(2,2) );
                        
                        
                        vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
                        triangle->GetPointIds()->SetId (0,triPtID+0);
                        triangle->GetPointIds()->SetId (1,triPtID+1);
                        triangle->GetPointIds()->SetId (2,triPtID+2);
                        
                        
                        const Eigen::Matrix<int,1,3>& regionClr(grainColors.at(regionID));
                        gbColors->InsertNextTuple3(regionClr(0),regionClr(1),regionClr(2)); // use this to assig color to each vertex
                        
                        
                        gbTriangles->InsertNextCell ( triangle );
                        
                        triPtID+=3;
                    }
                }
            }
            
            const double faceOffset(1.0e-1);
            
            // grain-boundaries
            if(showRegionBoundaries)
            {
                
                
                //                Eigen::MatrixXi regionColors(Eigen::MatrixXi::Random(mesh.regions().size(),3));
                
                unsigned char clr[3] = {0, 100, 100};
                //                gbColors->SetNumberOfComponents(3);
                //                gbColors->SetName("Colors");
                //                gbColors->InsertNextTypedTuple(red);
                //                size_t triPtID=0;
                for(const auto& meshTriangle : mesh.observer<2>())
                {
                    if(meshTriangle.second->isRegionBoundarySimplex())
                    {
                        const Eigen::Matrix<double,3,3> vM(meshTriangle.second->vertexPositionMatrix());
                        
                        const int regionID1=(*meshTriangle.second->parents().begin())->region->regionID;
                        const int regionID2=(*meshTriangle.second->parents().rbegin())->region->regionID;

                        const Eigen::Matrix<double,3,1> n1=-meshTriangle.second->outNormal(regionID1).normalized();
                        const Eigen::Matrix<double,3,1> n2=-meshTriangle.second->outNormal(regionID2).normalized();
                        
                        const Eigen::Matrix<int,1,3>& regionClr1(grainColors.at(regionID1));
                        const Eigen::Matrix<int,1,3>& regionClr2(grainColors.at(regionID2));
                        
                        gbPoints->InsertNextPoint ( vM(0,0)+faceOffset*n1(0), vM(1,0)+faceOffset*n1(1), vM(2,0)+faceOffset*n1(2) );
                        gbPoints->InsertNextPoint ( vM(0,1)+faceOffset*n1(0), vM(1,1)+faceOffset*n1(1), vM(2,1)+faceOffset*n1(2) );
                        gbPoints->InsertNextPoint ( vM(0,2)+faceOffset*n1(0), vM(1,2)+faceOffset*n1(1), vM(2,2)+faceOffset*n1(2) );
                        
                        vtkSmartPointer<vtkTriangle> triangle1 = vtkSmartPointer<vtkTriangle>::New();
                        triangle1->GetPointIds()->SetId (0,triPtID+0);
                        triangle1->GetPointIds()->SetId (1,triPtID+1);
                        triangle1->GetPointIds()->SetId (2,triPtID+2);
                        
                        
                        gbTriangles->InsertNextCell ( triangle1 );
                        
                         //                       gbColors->InsertNextTuple3(0,100,100); // use this to assig color to each vertex
//                        gbColors->InsertNextTuple3((regionClr1(0)+regionClr2(0))/2,(regionClr1(1)+regionClr2(1))/2,(regionClr1(2)+regionClr2(2))/2); // use this to assig color to each vertex
                        gbColors->InsertNextTuple3(regionClr1(0),regionClr1(1),regionClr1(2)); // use this to assig color to each vertex
                        
                        
                        triPtID+=3;
                        
                        gbPoints->InsertNextPoint ( vM(0,0)+faceOffset*n2(0), vM(1,0)+faceOffset*n2(1), vM(2,0)+faceOffset*n2(2) );
                        gbPoints->InsertNextPoint ( vM(0,1)+faceOffset*n2(0), vM(1,1)+faceOffset*n2(1), vM(2,1)+faceOffset*n2(2) );
                        gbPoints->InsertNextPoint ( vM(0,2)+faceOffset*n2(0), vM(1,2)+faceOffset*n2(1), vM(2,2)+faceOffset*n2(2) );
                        
                        vtkSmartPointer<vtkTriangle> triangle2 = vtkSmartPointer<vtkTriangle>::New();
                        triangle2->GetPointIds()->SetId (0,triPtID+0);
                        triangle2->GetPointIds()->SetId (1,triPtID+1);
                        triangle2->GetPointIds()->SetId (2,triPtID+2);
                        
                        
                        gbTriangles->InsertNextCell ( triangle2 );
                        
                        //                       gbColors->InsertNextTuple3(0,100,100); // use this to assig color to each vertex
                        gbColors->InsertNextTuple3(regionClr2(0),regionClr2(1),regionClr2(2)); // use this to assig color to each vertex
                        
                        
                        triPtID+=3;
                    }
                }
                
                // to show triangle edges:
                // http://stackoverflow.com/questions/7548966/how-to-display-only-triangle-boundaries-on-textured-surface-in-vtk
                // https://stackoverflow.com/questions/7548966/how-to-display-only-triangle-boundaries-on-textured-surface-in-vtk
                // or loop over SimplexObserver<3,1>
                //                for(const auto& edge : SimplexObserver<3,1>::simplices())
                //                {
                //                    if(edge.second->isRegionBoundarySimplex())
                //                    {
                //
                //                    }
                //                }
            }
            
            gbTrianglePolyData->SetPoints ( gbPoints );
            gbTrianglePolyData->SetPolys ( gbTriangles );
            gbTrianglePolyData->SetPolys ( gbTriangles );
            gbTrianglePolyData->GetCellData()->SetScalars(gbColors);
            
            //                gbTrianglePolyData->GetCellData()->SetScalars(gbColors); // use this to assig color to each vertex
            
            gbMapper->SetInputData(gbTrianglePolyData);
            gbActor->SetMapper(gbMapper);
            gbActor->GetProperty()->SetOpacity(0.05);
//            gbActor->GetProperty()->SetOpacity(1);
            
            //            gbActor->GetProperty()->SetColor(0.0,0.5,0.5);
            
//            clipper->SetInputConnection(gbTriangles->GetOutputPort());
            clipper->SetInputData(gbTrianglePolyData);

            clipper->SetClipFunction(clipPlane);
            clipper->SetValue(0);
            clipper->Update();

            clippedPolyData = clipper->GetOutput();
            clipMapper->SetInputData(clippedPolyData);
            clipActor->SetMapper(clipMapper);
            
            
            renderer->AddActor(gbActor);
//            renderer->AddActor(clipActor);
            
            
            // Update
            update(0/*,0,0.0,Eigen::Matrix<float,3,1>::Zero()*/);
            
            // Render
            renderer->AddActor(actor);
            
        }
        
        /**************************************************************************/
        void update(const long int& frameID
                    /*,const long int& lastFrameID,const float& anglePerStep, Eigen::Matrix<float,3,1> spinAxis*/
        )
        {
            
            if(DispContainerType::isGood(frameID,true))
            {
                dispFileIsGood=true;
                DispContainerType::read(frameID,true);
                
                modifyPts();
                
                
            }
            else
            {
                dispFileIsGood=false;
            }
            
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
            if(dispFileIsGood)
            {
                size_t k=0;
                for (const auto& edge : mesh.observer<1>())
                {
                    if(edge.second->isBoundarySimplex())
                    {
                        DispContainerType::const_iterator iterD1(DispContainerType::find(edge.second->child(0).xID(0)));
                        DispContainerType::const_iterator iterD2(DispContainerType::find(edge.second->child(1).xID(0)));
                        //                        VertexReader<'D',4,float>::const_iterator iterD3(DispContainerType::find(triangleID.second[2]));
                        
                        assert(iterD1!=DispContainerType::end() && "MESH NODE NOT FOUND IN D FILE");
                        assert(iterD2!=DispContainerType::end() && "MESH NODE NOT FOUND IN D FILE");
                        //                        assert(iterD3!=DispContainerType::end() && "MESH NODE NOT FOUND IN D FILE");
                        
                        
                        const float x1=edge.second->child(0).P0(0)+dispCorr*iterD1->second[0];
                        const float y1=edge.second->child(0).P0(1)+dispCorr*iterD1->second[1];
                        const float z1=edge.second->child(0).P0(2)+dispCorr*iterD1->second[2];
                        pts->SetPoint(k,x1,y1,z1);
                        
                        const float x2=edge.second->child(1).P0(0)+dispCorr*iterD2->second[0];
                        const float y2=edge.second->child(1).P0(1)+dispCorr*iterD2->second[1];
                        const float z2=edge.second->child(1).P0(2)+dispCorr*iterD2->second[2];
                        pts->SetPoint(k+1,x2,y2,z2);
                        
                        k=k+2;
                    }
                }
                
                pts->Modified();
            }
            
//            gbActor->GetProperty()->SetColor(0.0,0.5,0.5);
            
        }
        
    };
    
    double SimplicialMeshActor::dispCorr=1.0;
    bool SimplicialMeshActor::showGrainColors=true;
    bool SimplicialMeshActor::showRegionBoundaries=true;
    
} // namespace model
#endif







