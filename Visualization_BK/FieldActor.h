/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_FieldActor_h_
#define model_FieldActor_h_

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
//#include <VertexReader.h>
#include <IDreader.h>
#include <DDauxIO.h>
#include <MeshPlane.h>
#include <TriangularMesh.h>
#include <DislocationSegmentIO.h>
#include <DislocationNodeIO.h>
#include <StressStraight.h>
#include <DislocationStress.h>

namespace model
{
    
    struct FieldActorCallback : public vtkCommand
    {
        vtkPlane *plane;
        
        FieldActorCallback():
        /* init */,plane(0)
        {
        }
        
        static FieldActorCallback *New()
        {
            return new FieldActorCallback;
            
        }
        
        virtual void Execute(vtkObject *caller, unsigned long, void*)
        {
            vtkImplicitPlaneWidget2 *planeWidget = reinterpret_cast<vtkImplicitPlaneWidget2*>(caller);
            vtkImplicitPlaneRepresentation  *rep = reinterpret_cast<vtkImplicitPlaneRepresentation*>(planeWidget->GetRepresentation());
            rep->GetPlane(this->plane);
        }
    };
    
    class FieldActor : public vtkCommand
    {// https://lorensen.github.io/VTKExamples/site/Cxx/Meshes/ColoredElevationMap/
        // https://lorensen.github.io/VTKExamples/site/Cxx/PolyData/ColorCells/
        
        
    public:
        //        vtkClipPolyData* clipper;
        SimplicialMesh<3>* mesh;
        TriangularMesh planeTriangulation;
        
        /**********************************************************************/
        static FieldActor *New()
        {
            return new FieldActor;
            
        }
        
        /**********************************************************************/
        virtual void Execute(vtkObject *caller, unsigned long, void*)
        {
            vtkImplicitPlaneWidget2 *planeWidget = reinterpret_cast<vtkImplicitPlaneWidget2*>(caller);
            vtkImplicitPlaneRepresentation  *rep = reinterpret_cast<vtkImplicitPlaneRepresentation*>(planeWidget->GetRepresentation());
            rep->GetPlane(this->plane);
            
            //            if(clipper)
            //            {
            
            //                double     origin[3];
            Eigen::Vector3d origin;
            plane->GetOrigin(origin.data());
            //double     normal[3];
            Eigen::Vector3d normal;
            plane->GetNormal(normal.data());
            
            
            //                vtkPolyData* polydata = clipper->GetOutput();
            
            //                std::ofstream pointsFile("points.txt");
            
            //                std::cout<<"callback points="<<polydata->GetNumberOfPoints()<<std::endl;
            //                std::cout << "origin " << " : (" << origin[0] << " " << origin[1] << " " << origin[2] << ")" << std::endl;
            //                std::cout << "normal "  << " : (" << normal[0] << " " << normal[1] << " " << normal[2] << ")" << std::endl;
            //
            //                file<<origin.transpose()<<" "<<normal.transpose()<<std::endl;
            //
            //                Plane<3> projPlane(origin,normal);
            
            MeshPlane<3> meshPlane(*mesh,origin,normal);
            std::deque<Eigen::Matrix<double,2,1>> boundaryPts;
            std::deque<Eigen::Matrix<double,2,1>> internalPts;
            const double meshSize=100.0;
            for(const auto& bndLine : meshPlane.meshIntersections)
            {
                boundaryPts.push_back(meshPlane.localPosition(bndLine->P0));
                //                pointsFile<<meshPlane.localPosition(bndLine->P0).transpose()<<" "<<meshPlane.localPosition(bndLine->P1).transpose()<<std::endl;
            }
            
            //                        if(ddStyle)
            //                        {
            //                            if(ddStyle->ddSegments)
            //                            {
            //                            std::cout<<ddStyle->ddSegments->ddSegments.size()<<" segments"<<std::endl;
            //                            }
            //                        }
            
            //std::vector<std::pair<double,int>> Rs{{0.5,8},{1.0,8},{2.0,8},{4.0,8},{8.0,16},{16.0,16},{32.0,32}};
            
            //            for (const auto& segment : *ddSegments)
            //            {
            //                auto itSource(nodeMap->find(segment.second.sourceID)); //source
            //                assert(itSource!=nodeMap->end() && "SOURCE VERTEX NOT FOUND IN V-FILE");
            //                auto   itSink(nodeMap->find(segment.second.sinkID)); //sink
            //                assert(  itSink!=nodeMa->end() &&   "SINK VERTEX NOT FOUND IN V-FILE");
            //
            //                PlaneSegmentIntersection<3> psi(origin,normal,itSource->second->P,itSink->second->P);
            //                if(psi.type==PlaneSegmentIntersection<3>::INCIDENT)
            //                {
            //                    internalPts.push_back(meshPlane.localPosition(psi.x0));
            //
            ////                    for(int n=0;n<10;++n)
            ////                    {
            //////                        const double& r(Rs[n].first);
            //////                        const int& np(Rs[n].second);
            ////                        const double r(n);
            ////                        const int np(std::max(20.0,std::round(2.0*M_PI*r/20.0)));
            ////
            ////                        for(int k=0;k<np;++k)
            ////                        {
            ////                            double theta=k*2.0*M_PI/np;
            ////
            ////                            internalPts.push_back(meshPlane.localPosition(psi.x0)+(Eigen::Vector2d()<<r*cos(theta),r*sin(theta)).finished());
            ////                        }
            ////                    }
            //                }
            //            }
            
            
            planeTriangulation.reMesh(boundaryPts,internalPts,meshSize);
            
            
            //            std::cout<<"DONE MESHING"<<std::endl;
            
            
            vtkSmartPointer<vtkPoints> meshPts(vtkSmartPointer<vtkPoints>::New());
            
            //            std::ofstream pointsFile("points.txt");
            for(const auto& point2d : planeTriangulation.vertices())
            {
                //                pointsFile<<meshPlane.globalPosition(point2d).transpose()<<"\n";
                const auto point3d(meshPlane.globalPosition(point2d));
                meshPts->InsertNextPoint(point3d(0),point3d(1),point3d(2));
            }
            
            vtkSmartPointer<vtkCellArray> meshTriangles(vtkSmartPointer<vtkCellArray>::New());
            vtkSmartPointer<vtkUnsignedCharArray> meshColors(vtkSmartPointer<vtkUnsignedCharArray>::New());
            meshColors->SetNumberOfComponents(3);
            for(const auto& tri : planeTriangulation.triangles())
            {
                vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
                triangle->GetPointIds()->SetId (0,tri(0));
                triangle->GetPointIds()->SetId (1,tri(1));
                triangle->GetPointIds()->SetId (2,tri(2));
                meshTriangles->InsertNextCell ( triangle );
                const auto triColor(Eigen::Matrix<int,1,3>::Random()*255);
                meshColors->InsertNextTuple3(triColor(0),triColor(1),triColor(2)); // use this to assig color to each vertex
            }
            
            vtkSmartPointer<vtkUnsignedCharArray> fieldColors(vtkSmartPointer<vtkUnsignedCharArray>::New());
            fieldColors->SetNumberOfComponents(3);
            for(const auto& vtx : planeTriangulation.vertices())
            {
                Eigen::Matrix<double,3,3> stress(Eigen::Matrix<double,3,3>::Zero());
                for (const auto& segment : *ddSegments)
                {
                    auto itSource(nodeMap->find(segment.second.sourceID)); //source
                    assert(itSource!=nodeMap->end() && "SOURCE VERTEX NOT FOUND IN V-FILE");
                    auto   itSink(nodeMap->find(segment.second.sinkID)); //sink
                    assert(  itSink!=nodeMa->end() &&   "SINK VERTEX NOT FOUND IN V-FILE");
                    StressStraight<3> ss(itSource->second->P,itSink->second->P,segment.second.b);
                    stress+=ss.stress(meshPlane.globalPosition(vtx));
                }
                double dclr[3];
                lut->GetColor(stress(0.0)*1000, dclr);
                unsigned char cclr[3];
                for(unsigned int j = 0; j < 3; j++)
                {
                    cclr[j] = static_cast<unsigned char>(255.0 * dclr[j]);
                }
                fieldColors->InsertNextTypedTuple(cclr);
                
            }
            
            
            
            meshPolydata->SetPoints ( meshPts );
            meshPolydata->SetPolys ( meshTriangles );
            meshPolydata->GetPointData()->SetScalars(fieldColors);
            //meshPolydata->GetCellData()->SetScalars(meshColors);
            
            //            gbMapper->SetInputData(gbTrianglePolyData);
            //            gbActor->SetMapper(gbMapper);
            
            meshPolydata->Modified();
        }
        
        /**********************************************************************/
        FieldActor():
        /* init */ lut(vtkSmartPointer<vtkLookupTable>::New())
        /* init */,plane(0)
//        /* init */,actor(0)
        /* init */,meshPolydata(vtkSmartPointer<vtkPolyData>::New())
        /* init */,meshMapper(vtkSmartPointer<vtkPolyDataMapper>::New())
        /* init */,meshActor(vtkSmartPointer<vtkActor>::New())
        {
            
            //            vtkSmartPointer<vtkLookupTable> lut = vtkSmartPointer<vtkLookupTable>::New();
            //            int tableSize = std::max(resolution*resolution + 1, 10);
            lut->SetTableRange(-5.0, 5.0);
            if(false)
            {
                lut->SetNumberOfTableValues(3);
                vtkSmartPointer<vtkNamedColors> colors = vtkSmartPointer<vtkNamedColors>::New();
                //                lut->SetTableValue(0, colors->GetColor4d("Black").GetData());
                //                lut->SetTableValue(1, colors->GetColor4d("Banana").GetData());
                lut->SetTableValue(0, colors->GetColor4d("Tomato").GetData());
                lut->SetTableValue(1, colors->GetColor4d("White").GetData());
                lut->SetTableValue(2, colors->GetColor4d("Lavender").GetData());
                
                //                lut->SetTableValue(3, colors->GetColor4d("Wheat").GetData());
                //                lut->SetTableValue(4, colors->GetColor4d("Lavender").GetData());
                //                lut->SetTableValue(5, colors->GetColor4d("Flesh").GetData());
                //                lut->SetTableValue(6, colors->GetColor4d("Raspberry").GetData());
                //                lut->SetTableValue(7, colors->GetColor4d("Salmon").GetData());
                //                lut->SetTableValue(8, colors->GetColor4d("Mint").GetData());
                //                lut->SetTableValue(9, colors->GetColor4d("Peacock").GetData());
                
            }
            lut->Build();
            
            // Fill in a few known colors, the rest will be generated if needed
            
            meshPolydata->Allocate();
            meshMapper->SetInputData(meshPolydata);
            //            meshMapper->SetScalarRange(-1, 1);
            //           meshMapper->SetLookupTable(lut);
            
            meshActor->SetMapper ( meshMapper );
            //            meshActor->GetProperty()->SetLineWidth(0.5);
            //            meshActor->GetProperty()->SetColor(0.5,0.5,0.5); // Give some color to the mesh. (1,1,1) is white
            meshActor->GetProperty()->SetOpacity(0.8); //Make the mesh have some transparency.
        }
        
        vtkSmartPointer<vtkLookupTable> lut;
        vtkPlane *plane;
//        vtkActor *actor;
        std::map<std::pair<size_t,size_t>,DislocationSegmentIO<3>>* ddSegments;
        std::map<size_t,const DislocationNodeIO<3>* const>* nodeMap;
        //        DDinteractionStyle* ddStyle;
        //        vtkSmartPointer<vtkPoints> meshPts;
        vtkSmartPointer<vtkPolyData> meshPolydata;
        vtkSmartPointer<vtkPolyDataMapper> meshMapper;
        vtkSmartPointer<vtkActor> meshActor;
        
        
    };
}
#endif

