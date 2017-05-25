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

#include <model/Mesh/SimplicialMesh.h>
#include <model/Network/Readers/VertexReader.h>


namespace model
{

    struct SimplicialMeshActor : public VertexReader<'D',4,float>//: public vtkSmartPointer<vtkActor>
    {
        
        typedef VertexReader<'D',4,float> DispContainerType;

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
        
        
        SimplicialMesh<3> mesh;
        bool dispFileIsGood;
        
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
        /* init */ dispFileIsGood(false)
        {
            
        }
        
        /**************************************************************************/
        void init(const int& meshID,vtkRenderer* renderer)
        {
            mesh.readMesh(meshID);
            
            polydata->Allocate();
            
            size_t connectivityID=0;
            for (const auto& edge : SimplexObserver<3,1>::simplices())
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
            actor->GetProperty()->SetLineWidth(0.1);
            actor->GetProperty()->SetColor(0.5,0.5,0.5); // Give some color to the mesh. (1,1,1) is white
            actor->GetProperty()->SetOpacity(0.15); //Make the mesh have some transparency.

            
            // grain-boundaries
            bool showRegionBoundaries=true;
            if(showRegionBoundaries)
            {
                unsigned char clr[3] = {0, 100, 100};
                gbColors->SetNumberOfComponents(3);
                gbColors->SetName("Colors");
//                gbColors->InsertNextTypedTuple(red);
                size_t triPtID=0;
                for(const auto& meshTriangle : SimplexObserver<3,2>::simplices())
                {
                    if(meshTriangle.second->isRegionBoundarySimplex())
                    {
                        const Eigen::Matrix<double,3,3> vM(meshTriangle.second->vertexPositionMatrix());
                        
                        gbPoints->InsertNextPoint ( vM(0,0), vM(1,0), vM(2,0) );
                        gbPoints->InsertNextPoint ( vM(0,1), vM(1,1), vM(2,1) );
                        gbPoints->InsertNextPoint ( vM(0,2), vM(1,2), vM(2,2) );

                        vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
                        triangle->GetPointIds()->SetId (0,triPtID+0);
                        triangle->GetPointIds()->SetId (1,triPtID+1);
                        triangle->GetPointIds()->SetId (2,triPtID+2);

//                        gbColors->InsertNextTypedTuple(clr); // use this to assig color to each vertex
//                        gbColors->InsertNextTypedTuple(clr); // use this to assig color to each vertex
//                        gbColors->InsertNextTypedTuple(clr); // use this to assig color to each vertex

                        gbTriangles->InsertNextCell ( triangle );

                        triPtID+=3;
                    }
                }
                gbTrianglePolyData->SetPoints ( gbPoints );
                gbTrianglePolyData->SetPolys ( gbTriangles );
//                gbTrianglePolyData->GetCellData()->SetScalars(gbColors); // use this to assig color to each vertex

                gbMapper->SetInputData(gbTrianglePolyData);
                gbActor->SetMapper(gbMapper);
                gbActor->GetProperty()->SetOpacity(0.3);
                gbActor->GetProperty()->SetColor(0.0,0.4,0.4);


                renderer->AddActor(gbActor);
                
                // to show triangle edges:
                // http://stackoverflow.com/questions/7548966/how-to-display-only-triangle-boundaries-on-textured-surface-in-vtk
                // or loop over SimplexObserver<3,1>
//                for(const auto& edge : SimplexObserver<3,1>::simplices())
//                {
//                    if(edge.second->isRegionBoundarySimplex())
//                    {
//                        
//                    }
//                }
            }


            // Update
            update(meshID);
            
            // Render
            renderer->AddActor(actor);
            
        }
        
        /**************************************************************************/
        void update(const long int& frameID)
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

            
        }
        
        /**************************************************************************/
        void modifyPts()
        {
            if(dispFileIsGood)
            {
                size_t k=0;
                for (const auto& edge : SimplexObserver<3,1>::simplices())
                {
                    if(edge.second->isBoundarySimplex())
                    {
                        DispContainerType::const_iterator iterD1(DispContainerType::find(edge.second->child(0).xID(0)));
                        DispContainerType::const_iterator iterD2(DispContainerType::find(edge.second->child(1).xID(0)));
                        //                        VertexReader<'D',4,float>::const_iterator iterD3(DispContainerType::find(triangleID.second[2]));
                        
                        assert(iterD1!=DispContainerType::end() && "MESH NODE NOT FOUND IN D FILE");
                        assert(iterD2!=DispContainerType::end() && "MESH NODE NOT FOUND IN D FILE");
                        //                        assert(iterD3!=DispContainerType::end() && "MESH NODE NOT FOUND IN D FILE");
                        
                        
                        const float x1=edge.second->child(0).P0(0)+dispCorr*iterD1->second(0);
                        const float y1=edge.second->child(0).P0(1)+dispCorr*iterD1->second(1);
                        const float z1=edge.second->child(0).P0(2)+dispCorr*iterD1->second(2);
                        pts->SetPoint(k,x1,y1,z1);
                        
                        const float x2=edge.second->child(1).P0(0)+dispCorr*iterD2->second(0);
                        const float y2=edge.second->child(1).P0(1)+dispCorr*iterD2->second(1);
                        const float z2=edge.second->child(1).P0(2)+dispCorr*iterD2->second(2);
                        pts->SetPoint(k+1,x2,y2,z2);
                        
                        k=k+2;
                    }
                }
                
                pts->Modified();
            }

        }
        
    };
	
    double SimplicialMeshActor::dispCorr=1.0;
    
} // namespace model
#endif







