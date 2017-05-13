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

#include <model/Mesh/SimplicialMesh.h>


namespace model
{

    struct SimplicialMeshActor //: public vtkSmartPointer<vtkActor>
    {
        
        vtkSmartPointer<vtkPoints> pts;
        vtkSmartPointer<vtkPolyData> polydata;
        vtkSmartPointer<vtkPolyDataMapper> mapper;
        
        SimplicialMesh<3> mesh;
        
        /**************************************************************************/
//        SimplicialMeshActor(const SimplicialMesh<3>& mesh) :
        SimplicialMeshActor(const int& meshID) :
        //    vtkSmartPointer<vtkActor>(vtkSmartPointer<vtkActor>::New()),
        /* init */ pts(vtkSmartPointer<vtkPoints>::New()),
        /* init */ polydata(vtkSmartPointer<vtkPolyData>::New()),
        /* init */ mapper(vtkSmartPointer<vtkPolyDataMapper>::New())
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
            
            //        this->SetMapper ( mapper );
        }
        
        /**************************************************************************/
        vtkSmartPointer<vtkActor> actor() const
        {
            vtkSmartPointer<vtkActor> a =  vtkSmartPointer<vtkActor>::New();
            a->SetMapper ( mapper );
            a->GetProperty()->SetColor(0.5,0.5,0.5); // Give some color to the mesh
            a->GetProperty()->SetOpacity(0.5); //Make the tube have some transparency.
            return a;
        }
    };
	
	
} // namespace model
#endif







