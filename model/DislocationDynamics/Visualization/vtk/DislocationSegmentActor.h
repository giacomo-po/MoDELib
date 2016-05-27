/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationSegmentActor_H_
#define model_DislocationSegmentActor_H_

#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>

//#include <model/Mesh/SimplicialMesh.h>


namespace model
{
    struct DislocationSegmentActor
    {
        static constexpr int dim=3;
        typedef Eigen::Matrix<float,dim,1>  VectorDim;

        VectorDim planeNormal;
        VectorDim burgers;
        VectorDim chord;
        
        //    http://www.vtk.org/Wiki/VTK/Examples/Cxx/VisualizationAlgorithms/TubesWithVaryingRadiusAndColors
        
        vtkSmartPointer<vtkPoints> points;
        vtkSmartPointer<vtkCellArray> lines;
        vtkSmartPointer<vtkPolyData> polyData;
        vtkSmartPointer<vtkPolyDataMapper> lineMapper;
        vtkSmartPointer<vtkTubeFilter> tubeFilter;
        vtkSmartPointer<vtkPolyDataMapper> tubeMapper;
        
        
        DislocationSegmentActor(const Eigen::Matrix<float,dim,6>& P0T0P1T1BN) :
        /* init */ points(vtkSmartPointer<vtkPoints>::New()),
        /* init */ lines(vtkSmartPointer<vtkCellArray>::New()),
        /* init */ polyData(vtkSmartPointer<vtkPolyData>::New()),
        /* init */ lineMapper(vtkSmartPointer<vtkPolyDataMapper>::New()),
        /* init */ tubeFilter(vtkSmartPointer<vtkTubeFilter>::New()),
        /* init */ tubeMapper(vtkSmartPointer<vtkPolyDataMapper>::New())
        {
            
            float alpha=0.5;
            
            chord = P0T0P1T1BN.col(2)-P0T0P1T1BN.col(0);
            float g = std::pow(chord.norm(),alpha);
            
            unsigned int Np = 10;      // No. of vertices

            for (int k=0;k<Np;++k)
            {
                float u1=k*1.0/(Np-1);
                float u2=u1*u1;
                float u3=u2*u1;
                
                // Compute positions along axis
                VectorDim P =   ( 2.0f*u3-3.0f*u2+1.0f) * P0T0P1T1BN.col(0)
                /*************/ + g*(      u3-2.0f*u2+u1)   * P0T0P1T1BN.col(1)
                /*************/ +   (-2.0f*u3+3.0f*u2)      * P0T0P1T1BN.col(2)
                /*************/ + g*(      u3-u2)           * P0T0P1T1BN.col(3);
                
//                // Compute tangents along axis
//                tubeTangents.col(k)= (     ( 6.0f*u2-6.0f*u1)      * P0T0P1T1BN.col(0)
//                                      /*                  */ + g*( 3.0f*u2-4.0f*u1+1.0f) * P0T0P1T1BN.col(1)
//                                      /*                  */ +   (-6.0f*u2+6.0f*u1)      * P0T0P1T1BN.col(2)
//                                      /*                  */ + g*( 3.0f*u2-2.0f*u1)      * P0T0P1T1BN.col(3) ).normalized();
//                
//                // Compute unit vectors in radial direction orthogonal to axis
//                tubeCircles.push_back(getCircle(k));
                points->InsertPoint(k, P(0), P(1), P(2));

            }
            
            lines->InsertNextCell(Np);
            for (int i = 0; i < Np; i++)
            {
                lines->InsertCellPoint(i);
            }

            
            polyData->SetPoints(points);
            polyData->SetLines(lines);
            
            //        lineMapper->SetInputConnection(lines->GetOutputPort());
            lineMapper->SetInputData(polyData);
            
            if(false)
            {
            tubeFilter->SetInputData(polyData);
            tubeFilter->SetRadius(10); //default is .5
            tubeFilter->SetNumberOfSides(50);
            tubeFilter->Update();
            tubeMapper->SetInputConnection(tubeFilter->GetOutputPort());
            tubeMapper->ScalarVisibilityOn();
            }
            
            //        tubeMapper->SetScalarModeToUsePointFieldData();
            //        tubeMapper->SelectColorArray("Colors");
            
        }
        
        /**************************************************************************/
        vtkSmartPointer<vtkActor> lineActor() const
        {
            vtkSmartPointer<vtkActor> a =  vtkSmartPointer<vtkActor>::New();
            a->GetProperty()->SetColor(0.0,1.0,0.0); // Give some color to the tube
            a->SetMapper ( lineMapper );
            return a;
        }
        
        /**************************************************************************/
        vtkSmartPointer<vtkActor> tubeActor() const
        {
            vtkSmartPointer<vtkActor> a = vtkSmartPointer<vtkActor>::New();
            a->SetMapper(tubeMapper);
            a->GetProperty()->SetColor(1.0,0.0,0.0); // Give some color to the tube
            a->GetProperty()->SetOpacity(0.5); //Make the tube have some transparency.
            return a;
        }
        
    };
    
	
} // namespace model
#endif







