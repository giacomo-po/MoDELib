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
#include <vtkMath.h>
#include <vtkProperty.h>
#include <vtkTubeFilter.h>

//#include <model/DislocationDynamics/Visualization/vtk/vtkActorWrapper.h>

//#include <model/Mesh/SimplicialMesh.h>


namespace model
{
    struct DislocationSegmentActor
    {
        static constexpr int dim=3;
        typedef Eigen::Matrix<float,dim,1>  VectorDim;

        enum ColorScheme {colorBurgers=0,colorSessile=1,colorNormal=2,colorEdgeScrew=3,colorComponent=4};

        
        VectorDim planeNormal;
        VectorDim burgers;
        VectorDim chord;
        VectorDim colorVector;
        
        //    http://www.vtk.org/Wiki/VTK/Examples/Cxx/VisualizationAlgorithms/TubesWithVaryingRadiusAndColors
        
        vtkSmartPointer<vtkPoints> points;
        vtkSmartPointer<vtkCellArray> lines;
        vtkSmartPointer<vtkPolyData> polyData;
        vtkSmartPointer<vtkPolyDataMapper> lineMapper;
        vtkSmartPointer<vtkTubeFilter> tubeFilter;
        vtkSmartPointer<vtkPolyDataMapper> tubeMapper;
        
        vtkSmartPointer<vtkActor> line;
        vtkSmartPointer<vtkActor> tube;
//        vtkSmartPointer<vtkActorWrapper<double>> tube;
        
        static float tubeRadius;
        static ColorScheme clr;

        /*********************************************************************/
        void computeColor()
        {
            
            switch (clr)
            {
                    //                case colorSessile:
                    //                    colorVector(0)= isSessile? 1.0 : 0.1;
                    //                    colorVector(1)= isSessile? 0.5 : 0.4;
                    //                    colorVector(2)= isSessile? 0.0 : 0.9;
                    //                    break;
                    
                case colorNormal:
                    colorVector = planeNormal;
//                    flipColor(colorVector);
                    break;
                    
                    //                case colorComponent:
                    //                {
                    //                    RGBcolor rgb(RGBmap::getColor(ids,sIDmin,sIDmax));
                    //                    colorVector << rgb.r, rgb.g, rgb.b;
                    //                }
                    //                    break;
                    
                    //                case colorEdgeScrew:
                    //                {
                    //                    const float u = std::fabs(tubeTangents.col(k).normalized().dot(burgers.normalized()));
                    //                    //                            RGBcolor rgb(RGBmap::getColor(std::fabs(tubeTangents.col(k).normalized().dot(burgers.normalized())),0,1));
                    //                    //                            colorVector << rgb.r, rgb.g, rgb.b;
                    //                    colorVector=screwColor*u+edgeColor*(1-u);
                    //                }
                    //                    break;
                    
                default:
                    colorVector = burgers.normalized();
//                    flipColor(colorVector);
                    break;
            }
            
            if(colorVector(0)<0.0)
            {
                colorVector*=-1;
            }
            else if(colorVector(0)==0.0)
            {
                if(colorVector(1)<0.0)
                {
                    colorVector*=-1;
                }
                else if(colorVector(1)==0.0)
                {
                    if(colorVector(2)<0.0)
                    {
                        colorVector*=-1;
                    }
                }
            }
            
            //			VectorDim colorVector = burgers + VectorDim::Ones(dim) * burgers.norm();
            colorVector = (colorVector + VectorDim::Ones(dim) * colorVector.norm()).eval();
            
            //		colorVector << 0.0f,0.6f,0.4f;
            colorVector.normalize();
            

        }
        
        /************************************************************************/
        DislocationSegmentActor(const Eigen::Matrix<float,dim,6>& P0T0P1T1BN) :
        /* init */ points(vtkSmartPointer<vtkPoints>::New()),
        /* init */ lines(vtkSmartPointer<vtkCellArray>::New()),
        /* init */ polyData(vtkSmartPointer<vtkPolyData>::New()),
        /* init */ lineMapper(vtkSmartPointer<vtkPolyDataMapper>::New()),
        /* init */ tubeFilter(vtkSmartPointer<vtkTubeFilter>::New()),
        /* init */ tubeMapper(vtkSmartPointer<vtkPolyDataMapper>::New()),
        /* init */ line(vtkSmartPointer<vtkActor>::New()),
        /* init */ tube(vtkSmartPointer<vtkActor>::New())
        {
            
            float alpha=0.5;
            
            chord = P0T0P1T1BN.col(2)-P0T0P1T1BN.col(0);
            burgers=P0T0P1T1BN.col(4);
            planeNormal=P0T0P1T1BN.col(5);
            float g = std::pow(chord.norm(),alpha);
            
            unsigned int Np = 10;      // No. of vertices

            for (int k=0;k<Np;++k) // this may have to go to Np+1
            {
                const float u1=k*1.0/(Np-1);
                const float u2=u1*u1;
                const float u3=u2*u1;
                
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
            line->GetProperty()->SetColor(0.0,1.0,0.0); // Give some color to the tube
            line->SetMapper ( lineMapper );

            if(true)
            {
                tubeFilter->SetInputData(polyData);
                tubeFilter->SetRadius(tubeRadius); // this must be a function similar to setColor
                tubeFilter->SetNumberOfSides(10);
                tubeFilter->Update();
                tubeMapper->SetInputConnection(tubeFilter->GetOutputPort());
                tubeMapper->ScalarVisibilityOn();
                tube->SetMapper(tubeMapper);
                tube->GetProperty()->SetOpacity(1.0); //Make the tube have some transparency.
                computeColor();
                tube->GetProperty()->SetColor(colorVector(0),colorVector(1),colorVector(2)); // Give some color to the tube

            }
            
        }
        
        /**************************************************************************/
        vtkSmartPointer<vtkActor> lineActor() const
        {
            return line;
        }
        
        /**************************************************************************/
        vtkSmartPointer<vtkActor> tubeActor() const
        {
             return tube;
        }
        
        void modify()
        {
            line->GetProperty()->SetColor(0.0,1.0,0.0); // Give some color to the tube
//            line->Modified();
            
            tubeFilter->SetRadius(tubeRadius); // this must be a function similar to setColor
            computeColor();
            tube->GetProperty()->SetColor(colorVector(0),colorVector(1),colorVector(2)); // Give some color to the tube
//            tube->Modified();

        }
        
    };
    
    float DislocationSegmentActor::tubeRadius=5.0;
    DislocationSegmentActor::ColorScheme DislocationSegmentActor::clr=DislocationSegmentActor::colorBurgers;
} // namespace model
#endif







