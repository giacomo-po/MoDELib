/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationSegmentActor_H_
#define model_DislocationSegmentActor_H_

#include <deque>
#include <string>

#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkMath.h>
#include <vtkProperty.h>
#include <vtkTubeFilter.h>
#include <vtkPolyLine.h>
#include <vtkSphereSource.h>
#include <vtkArrowSource.h>
#include <vtkGlyph3D.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkLabeledDataMapper.h>
#include <vtkFloatArray.h>
#include <vtkLine.h>

//#include <EdgeReader.h>
//#include <vertexReader.h>
#include <IDreader.h>
#include <PlanarPolygon.h>
#include <DDconfigIO.h>
#include <MeshPlane.h>
#include <DDconfigVtkBase.h>

// VTK documentation
// http://www.vtk.org/Wiki/VTK/Examples/Cxx/GeometricObjects/PolyLine
// http://www.vtk.org/Wiki/VTK/Examples/Cxx/VisualizationAlgorithms/TubesWithVaryingRadiusAndColors
// https://stackoverflow.com/questions/45297617/vtk-qt-scene-with-50-moving-actors
namespace model
{
    struct DislocationSegmentActor : public DDconfigVtkBase
    /*                            */,public std::map<std::pair<size_t,size_t>,DislocationSegmentIO<3>>
    {
        static constexpr int dim=3;
        enum ColorScheme {colorBurgers=0,colorSessile=1,colorNormal=2,colorEdgeScrew=3,colorComponent=4};
        typedef Eigen::Matrix<double,dim,1>  VectorDim;
        
        vtkSmartPointer<vtkPolyData> polyData;
        vtkSmartPointer<vtkPolyData> polyDataBnd;
        vtkSmartPointer<vtkPolyData> polyData0;
        vtkSmartPointer<vtkTubeFilter> tubeFilter;
        vtkSmartPointer<vtkTubeFilter> tubeFilterBnd;
        vtkSmartPointer<vtkTubeFilter> tubeFilter0;
        vtkSmartPointer<vtkPolyDataMapper> tubeMapper;
        vtkSmartPointer<vtkActor> tubeActor;
        vtkSmartPointer<vtkPolyDataMapper> tubeMapperBnd;
        vtkSmartPointer<vtkActor> tubeActorBnd;
        vtkSmartPointer<vtkPolyDataMapper> tubeMapper0;
        vtkSmartPointer<vtkActor> tubeActor0;

        std::map<std::pair<size_t,size_t>,DislocationSegmentIO<dim>>& segments()
        {
            return *this;
        }

        const std::map<std::pair<size_t,size_t>,DislocationSegmentIO<dim>>& segments() const
        {
            return *this;
        }

        
        /**********************************************************************/
        DislocationSegmentActor(vtkRenderer* const renderer) :
        /* init */ polyData(vtkSmartPointer<vtkPolyData>::New())
        /* init */,polyDataBnd(vtkSmartPointer<vtkPolyData>::New())
        /* init */,polyData0(vtkSmartPointer<vtkPolyData>::New())
        /* init */,tubeFilter(vtkSmartPointer<vtkTubeFilter>::New())
        /* init */,tubeFilterBnd(vtkSmartPointer<vtkTubeFilter>::New())
        /* init */,tubeMapper(vtkSmartPointer<vtkPolyDataMapper>::New())
        /* init */,tubeActor(vtkSmartPointer<vtkActor>::New())
        /* init */,tubeMapperBnd(vtkSmartPointer<vtkPolyDataMapper>::New())
        /* init */,tubeActorBnd(vtkSmartPointer<vtkActor>::New())
        /* init */,tubeFilter0(vtkSmartPointer<vtkTubeFilter>::New())
        /* init */,tubeMapper0(vtkSmartPointer<vtkPolyDataMapper>::New())
        /* init */,tubeActor0(vtkSmartPointer<vtkActor>::New())
        {
            // Segments
            tubeFilter->SetInputData(polyData);
            tubeFilter->SetRadius(tubeRadius); // this must be a function similar to setColor
            if(scaleRadiusByBurgers)
            {
                tubeFilter->SetVaryRadiusToVaryRadiusByAbsoluteScalar();
            }
            tubeFilter->SetNumberOfSides(10);
            tubeFilter->Update();
            tubeMapper->SetInputConnection(tubeFilter->GetOutputPort());
            tubeMapper->ScalarVisibilityOn();
            tubeMapper->SetScalarModeToUseCellFieldData();
            tubeMapper->SelectColorArray("Colors");
            tubeActor->SetMapper(tubeMapper);
            
            // Boundary segments
            tubeFilterBnd->SetInputData(polyDataBnd);
            tubeFilterBnd->SetRadius(tubeRadius); // this must be a function similar to setColor
            tubeFilterBnd->SetNumberOfSides(10);
            tubeFilterBnd->Update();
            tubeMapperBnd->SetInputConnection(tubeFilterBnd->GetOutputPort());
            tubeMapperBnd->ScalarVisibilityOn();
            tubeActorBnd->SetMapper(tubeMapperBnd);
            tubeActorBnd->GetProperty()->SetOpacity(0.3); //(R,G,B)
            
            // Zero-Burgers Segments
            tubeFilter0->SetInputData(polyData0);
            tubeFilter0->SetRadius(tubeRadius); // this must be a function similar to setColor
            tubeFilter0->SetNumberOfSides(10);
            tubeFilter0->Update();
            tubeMapper0->SetInputConnection(tubeFilter0->GetOutputPort());
            tubeMapper0->ScalarVisibilityOn();
            tubeActor0->SetMapper(tubeMapper0);
            tubeActor0->GetProperty()->SetColor(0.5, 0.5, 0.5); //(R,G,B)
            tubeActor0->GetProperty()->SetOpacity(0.3); //(R,G,B)

            renderer->AddActor(tubeActor);
            renderer->AddActor(tubeActorBnd);
            renderer->AddActor(tubeActor0);

            modify();
        }
        
        
        
        /**********************************************************************/
        void updateConfiguration(const DDconfigIO<3>& configIO, vtkPolyData* const nodePolyData)
        {
            std::cout<<"Updating segments..."<<std::flush;
            const auto t2= std::chrono::system_clock::now();
            segments()=configIO.segments(); // compute and store segments from configIO
            
            vtkSmartPointer<vtkFloatArray> radii(vtkSmartPointer<vtkFloatArray>::New());
            vtkSmartPointer<vtkUnsignedCharArray> colors(vtkSmartPointer<vtkUnsignedCharArray>::New());
            vtkSmartPointer<vtkUnsignedCharArray> colorsBnd(vtkSmartPointer<vtkUnsignedCharArray>::New());
            radii->SetName("TubeRadius");
            colors->SetName("Colors");
            colorsBnd->SetNumberOfComponents(3);
            colors->SetNumberOfComponents(3);
            
            vtkSmartPointer<vtkCellArray> cells(vtkSmartPointer<vtkCellArray>::New());
            vtkSmartPointer<vtkCellArray> cellsBnd(vtkSmartPointer<vtkCellArray>::New());
            vtkSmartPointer<vtkCellArray> cells0(vtkSmartPointer<vtkCellArray>::New());
            
            size_t ptID=0;
            for (const auto& segment : segments())
            {
                auto itSource(configIO.nodeMap().find(segment.second.sourceID)); //source
                assert(itSource!=configIO.nodeMap().end() && "SOURCE VERTEX NOT FOUND IN V-FILE");
                auto   itSink(configIO.nodeMap().find(segment.second.sinkID)); //sink
                assert(  itSink!=configIO.nodeMap().end() &&   "SINK VERTEX NOT FOUND IN V-FILE");

                vtkSmartPointer<vtkLine> line(vtkSmartPointer<vtkLine>::New());
                line->GetPointIds()->SetId(0, std::distance(configIO.nodeMap().begin(),itSource)); // the second 0 is the index of the Origin in linesPolyData's points
                line->GetPointIds()->SetId(1, std::distance(configIO.nodeMap().begin(),itSink));
                
                const auto chord(itSink->second->P-itSource->second->P);
                const double burgersNorm(segment.second.b.norm());
                if(burgersNorm>FLT_EPSILON)
                {
                    Eigen::Matrix<int,3,1> colorVector=computeColor(segment.second.b,chord,segment.second.n);
                    unsigned char lineClr[3]={(unsigned char) colorVector(0),(unsigned char) colorVector(1),(unsigned char) colorVector(2)};
                    
                    if(segment.second.meshLocation!=0)
                    {// 0=MeshLocation::insideMesh
                        cellsBnd->InsertNextCell(line);
//                        polyDataBnd->InsertNextCell(VTK_LINE,2,connectivity); //Connects the first and fourth point we inserted into a line
                        colorsBnd->InsertNextTypedTuple(lineClr);
                    }
                    else
                    {
//                        polyData->InsertNextCell(VTK_LINE,2,connectivity); //Connects the first and fourth point we inserted into a line
                        cells->InsertNextCell(line);
                        if(segment.second.meshLocation==2 && blackGrainBoundarySegments)
                        {
                            unsigned char lineClr1[3]={1,1,1};
                            colors->InsertNextTypedTuple(lineClr1);
                        }
                        else
                        {
                            colors->InsertNextTypedTuple(lineClr);
                        }
                        radii->InsertNextValue(burgersNorm*tubeRadius);
                    }
                }
                else
                {
                    cells0->InsertNextCell(line);
//                    polyData0->InsertNextCell(VTK_LINE,2,connectivity); //Connects the first and fourth point we inserted into a line
                }
                
            }
            
            polyData->SetPoints(nodePolyData->GetPoints());
            polyData->SetLines(cells);
            polyData->GetCellData()->AddArray(colors);
            polyData->GetPointData()->AddArray(radii);
            polyData->GetPointData()->SetActiveScalars("TubeRadius");
            polyData->Modified();
            
            polyDataBnd->SetPoints(nodePolyData->GetPoints());
            polyDataBnd->SetLines(cellsBnd);
            polyDataBnd->GetCellData()->SetScalars(colorsBnd);
            polyDataBnd->Modified();
            
            polyData0->SetPoints(nodePolyData->GetPoints());
            polyData0->SetLines(cells0);
            polyData0->Modified();

            model::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t2)).count()<<" sec]"<<defaultColor<<std::endl;
        }
        
        /**********************************************************************/
        void modify()
        {
            
            tubeFilter->SetRadius(tubeRadius); // this must be a function similar to setColor
            tubeFilterBnd->SetRadius(tubeRadius); // this must be a function similar to setColor
            tubeFilter0->SetRadius(tubeRadius); // this must be a function similar to setColor
            
            if(showBoundarySegments)
            {
                tubeActorBnd->VisibilityOn();
                
            }
            else
            {
                tubeActorBnd->VisibilityOff();
            }
            
            if(showZeroBuergers)
            {
                tubeActor0->VisibilityOn();
                
            }
            else
            {
                tubeActor0->VisibilityOff();
            }
        }
        
        /*********************************************************************/
        Eigen::Matrix<int,3,1> computeColor(const VectorDim& burgers, const VectorDim& chord, const VectorDim& planeNormal)
        {
            
            VectorDim clrVector(VectorDim::Zero());
            
            switch (clr)
            {
                    
                    
                case colorSessile:
                {
                    const bool isSessile= planeNormal.squaredNorm()<FLT_EPSILON || fabs(burgers.dot(planeNormal))>FLT_EPSILON;
                    clrVector(0)= isSessile? 1.0 : 0.1;
                    clrVector(1)= isSessile? 0.5 : 0.4;
                    clrVector(2)= isSessile? 0.0 : 0.9;
                    break;
                }
                    
                case colorNormal:
                {
                    clrVector=planeNormal;
                    break;
                }
                    
                case colorEdgeScrew:
                {
                    const float u = acos(std::fabs(chord.normalized().dot(burgers.normalized())))*2.0/M_PI;
                    //                            RGBcolor rgb(RGBmap::getColor(std::fabs(tubeTangents.col(k).normalized().dot(burgers.normalized())),0,1));
                    //                            colorVector << rgb.r, rgb.g, rgb.b;
                    clrVector=(VectorDim()<<1.0,0.647,0.0).finished()*u+VectorDim::UnitZ()*(1-u);
                    break;
                }
                    //                    break;
                    
                default:
                    clrVector = burgers.normalized();
                    //                    flipColor(colorVector);
                    break;
            }
            
            float clrTol=100.0*FLT_EPSILON;
            if(clrVector(0)<-clrTol)
            {// first component not zero but begative, flip color
                clrVector*=-1.0;
            }
            else if(fabs(clrVector(0))<=clrTol)
            {// first component is zero, use second component
                if(clrVector(1)<-clrTol)
                {// second component not zero but begative, flip color
                    clrVector*=-1.0;
                }
                else if(fabs(clrVector(1))<=clrTol)
                {// second component is zero, use third component
                    if(clrVector(2)<-clrTol)
                    {
                        clrVector*=-1.0;
                    }
                }
            }
            
            clrVector = (clrVector + VectorDim::Ones(dim) * clrVector.norm()).eval();
            clrVector.normalize();
            return (clrVector*255).cast<int>();
        }
        
    };
    
} // namespace model
#endif
