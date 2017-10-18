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

//#include <model/IO/EdgeReader.h>
//#include <model/IO/vertexReader.h>
#include <model/IO/IDreader.h>

// VTK documentation
// http://www.vtk.org/Wiki/VTK/Examples/Cxx/GeometricObjects/PolyLine
// http://www.vtk.org/Wiki/VTK/Examples/Cxx/VisualizationAlgorithms/TubesWithVaryingRadiusAndColors

namespace model
{
    struct DislocationSegmentActor :
    /* inherits from   */ public IDreader<'V',1,10, double>,
    /* inherits from   */ public IDreader<'K',2,13,double>
    {
        
        //    public:
        
        static constexpr int dim=3;
        enum ColorScheme {colorBurgers=0,colorSessile=1,colorNormal=2,colorEdgeScrew=3,colorComponent=4};
        //        typedef VertexReader<'V',10,double> VertexReaderType;
        //        typedef EdgeReader<'K',15,double>   EdgeReaderType;
        typedef IDreader<'V',1,10, double> VertexReaderType;
        typedef IDreader<'K',2,13,double> EdgeReaderType;
        
        typedef Eigen::Matrix<float,dim,1>  VectorDim;
        
        static float alpha;
        static float tubeRadius;
        static ColorScheme clr;
        static size_t Np;      // No. of vertices per line
        static bool plotBoundarySegments;
        static bool showVelocities;
        static bool showNodeIDs;
        static float velocityFactor;
        static bool showZeroBuergers;
        
        vtkRenderer* const renderer;
        
        vtkSmartPointer<vtkActor> lineActor;
        vtkSmartPointer<vtkActor> tubeActor;
        
        vtkSmartPointer<vtkActor> lineActor0;
        vtkSmartPointer<vtkActor> tubeActor0;
        
        
        //    private:
        
        VectorDim planeNormal;
        VectorDim burgers;
        VectorDim chord;
        VectorDim colorVector;
        
        
        // segments objects
        vtkSmartPointer<vtkPoints> points;
        //        std::deque<vtkSmartPointer<vtkPolyLine>> lines;
        vtkSmartPointer<vtkCellArray> cells;
        vtkSmartPointer<vtkCellArray> cells0;
        vtkSmartPointer<vtkPolyData> polyData;
        vtkSmartPointer<vtkPolyData> polyData0;
        vtkSmartPointer<vtkUnsignedCharArray> colors;
        vtkSmartPointer<vtkPolyDataMapper> lineMapper;
        vtkSmartPointer<vtkPolyDataMapper> lineMapper0;
        vtkSmartPointer<vtkTubeFilter> tubeFilter;
        vtkSmartPointer<vtkTubeFilter> tubeFilter0;
        vtkSmartPointer<vtkPolyDataMapper> tubeMapper;
        vtkSmartPointer<vtkPolyDataMapper> tubeMapper0;
        
        // node objects
        vtkSmartPointer<vtkPoints> nodePoints;
        //        vtkSmartPointer<vtkStringArray> nodeLabels;
        vtkSmartPointer<vtkSphereSource> sphereSource;
        vtkSmartPointer<vtkPolyData> nodeData;
        vtkSmartPointer<vtkGlyph3D> nodeGlyphs;
        vtkSmartPointer<vtkPolyDataMapper> nodeMapper;
        vtkSmartPointer<vtkActor> nodeActor;
        
        // velocity objects
        //        vtkSmartPointer<vtkPoints> velocityPoints;
        vtkSmartPointer<vtkDoubleArray> velocityVectors;
        vtkSmartPointer<vtkUnsignedCharArray> velocityColors;
        vtkSmartPointer<vtkPolyData> velocityPolyData;
        vtkSmartPointer<vtkArrowSource> velocityArrowSource;
        vtkSmartPointer<vtkGlyph3D> velocityGlyphs;
        vtkSmartPointer<vtkPolyDataMapper> velocityMapper;
        vtkSmartPointer<vtkActor> velocityActor;
        
        vtkSmartPointer<vtkPolyData> labelPolyData;
        vtkSmartPointer<vtkDoubleArray> labelScalars;
        vtkSmartPointer<vtkLabeledDataMapper> labelMapper;
        vtkSmartPointer<vtkActor2D> labelActor;
        
        
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
        
        
        
        //    public:
        
        /**********************************************************************/
        VertexReaderType& vertexReader()
        {
            return *this;
        }
        
        /**********************************************************************/
        EdgeReaderType& edgeReader()
        {
            return *this;
        }
        
        /**********************************************************************/
        void readNodes(const size_t& frameID)
        {
            if (vertexReader().isGood(frameID,false)) // bin format
            {
                vertexReader().read(frameID,false);
            }
            else // txt format
            {
                vertexReader().read(frameID,true);
            }
            
            //            nodeLabels->SetNumberOfValues(vertexReader().size());
            //            size_t labelID=0;
            for(const auto& node : vertexReader())
            {
                Eigen::Map<const Eigen::Matrix<double,1,10>> row(node.second.data());
                
                nodePoints->InsertNextPoint(row.template segment<dim>(0).data());
                velocityVectors->InsertNextTuple(row.template segment<dim>(dim).data()); // arrow vactor
                unsigned char clr[3]={255,0,255};
                velocityColors->InsertNextTypedTuple(clr);
                
                
                labelScalars->InsertNextTuple1(node.first);
                
                //                nodeLabels->SetValue(labelID, std::to_string(node.first));
                //                labelID++;
                
            }
            
            nodeData->SetPoints(nodePoints);
            //            nodeData->GetOutput()->GetPointData()->AddArray(labels);
            
            //nodeData->GetPointData()->SetVectors(vectors);
            nodeData->Modified();
            
            velocityPolyData->SetPoints(nodePoints);
            velocityPolyData->GetPointData()->SetVectors(velocityVectors);
            velocityPolyData->Modified();
            velocityPolyData->GetCellData()->SetScalars(velocityColors);
            velocityPolyData->Modified();
            
            //            labels->SetValue(1, "Priority 7");
            //            labels->SetValue(2, "Priority 6");
            //            labels->SetValue(3, "Priority 4");
            //            labels->SetValue(4, "Priority 4");
            //            labels->SetValue(5, "Priority 4");
            //            nodeData->GetOutput()->GetPointData()->AddArray(labels);
            
            labelPolyData->SetPoints(nodePoints);
            labelPolyData->GetPointData()->SetScalars(labelScalars);
        }
        
        
        /**********************************************************************/
        void readSegments(const size_t& frameID)
        {
            if (edgeReader().isGood(frameID,false)) // bin format
            {
                edgeReader().read(frameID,false);
            }
            else // txt format
            {
                edgeReader().read(frameID,true);
            }
            
            size_t ptID=0;
            for (const auto& edge : edgeReader())
            {
                
                VertexReaderType::const_iterator itSource(vertexReader().find(edge.first[0])); //source
                assert(itSource!=vertexReader().end() && "SOURCE VERTEX NOT FOUND IN V-FILE");
                VertexReaderType::const_iterator   itSink(vertexReader().find(edge.first[1])); //sink
                assert(  itSink!=vertexReader().end() &&   "SINK VERTEX NOT FOUND IN V-FILE");
                
                Eigen::Map<const Eigen::Matrix<double,1,6>> sourceRow(itSource->second.data());
                Eigen::Map<const Eigen::Matrix<double,1,6>>   sinkRow(  itSink->second.data());
                Eigen::Map<const Eigen::Matrix<double,1,13>>   edgeRow(edge.second.data());
                
                const int   snID(edgeRow(2*dim+2));
                const bool sourceOnBoundary(sourceRow(2*dim+1));
                const bool   sinkOnBoundary(  sinkRow(2*dim+1));
                
                if(!(sourceOnBoundary && sinkOnBoundary) || plotBoundarySegments)
                {
                    
                    
                    Eigen::Matrix<float,dim,6> P0T0P1T1BN;
                    
                    P0T0P1T1BN.col(0) = sourceRow.segment<dim>(0*dim).transpose().template cast<float>();	// source position
                    P0T0P1T1BN.col(2) =   sinkRow.segment<dim>(0*dim).transpose().template cast<float>();	// sink position
                    //                    P0T0P1T1BN.col(1) = sourceTfactor*(itSource->second.segment<dim>(1*dim).transpose().template cast<float>());	// source tangent
                    //                    P0T0P1T1BN.col(3) =  -sinkTfactor*(  itSink->second.segment<dim>(1*dim).transpose().template cast<float>());	// sink tangent
                    P0T0P1T1BN.col(1) = edgeRow.segment<dim>(2*dim).transpose().template cast<float>();	// source tangent
                    P0T0P1T1BN.col(3) = edgeRow.segment<dim>(3*dim).transpose().template cast<float>();	// sink tangent
                    P0T0P1T1BN.col(4) = edgeRow.segment<dim>(0*dim).transpose().template cast<float>();		// Burgers vector
                    P0T0P1T1BN.col(5) = edgeRow.segment<dim>(1*dim).transpose().template cast<float>();		// plane normal
                    chord = P0T0P1T1BN.col(2)-P0T0P1T1BN.col(0);
                    burgers=P0T0P1T1BN.col(4);
                    planeNormal=P0T0P1T1BN.col(5);
                    const float g = std::pow(chord.norm(),alpha);
                    
                    
                    //                    lines.push_back(vtkSmartPointer<vtkPolyLine>::New());
                    //                    auto& line(*lines.rbegin());
                    vtkSmartPointer<vtkPolyLine> line=vtkSmartPointer<vtkPolyLine>::New();
                    line->GetPointIds()->SetNumberOfIds(Np);
                    unsigned char clr[3]={51,153,255};
                    //                    unsigned char clr0[3]={255,255,255};
                    
                    colors->InsertNextTypedTuple(clr);
                    
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
                        
                        points->InsertNextPoint(P.data());
                        line->GetPointIds()->SetId(k,ptID);
                        
                        ptID++;
                    }
                    
                    if(burgers.squaredNorm()>FLT_EPSILON)
                    {
                        cells->InsertNextCell(line);
                    }
                    else
                    {
                        cells0->InsertNextCell(line);
                    }
                    
                }
            }
            
            polyData->SetPoints(points);
            polyData->SetLines(cells);
            polyData->GetCellData()->SetScalars(colors);
            
            polyData0->SetPoints(points);
            polyData0->SetLines(cells0);
            //            polyData0->GetCellData()->SetScalars(colors);
            
        }
        
        /**********************************************************************/
        DislocationSegmentActor(const size_t& frameID,vtkRenderer* const ren) :
        /* init */ renderer(ren),
        //        /* init */ vertexReader(vertexReader_in),
        //        /* init */ edgeReader(edgeReader_in),
        /* init */ lineActor(vtkSmartPointer<vtkActor>::New()),
        /* init */ tubeActor(vtkSmartPointer<vtkActor>::New()),
        /* init */ lineActor0(vtkSmartPointer<vtkActor>::New()),
        /* init */ tubeActor0(vtkSmartPointer<vtkActor>::New()),
        
        //        /* init */ ptID(0),
        /* init */ points(vtkSmartPointer<vtkPoints>::New()),
        /* init */ cells(vtkSmartPointer<vtkCellArray>::New()),
        /* init */ polyData(vtkSmartPointer<vtkPolyData>::New()),
        /* init */ cells0(vtkSmartPointer<vtkCellArray>::New()),
        /* init */ polyData0(vtkSmartPointer<vtkPolyData>::New()),
        /* init */ colors(vtkSmartPointer<vtkUnsignedCharArray>::New()),
        /* init */ lineMapper(vtkSmartPointer<vtkPolyDataMapper>::New()),
        /* init */ tubeFilter(vtkSmartPointer<vtkTubeFilter>::New()),
        /* init */ tubeMapper(vtkSmartPointer<vtkPolyDataMapper>::New()),
        /* init */ lineMapper0(vtkSmartPointer<vtkPolyDataMapper>::New()),
        /* init */ tubeFilter0(vtkSmartPointer<vtkTubeFilter>::New()),
        /* init */ tubeMapper0(vtkSmartPointer<vtkPolyDataMapper>::New()),
        /* init */ nodePoints(vtkSmartPointer<vtkPoints>::New()),
        //        /* init */ nodeLabels(vtkSmartPointer<vtkStringArray>::New()),
        /* init */ sphereSource(vtkSmartPointer<vtkSphereSource>::New()),
        /* init */ nodeData(vtkSmartPointer<vtkPolyData>::New()),
        /* init */ nodeGlyphs(vtkSmartPointer<vtkGlyph3D>::New()),
        /* init */ nodeMapper(vtkSmartPointer<vtkPolyDataMapper>::New()),
        /* init */ nodeActor(vtkSmartPointer<vtkActor>::New()),
        //        /* init */ velocityPoints(vtkSmartPointer<vtkPoints>::New()),
        /* init */ velocityVectors(vtkSmartPointer<vtkDoubleArray>::New()),
        /* init */ velocityColors(vtkSmartPointer<vtkUnsignedCharArray>::New()),
        /* init */ velocityPolyData(vtkSmartPointer<vtkPolyData>::New()),
        /* init */ velocityArrowSource(vtkSmartPointer<vtkArrowSource>::New()),
        /* init */ velocityGlyphs(vtkSmartPointer<vtkGlyph3D>::New()),
        /* init */ velocityMapper(vtkSmartPointer<vtkPolyDataMapper>::New()),
        /* init */ velocityActor(vtkSmartPointer<vtkActor>::New()),
        /* init */ labelPolyData(vtkSmartPointer<vtkPolyData>::New()),
        /* init */ labelScalars(vtkSmartPointer<vtkDoubleArray>::New()),
        /* init */ labelMapper(vtkSmartPointer<vtkLabeledDataMapper>::New()),
        /* init */ labelActor(vtkSmartPointer<vtkActor2D>::New())
        {
            
            colors->SetNumberOfComponents(3);
            velocityVectors->SetNumberOfComponents(3);
            velocityVectors->SetName("nodeVelocity");
            velocityColors->SetNumberOfComponents(3);
            
            labelScalars->SetNumberOfComponents(1);
            labelScalars->SetName("node IDs");
            
            
            //            nodeLabels->SetName("node IDs");
            
            
            readNodes(frameID);
            readSegments(frameID);
            
            // Populate polyData
            
            // tube filter
            tubeFilter->SetInputData(polyData);
            tubeFilter->SetRadius(tubeRadius); // this must be a function similar to setColor
            tubeFilter->SetNumberOfSides(10);
            tubeFilter->Update();
            tubeFilter0->SetInputData(polyData0);
            tubeFilter0->SetRadius(tubeRadius); // this must be a function similar to setColor
            tubeFilter0->SetNumberOfSides(10);
            tubeFilter0->Update();
            
            // Mappers
            tubeMapper->SetInputConnection(tubeFilter->GetOutputPort());
            tubeMapper->ScalarVisibilityOn();
            lineMapper->SetInputData(polyData);
            tubeMapper0->SetInputConnection(tubeFilter0->GetOutputPort());
            tubeMapper0->ScalarVisibilityOn();
            lineMapper0->SetInputData(polyData0);
            
            // Actors
            tubeActor->SetMapper(tubeMapper);
            tubeActor0->SetMapper(tubeMapper0);
            tubeActor0->GetProperty()->SetColor(0.5, 0.5, 0.5); //(R,G,B)
            tubeActor0->GetProperty()->SetOpacity(0.3); //(R,G,B)
            
            
            //            tubeActor->GetProperty()->SetOpacity(1.0); //Make the tube have some transparency.
            //            computeColor();
            //tube->GetProperty()->SetColor(colorVector(0),colorVector(1),colorVector(2)); // Give some color to the tube
            lineActor->SetMapper(lineMapper);
            lineActor0->SetMapper(lineMapper0);
            
            // Add actors to renderer
            renderer->AddActor(tubeActor);
            renderer->AddActor(tubeActor0);
            
            nodeGlyphs->SetSourceConnection(sphereSource->GetOutputPort());
            nodeGlyphs->SetInputData(nodeData);
            nodeGlyphs->ScalingOn();
            //            nodeGlyphs->SetScaleModeToScaleByVector();
            nodeGlyphs->SetScaleFactor(2.0*tubeRadius*1.2);
            nodeGlyphs->OrientOn();
            nodeGlyphs->ClampingOff();
            nodeGlyphs->SetVectorModeToUseVector();
            nodeGlyphs->SetIndexModeToOff();
            
            nodeMapper->SetInputConnection(nodeGlyphs->GetOutputPort());
            nodeMapper->ScalarVisibilityOff();
            
            // Set up actor
            nodeActor->SetMapper(nodeMapper);
            
            // Add actor to renderer
            renderer->AddActor(nodeActor);
            
            
            velocityGlyphs->SetSourceConnection(velocityArrowSource->GetOutputPort());
            velocityGlyphs->SetInputData(velocityPolyData);
            velocityGlyphs->ScalingOn();
            velocityGlyphs->SetScaleModeToScaleByVector();
            //            velocityGlyphs->SetColorModeToColorByScale();
            //            velocityGlyphs->SetColorModeToColorByScalar();
            velocityGlyphs->SetColorModeToColorByVector();
            //            velocityGlyphs->SetScaleFactor(velocityFactor);
            velocityGlyphs->OrientOn();
            velocityGlyphs->ClampingOff();
            velocityGlyphs->SetVectorModeToUseVector();
            velocityGlyphs->SetIndexModeToOff();
            
            // Velocities
            velocityMapper->SetInputConnection(velocityGlyphs->GetOutputPort());
            velocityMapper->ScalarVisibilityOff();
            velocityActor->SetMapper(velocityMapper);
            velocityActor->GetProperty()->SetColor(1.0, 0.0, 1.0); //(R,G,B)
            renderer->AddActor(velocityActor);
            
            // Labels
            labelMapper->SetInputData(labelPolyData);
            labelMapper->SetLabelModeToLabelScalars();
            labelMapper->SetLabelFormat("%1.0f");
            labelActor->SetMapper(labelMapper);
            labelActor->GetProperty()->SetColor(0.0, 0.0, 0.0); //(R,G,B)
            renderer->AddActor(labelActor);
            
            modify();
        }
        
        /**********************************************************************/
        ~DislocationSegmentActor()
        {
            renderer->RemoveActor(tubeActor);
            renderer->RemoveActor(tubeActor0);
            renderer->RemoveActor(nodeActor);
            renderer->RemoveActor(velocityActor);
            renderer->RemoveActor(labelActor);
        }
        
        /**********************************************************************/
        void modify()
        {
            //            line->GetProperty()->SetColor(0.0,1.0,0.0); // Give some color to the tube
            //            line->Modified();
            
            //            tubeFilter->SetRadius(tubeRadius); // this must be a function similar to setColor
            //            computeColor();
            //            tubeActor->GetProperty()->SetColor(colorVector(0),colorVector(1),colorVector(2)); // Give some color to the tube
            //            tube->Modified();
            
            tubeFilter->SetRadius(tubeRadius); // this must be a function similar to setColor
            tubeFilter0->SetRadius(tubeRadius); // this must be a function similar to setColor
            
            if(showZeroBuergers)
            {
                tubeActor0->VisibilityOn();
                
            }
            else
            {
                tubeActor0->VisibilityOff();
            }
            
            nodeGlyphs->SetScaleFactor(2.0*tubeRadius*1.2);
            
            if(showVelocities)
            {
                velocityActor->VisibilityOn();
                
            }
            else
            {
                velocityActor->VisibilityOff();
            }
            
            if(showNodeIDs)
            {
                labelActor->VisibilityOn();
                
            }
            else
            {
                labelActor->VisibilityOff();
            }
            
            velocityGlyphs->SetScaleFactor(velocityFactor);
            
            
        }
        
    };
    
    // Static data members
    float DislocationSegmentActor::tubeRadius=5.0;
    float DislocationSegmentActor::alpha=0.5;
    DislocationSegmentActor::ColorScheme DislocationSegmentActor::clr=DislocationSegmentActor::colorBurgers;
    size_t DislocationSegmentActor::Np=10;
    bool DislocationSegmentActor::plotBoundarySegments=true;
    bool DislocationSegmentActor::showVelocities=false;
    bool DislocationSegmentActor::showNodeIDs=false;
    float DislocationSegmentActor::velocityFactor=100.0;
    bool DislocationSegmentActor::showZeroBuergers=false;
    
    
} // namespace model
#endif






