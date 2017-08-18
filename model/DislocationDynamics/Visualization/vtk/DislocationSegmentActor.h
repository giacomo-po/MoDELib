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
#include <vtkGlyph3D.h>

//#include <model/IO/EdgeReader.h>
//#include <model/IO/vertexReader.h>
#include <model/IO/IDreader.h>

// VTK documentation
// http://www.vtk.org/Wiki/VTK/Examples/Cxx/GeometricObjects/PolyLine
// http://www.vtk.org/Wiki/VTK/Examples/Cxx/VisualizationAlgorithms/TubesWithVaryingRadiusAndColors

namespace model
{
    struct DislocationSegmentActor :
    /* inherits from   */ public IDreader<'V',1,6, double>,
    /* inherits from   */ public IDreader<'K',2,13,double>
    {
        
//    public:

        static constexpr int dim=3;
        enum ColorScheme {colorBurgers=0,colorSessile=1,colorNormal=2,colorEdgeScrew=3,colorComponent=4};
//        typedef VertexReader<'V',10,double> VertexReaderType;
//        typedef EdgeReader<'K',15,double>   EdgeReaderType;
        typedef IDreader<'V',1,6, double> VertexReaderType;
        typedef IDreader<'K',2,13,double> EdgeReaderType;

        typedef Eigen::Matrix<float,dim,1>  VectorDim;
        
        static float alpha;
        static float tubeRadius;
        static ColorScheme clr;
        static size_t Np;      // No. of vertices per line
        static bool plotBoundarySegments;
        
        vtkSmartPointer<vtkActor> lineActor;
        vtkSmartPointer<vtkActor> tubeActor;
        
        
//    private:
        
        VectorDim planeNormal;
        VectorDim burgers;
        VectorDim chord;
        VectorDim colorVector;
        
        
        // segments objects
        vtkSmartPointer<vtkPoints> points;
//        std::deque<vtkSmartPointer<vtkPolyLine>> lines;
        vtkSmartPointer<vtkCellArray> cells;
        vtkSmartPointer<vtkPolyData> polyData;
        vtkSmartPointer<vtkUnsignedCharArray> colors;
        vtkSmartPointer<vtkPolyDataMapper> lineMapper;
        vtkSmartPointer<vtkTubeFilter> tubeFilter;
        vtkSmartPointer<vtkPolyDataMapper> tubeMapper;
        
        // node objects
        vtkSmartPointer<vtkPoints> nodePoints;
        vtkSmartPointer<vtkSphereSource> sphereSource;
        vtkSmartPointer<vtkPolyData> nodeData;
        vtkSmartPointer<vtkGlyph3D> nodeGlyphs;
        vtkSmartPointer<vtkPolyDataMapper> nodeMapper;
        vtkSmartPointer<vtkActor> nodeActor;

        
        
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
            
            for(const auto& node : vertexReader())
            {
                Eigen::Map<const Eigen::Matrix<double,1,6>> row(node.second.data());
                
                nodePoints->InsertNextPoint(row.template segment<dim>(0).data());

            }
         
            nodeData->SetPoints(nodePoints);
            //nodeData->GetPointData()->SetVectors(vectors);
            nodeData->Modified();

            
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
                    unsigned char clr[3]={0,255,255};
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
                    
                    cells->InsertNextCell(line);
                    
                }
            }
            
            polyData->SetPoints(points);
            polyData->SetLines(cells);
            polyData->GetCellData()->SetScalars(colors);
        }
        
        /**********************************************************************/
        DislocationSegmentActor(const size_t& frameID,vtkRenderer* renderer) :
//        /* init */ vertexReader(vertexReader_in),
//        /* init */ edgeReader(edgeReader_in),
        /* init */ lineActor(vtkSmartPointer<vtkActor>::New()),
        /* init */ tubeActor(vtkSmartPointer<vtkActor>::New()),
        //        /* init */ ptID(0),
        /* init */ points(vtkSmartPointer<vtkPoints>::New()),
        /* init */ cells(vtkSmartPointer<vtkCellArray>::New()),
        /* init */ polyData(vtkSmartPointer<vtkPolyData>::New()),
        /* init */ colors(vtkSmartPointer<vtkUnsignedCharArray>::New()),
        /* init */ lineMapper(vtkSmartPointer<vtkPolyDataMapper>::New()),
        /* init */ tubeFilter(vtkSmartPointer<vtkTubeFilter>::New()),
        /* init */ tubeMapper(vtkSmartPointer<vtkPolyDataMapper>::New()),
        /* init */ nodePoints(vtkSmartPointer<vtkPoints>::New()),
        /* init */ sphereSource(vtkSmartPointer<vtkSphereSource>::New()),
        /* init */ nodeData(vtkSmartPointer<vtkPolyData>::New()),
        /* init */ nodeGlyphs(vtkSmartPointer<vtkGlyph3D>::New()),
        /* init */ nodeMapper(vtkSmartPointer<vtkPolyDataMapper>::New()),
        /* init */ nodeActor(vtkSmartPointer<vtkActor>::New())

        {
            
            colors->SetNumberOfComponents(3);

            readNodes(frameID);
            readSegments(frameID);
            
            // Populate polyData
            
            // tube filter
            tubeFilter->SetInputData(polyData);
            tubeFilter->SetRadius(tubeRadius); // this must be a function similar to setColor
            tubeFilter->SetNumberOfSides(10);
            tubeFilter->Update();
            
            // Mappers
            tubeMapper->SetInputConnection(tubeFilter->GetOutputPort());
            tubeMapper->ScalarVisibilityOn();
            lineMapper->SetInputData(polyData);
            
            // Actors
            tubeActor->SetMapper(tubeMapper);
//            tubeActor->GetProperty()->SetOpacity(1.0); //Make the tube have some transparency.
            //            computeColor();
            //tube->GetProperty()->SetColor(colorVector(0),colorVector(1),colorVector(2)); // Give some color to the tube
            lineActor->SetMapper(lineMapper);
            
            // Add actors to renderer
            renderer->AddActor(tubeActor);
            
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
            
        }
        
    };
    
    // Static data members
    float DislocationSegmentActor::tubeRadius=5.0;
    float DislocationSegmentActor::alpha=0.5;
    DislocationSegmentActor::ColorScheme DislocationSegmentActor::clr=DislocationSegmentActor::colorBurgers;
    size_t DislocationSegmentActor::Np=10;
    bool DislocationSegmentActor::plotBoundarySegments=true;
} // namespace model
#endif







