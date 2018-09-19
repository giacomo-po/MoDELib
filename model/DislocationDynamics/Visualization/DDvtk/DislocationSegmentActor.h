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

//#include <model/IO/EdgeReader.h>
//#include <model/IO/vertexReader.h>
#include <model/IO/IDreader.h>
#include <model/Geometry/PlanarPolygon.h>
#include <model/DislocationDynamics/IO/EVLio.h>

// VTK documentation
// http://www.vtk.org/Wiki/VTK/Examples/Cxx/GeometricObjects/PolyLine
// http://www.vtk.org/Wiki/VTK/Examples/Cxx/VisualizationAlgorithms/TubesWithVaryingRadiusAndColors

namespace model
{
    struct DislocationSegmentActor : public EVLio<3>,
    //    public std::map<size_t,const DislocationNodeIO<3>* const>,
    //        public std::map<size_t,const DislocationLoopIO<3>* const>,
    public std::map<std::pair<size_t,size_t>,std::set<const DislocationEdgeIO<3>*>>
    //    :
    //    /* inherits from   */ public IDreader<'V',1,10, double>,
    //    /* inherits from   */ public IDreader<'K',2,13,double>,
    //    /* inherits from   */ IDreader<'E',3,0,double>,
    //    /* inherits from   */ IDreader<'L',1,13,double>
    {
        
        //    public:
        
        static constexpr int dim=3;
        enum ColorScheme {colorBurgers=0,colorSessile=1,colorNormal=2,colorEdgeScrew=3,colorComponent=4};
        //        typedef VertexReader<'V',10,double> VertexReaderType;
        //        typedef EdgeReader<'K',15,double>   EdgeReaderType;
        typedef IDreader<'V',1,10, double> VertexReaderType;
        typedef IDreader<'K',2,13,double> EdgeReaderType;
        typedef IDreader<'E',3,0,double> LoopLinkReaderType;
        typedef IDreader<'L',1,13,double> LoopReaderType;
        
        
        typedef Eigen::Matrix<float,dim,1>  VectorDim;
        
        static float alpha;
        static float tubeRadius;
        static bool scaleRadiusByBurgers;
        static ColorScheme clr;
        static size_t Np;      // No. of vertices per line
        static bool showBoundarySegments;
        static bool blackGrainBoundarySegments;
        static bool showVelocities;
        static bool showNodeIDs;
        static float velocityFactor;
        static bool showZeroBuergers;
        static bool showSingleNode;
        static size_t singleNodeID;
        static bool showNodes;
        static bool showSlippedArea;
        static float slippedAreaOpacity;
        
        vtkRenderer* const renderer;
        
        vtkSmartPointer<vtkActor> lineActor;
        vtkSmartPointer<vtkActor> tubeActor;
        vtkSmartPointer<vtkActor> tubeActorBnd;
        
        vtkSmartPointer<vtkActor> lineActor0;
        vtkSmartPointer<vtkActor> tubeActor0;
        
        
        //    private:
        
        //        VectorDim planeNormal;
        //        VectorDim burgers;
        //        VectorDim chord;
        //        Eigen::Matrix<int,dim,1> colorVector;
        
        
        // segments objects
        vtkSmartPointer<vtkPoints> points;
        //        std::deque<vtkSmartPointer<vtkPolyLine>> lines;
        vtkSmartPointer<vtkCellArray> cells;
        vtkSmartPointer<vtkCellArray> cellsBnd;
        vtkSmartPointer<vtkCellArray> cells0;
        vtkSmartPointer<vtkPolyData> polyData;
        vtkSmartPointer<vtkPolyData> polyDataBnd;
        vtkSmartPointer<vtkPolyData> polyData0;
        vtkSmartPointer<vtkFloatArray> radii;
        vtkSmartPointer<vtkUnsignedCharArray> colors;
        vtkSmartPointer<vtkUnsignedCharArray> colorsBnd;
        vtkSmartPointer<vtkPolyDataMapper> lineMapper;
        //        vtkSmartPointer<vtkPolyDataMapper> lineMapperBnd;
        vtkSmartPointer<vtkPolyDataMapper> lineMapper0;
        vtkSmartPointer<vtkTubeFilter> tubeFilter;
        vtkSmartPointer<vtkTubeFilter> tubeFilterBnd;
        vtkSmartPointer<vtkTubeFilter> tubeFilter0;
        vtkSmartPointer<vtkPolyDataMapper> tubeMapper;
        vtkSmartPointer<vtkPolyDataMapper> tubeMapperBnd;
        vtkSmartPointer<vtkPolyDataMapper> tubeMapper0;
        
        // node objects
        vtkSmartPointer<vtkPoints> nodePoints;
        //        vtkSmartPointer<vtkStringArray> nodeLabels;
        vtkSmartPointer<vtkSphereSource> sphereSource;
        vtkSmartPointer<vtkPolyData> nodeData;
        vtkSmartPointer<vtkUnsignedCharArray> nodeColors;
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
        
        
        vtkSmartPointer<vtkPoints> singleNodePoint;
        vtkSmartPointer<vtkPolyData> singleNodeLabelPolyData;
        vtkSmartPointer<vtkDoubleArray> singleNodeLabelScalars;
        vtkSmartPointer<vtkLabeledDataMapper> singleNodeLabelMapper;
        vtkSmartPointer<vtkActor2D> singleNodeLabelActor;
        
        
        // Slipped Area
        vtkSmartPointer<vtkCellArray> triangles;
        vtkSmartPointer<vtkPolyData> trianglePolyData;
        vtkSmartPointer<vtkPolyDataMapper> triangleMapper;
        vtkSmartPointer<vtkActor> triangleActor;
        
        
        //        const std::map<size_t,const DislocationNodeIO<dim>* const>& this->nodeMap() const
        //        {
        //            return *this;
        //        }
        //
        //        std::map<size_t,const DislocationNodeIO<dim>* const>& this->nodeMap()
        //        {
        //            return *this;
        //        }
        //
        //        const std::map<size_t,const DislocationLoopIO<dim>* const>& loopMap() const
        //        {
        //            return *this;
        //        }
        //
        //        std::map<size_t,const DislocationLoopIO<dim>* const>& loopMap()
        //        {
        //            return *this;
        //        }
        
        const std::map<std::pair<size_t,size_t>,std::set<const DislocationEdgeIO<3>*>>& edgeMap() const
        {
            return *this;
        }
        
        std::map<std::pair<size_t,size_t>,std::set<const DislocationEdgeIO<3>*>>& edgeMap()
        {
            return *this;
        }
        
        
        /*********************************************************************/
        Eigen::Matrix<int,3,1> computeColor(const VectorDim& burgers, const VectorDim& chord, const VectorDim& planeNormal)
        {
            
            VectorDim clrVector(VectorDim::Zero());
            
            switch (clr)
            {

                    
                case colorSessile:
                {
                    const bool isSessile=fabs(burgers.dot(planeNormal))>FLT_EPSILON;
                    clrVector(0)= isSessile? 1.0 : 0.1;
                    clrVector(1)= isSessile? 0.5 : 0.4;
                    clrVector(2)= isSessile? 0.0 : 0.9;
                    break;
                }
                    
                    //                case colorComponent:
                    //                {
                    //                    RGBcolor rgb(RGBmap::getColor(ids,sIDmin,sIDmax));
                    //                    colorVector << rgb.r, rgb.g, rgb.b;
                    //                }
                    //                    break;
                    
//                                    case colorEdgeScrew:
//                                    {
//                                        const float u = std::fabs(tubeTangents.col(k).normalized().dot(burgers.normalized()));
//                                        //                            RGBcolor rgb(RGBmap::getColor(std::fabs(tubeTangents.col(k).normalized().dot(burgers.normalized())),0,1));
//                                        //                            colorVector << rgb.r, rgb.g, rgb.b;
//                                        colorVector=screwColor*u+edgeColor*(1-u);
//                                    }
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
            
            //			VectorDim colorVector = burgers + VectorDim::Ones(dim) * burgers.norm();
            clrVector = (clrVector + VectorDim::Ones(dim) * clrVector.norm()).eval();
            
            //		colorVector << 0.0f,0.6f,0.4f;
            clrVector.normalize();
            
//            colorVector=(clrVector*255).cast<int>();
            
            return (clrVector*255).cast<int>();
        }
        
        
        
        //    public:
        
        //        /**********************************************************************/
        //        VertexReaderType& vertexReader()
        //        {
        //            return *this;
        //        }
        //
        //        /**********************************************************************/
        //        EdgeReaderType& edgeReader()
        //        {
        //            return *this;
        //        }
        //
        //        /**********************************************************************/
        //        LoopLinkReaderType& loopLinkReader()
        //        {
        //            return *this;
        //        }
        //
        //        LoopReaderType& loopReader()
        //        {
        //            return *this;
        //        }
        
        
        /**********************************************************************/
        void createNodes()
        {
            //            if (vertexReader().isGood(frameID,false)) // bin format
            //            {
            //                vertexReader().read(frameID,false);
            //            }
            //            else // txt format
            //            {
            //                vertexReader().read(frameID,true);
            //            }
            
            //            nodeLabels->SetNumberOfValues(vertexReader().size());
            //            size_t labelID=0;
            //            for(const auto& node : this->nodes())
            //            {
            //                this->nodeMap().emplace(node.sID,&node);
            //            }
            
            for(const auto& node : this->nodeMap())
            {
                
                //                this->nodeMap().emplace(node.sID,&node);
                
                //                Eigen::Map<const Eigen::Matrix<double,1,10>> row(node.second.data());
                //                std::cout<<node.second->P.transpose()<<std::endl;
                nodePoints->InsertNextPoint(node.second->P.data());
                
                //                const int& meshLocation(node.(8));
                switch (node.second->meshLocation)
                {
                    case 0:
                    {
                        unsigned char nodeClr[3]={100,100,100};
                        nodeColors->InsertNextTypedTuple(nodeClr);
                        break;
                    }
                        
                    case 1:
                    {
                        unsigned char nodeClr[3]={0,255,255};
                        nodeColors->InsertNextTypedTuple(nodeClr);
                        break;
                    }
                        
                    case 2:
                    {
                        unsigned char nodeClr[3]={255,0,255};
                        nodeColors->InsertNextTypedTuple(nodeClr);
                        break;
                    }
                        
                    default:
                    {
                        unsigned char nodeClr[3]={1,1,1};
                        nodeColors->InsertNextTypedTuple(nodeClr);
                        break;
                    }
                }
                
                // Velocity
                velocityVectors->InsertNextTuple(node.second->V.data()); // arrow vactor
                unsigned char velClr[3]={255,0,255};
                velocityColors->InsertNextTypedTuple(velClr);
                
                
                
                // Labels
                labelScalars->InsertNextTuple1(node.second->sID);
                
                // Single node
                if(node.second->sID==singleNodeID)
                {
                    singleNodePoint->InsertNextPoint(node.second->P.data());
                    singleNodeLabelScalars->InsertNextTuple1(node.second->sID);
                }
                
            }
            
            nodeData->SetPoints(nodePoints);
            nodeData->GetPointData()->SetScalars(nodeColors);
            nodeData->Modified();
            
            velocityPolyData->SetPoints(nodePoints);
            velocityPolyData->GetPointData()->SetVectors(velocityVectors);
            //            velocityPolyData->Modified();
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
            
            singleNodeLabelPolyData->SetPoints(singleNodePoint);
            singleNodeLabelPolyData->GetPointData()->SetScalars(singleNodeLabelScalars);
        }
        

        /**********************************************************************/
        void createSegments()
        {
            
            
                        std::map<std::pair<size_t,size_t>,DislocationSegmentIO<dim>> segments=this->segments();
            
            size_t ptID=0;
            for (const auto& segment : segments)
            {
//                assert(edge.first.second.size());
                
                auto itSource(this->nodeMap().find(segment.second.sourceID)); //source
                assert(itSource!=this->nodeMap().end() && "SOURCE VERTEX NOT FOUND IN V-FILE");
                auto   itSink(this->nodeMap().find(segment.second.sinkID)); //sink
                assert(  itSink!=this->nodeMap().end() &&   "SINK VERTEX NOT FOUND IN V-FILE");
                
                //                Eigen::Map<const Eigen::Matrix<double,1,6>> sourceRow(itSource->second.data());
                //                Eigen::Map<const Eigen::Matrix<double,1,6>>   sinkRow(  itSink->second.data());
                //                Eigen::Map<const Eigen::Matrix<double,1,13>>   edgeRow(edge.second.data());
                
                //                const int   snID(edgeRow(2*dim+2));
                //                const bool sourceOnBoundary(sourceRow(2*dim+1));
                //                const bool   sinkOnBoundary(  sinkRow(2*dim+1));
                //
                //                if(!(sourceOnBoundary && sinkOnBoundary) || showBoundarySegments)
                //                {
                
//                int meshLocation=(*edge.second.begin())->meshLocation;
                
                Eigen::Matrix<float,dim,6> P0T0P1T1BN;
                
                //                P0T0P1T1BN.col(0) = sourceRow.segment<dim>(0*dim).transpose().template cast<float>();	// source position
                //                P0T0P1T1BN.col(2) =   sinkRow.segment<dim>(0*dim).transpose().template cast<float>();	// sink position
                //                //                    P0T0P1T1BN.col(1) = sourceTfactor*(itSource->second.segment<dim>(1*dim).transpose().template cast<float>());	// source tangent
                //                //                    P0T0P1T1BN.col(3) =  -sinkTfactor*(  itSink->second.segment<dim>(1*dim).transpose().template cast<float>());	// sink tangent
                //                P0T0P1T1BN.col(1) = edgeRow.segment<dim>(2*dim).transpose().template cast<float>();	// source tangent
                //                P0T0P1T1BN.col(3) = edgeRow.segment<dim>(3*dim).transpose().template cast<float>();	// sink tangent
                //                P0T0P1T1BN.col(4) = edgeRow.segment<dim>(0*dim).transpose().template cast<float>();		// Burgers vector
                //                P0T0P1T1BN.col(5) = edgeRow.segment<dim>(1*dim).transpose().template cast<float>();		// plane normal
                P0T0P1T1BN.col(0) = itSource->second->P.template cast<float>();	// source position
                P0T0P1T1BN.col(2) = itSink->second->P.template cast<float>();	// sink position
                //                    P0T0P1T1BN.col(1) = sourceTfactor*(itSource->second.segment<dim>(1*dim).transpose().template cast<float>());	// source tangent
                //                    P0T0P1T1BN.col(3) =  -sinkTfactor*(  itSink->second.segment<dim>(1*dim).transpose().template cast<float>());	// sink tangent
                //                P0T0P1T1BN.col(1) = edgeRow.segment<dim>(2*dim).transpose().template cast<float>();	// source tangent
                //                P0T0P1T1BN.col(3) = edgeRow.segment<dim>(3*dim).transpose().template cast<float>();	// sink tangent
                P0T0P1T1BN.col(4)=segment.second.b.template cast<float>();		// Burgers vector
//                for(auto& loopLink : edge.second)
//                {
//                    auto loopIter(this->loopMap().find(loopLink->loopID));
//                    assert(loopIter!=this->loopMap().end());
//                    if(loopLink->sourceID==edge.first.first && loopLink->sinkID==edge.first.second)
//                    {
//                        P0T0P1T1BN.col(4)+=loopIter->second->B.template cast<float>();
//                    }
//                    else if(loopLink->sourceID==edge.first.second && loopLink->sinkID==edge.first.first)
//                    {
//                        P0T0P1T1BN.col(4)-=loopIter->second->B.template cast<float>();
//                    }
//                    else
//                    {
//                        assert(0);
//                    }
//                }
                //                P0T0P1T1BN.col(5) = edgeRow.segment<dim>(1*dim).transpose().template cast<float>();		// plane normal
                P0T0P1T1BN.col(5)=segment.second.n.template cast<float>();		// plane normal
                
                VectorDim chord = P0T0P1T1BN.col(2)-P0T0P1T1BN.col(0);
                const VectorDim burgers=segment.second.b.template cast<float>();
                const VectorDim planeNormal=segment.second.n.template cast<float>();
                
                const float g = std::pow(chord.norm(),alpha);
                
                
                //                    lines.push_back(vtkSmartPointer<vtkPolyLine>::New());
                //                    auto& line(*lines.rbegin());
                vtkSmartPointer<vtkPolyLine> line=vtkSmartPointer<vtkPolyLine>::New();
                line->GetPointIds()->SetNumberOfIds(Np);
                
                //                    unsigned char clr0[3]={255,255,255};
                
                const float burgersNorm(burgers.norm());
                
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
                    radii->InsertNextValue(burgersNorm*tubeRadius);
                    line->GetPointIds()->SetId(k,ptID);
                    
                    ptID++;
                }
                
                
                if(burgers.squaredNorm()>FLT_EPSILON)
                {
                    
                    Eigen::Matrix<int,3,1> colorVector=computeColor(burgers,chord,planeNormal);
                    
                    //                    unsigned char lineClr[3]={51,153,255};
                    unsigned char lineClr[3]={(unsigned char) colorVector(0),(unsigned char) colorVector(1),(unsigned char) colorVector(2)};
                    
                    if(segment.second.meshLocation>=1)
                    {
                        cellsBnd->InsertNextCell(line);
                        colorsBnd->InsertNextTypedTuple(lineClr);
                        
                    }
                    else
                    {
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
                    }
                    
                    //                        if(meshLocation==1)
                    //                        {
                    //                            cellsBnd->InsertNextCell(line);
                    //                            colorsBnd->InsertNextTypedTuple(lineClr);
                    //
                    //                        }
                    //                        else
                    //                        {
                    //                            cells->InsertNextCell(line);
                    //
                    //                            if(meshLocation==2 && blackGrainBoundarySegments)
                    //                            {
                    //                                unsigned char lineClr1[3]={1,1,1};
                    //                                colors->InsertNextTypedTuple(lineClr1);
                    //                            }
                    //                            else
                    //                            {
                    //                                colors->InsertNextTypedTuple(lineClr);
                    //                            }
                    //                        }
                }
                else
                {
                    cells0->InsertNextCell(line);
                }
                
                //                }
            }
            
            polyData->SetPoints(points);
            polyData->SetLines(cells);
            //            polyData->GetCellData()->SetScalars(colors);
            
            //            polyData->GetPointData()->SetScalars(radii);
            //            polyData->GetCellData()->SetScalars(radii);
            polyData->GetCellData()->AddArray(colors);
            polyData->GetPointData()->AddArray(radii);
            polyData->GetPointData()->SetActiveScalars("TubeRadius");
            
            
            polyDataBnd->SetPoints(points);
            polyDataBnd->SetLines(cellsBnd);
            polyDataBnd->GetCellData()->SetScalars(colorsBnd);
            
            
            polyData0->SetPoints(points);
            polyData0->SetLines(cells0);
            //            polyData0->GetCellData()->SetScalars(colors);
            
        }
        
//        /**********************************************************************/
//        void createSegments()
//        {
//            //            if (edgeReader().isGood(frameID,false)) // bin format
//            //            {
//            //                edgeReader().read(frameID,false);
//            //            }
//            //            else // txt format
//            //            {
//            //                edgeReader().read(frameID,true);
//            //            }
//            
//            
////            std::map<std::pair<size_t,size_t>,DislocationSegmentIO<dim>> segments=this->segments();
//            
//            size_t ptID=0;
//            for (const auto& edge : edgeMap())
//            {
//                assert(edge.first.second.size());
//                
//                auto itSource(this->nodeMap().find(edge.first.first)); //source
//                assert(itSource!=this->nodeMap().end() && "SOURCE VERTEX NOT FOUND IN V-FILE");
//                auto   itSink(this->nodeMap().find(edge.first.second)); //sink
//                assert(  itSink!=this->nodeMap().end() &&   "SINK VERTEX NOT FOUND IN V-FILE");
//                
//                //                Eigen::Map<const Eigen::Matrix<double,1,6>> sourceRow(itSource->second.data());
//                //                Eigen::Map<const Eigen::Matrix<double,1,6>>   sinkRow(  itSink->second.data());
//                //                Eigen::Map<const Eigen::Matrix<double,1,13>>   edgeRow(edge.second.data());
//                
//                //                const int   snID(edgeRow(2*dim+2));
//                //                const bool sourceOnBoundary(sourceRow(2*dim+1));
//                //                const bool   sinkOnBoundary(  sinkRow(2*dim+1));
//                //
//                //                if(!(sourceOnBoundary && sinkOnBoundary) || showBoundarySegments)
//                //                {
//                
//                int meshLocation=(*edge.second.begin())->meshLocation;
//                
//                Eigen::Matrix<float,dim,6> P0T0P1T1BN;
//                
//                //                P0T0P1T1BN.col(0) = sourceRow.segment<dim>(0*dim).transpose().template cast<float>();	// source position
//                //                P0T0P1T1BN.col(2) =   sinkRow.segment<dim>(0*dim).transpose().template cast<float>();	// sink position
//                //                //                    P0T0P1T1BN.col(1) = sourceTfactor*(itSource->second.segment<dim>(1*dim).transpose().template cast<float>());	// source tangent
//                //                //                    P0T0P1T1BN.col(3) =  -sinkTfactor*(  itSink->second.segment<dim>(1*dim).transpose().template cast<float>());	// sink tangent
//                //                P0T0P1T1BN.col(1) = edgeRow.segment<dim>(2*dim).transpose().template cast<float>();	// source tangent
//                //                P0T0P1T1BN.col(3) = edgeRow.segment<dim>(3*dim).transpose().template cast<float>();	// sink tangent
//                //                P0T0P1T1BN.col(4) = edgeRow.segment<dim>(0*dim).transpose().template cast<float>();		// Burgers vector
//                //                P0T0P1T1BN.col(5) = edgeRow.segment<dim>(1*dim).transpose().template cast<float>();		// plane normal
//                P0T0P1T1BN.col(0) = itSource->second->P.template cast<float>();	// source position
//                P0T0P1T1BN.col(2) = itSink->second->P.template cast<float>();	// sink position
//                //                    P0T0P1T1BN.col(1) = sourceTfactor*(itSource->second.segment<dim>(1*dim).transpose().template cast<float>());	// source tangent
//                //                    P0T0P1T1BN.col(3) =  -sinkTfactor*(  itSink->second.segment<dim>(1*dim).transpose().template cast<float>());	// sink tangent
//                //                P0T0P1T1BN.col(1) = edgeRow.segment<dim>(2*dim).transpose().template cast<float>();	// source tangent
//                //                P0T0P1T1BN.col(3) = edgeRow.segment<dim>(3*dim).transpose().template cast<float>();	// sink tangent
//                P0T0P1T1BN.col(4).setZero();		// Burgers vector
//                for(auto& loopLink : edge.second)
//                {
//                    auto loopIter(this->loopMap().find(loopLink->loopID));
//                    assert(loopIter!=this->loopMap().end());
//                    if(loopLink->sourceID==edge.first.first && loopLink->sinkID==edge.first.second)
//                    {
//                        P0T0P1T1BN.col(4)+=loopIter->second->B.template cast<float>();
//                    }
//                    else if(loopLink->sourceID==edge.first.second && loopLink->sinkID==edge.first.first)
//                    {
//                        P0T0P1T1BN.col(4)-=loopIter->second->B.template cast<float>();
//                    }
//                    else
//                    {
//                        assert(0);
//                    }
//                }
//                //                P0T0P1T1BN.col(5) = edgeRow.segment<dim>(1*dim).transpose().template cast<float>();		// plane normal
//                P0T0P1T1BN.col(5).setZero();		// plane normal
//                
//                VectorDim chord = P0T0P1T1BN.col(2)-P0T0P1T1BN.col(0);
//                VectorDim burgers=P0T0P1T1BN.col(4);
//                VectorDim planeNormal=P0T0P1T1BN.col(5);
//                
//                const float g = std::pow(chord.norm(),alpha);
//                
//                
//                //                    lines.push_back(vtkSmartPointer<vtkPolyLine>::New());
//                //                    auto& line(*lines.rbegin());
//                vtkSmartPointer<vtkPolyLine> line=vtkSmartPointer<vtkPolyLine>::New();
//                line->GetPointIds()->SetNumberOfIds(Np);
//                
//                //                    unsigned char clr0[3]={255,255,255};
//                
//                const float burgersNorm(burgers.norm());
//                
//                for (int k=0;k<Np;++k) // this may have to go to Np+1
//                {
//                    const float u1=k*1.0/(Np-1);
//                    const float u2=u1*u1;
//                    const float u3=u2*u1;
//                    
//                    // Compute positions along axis
//                    VectorDim P =   ( 2.0f*u3-3.0f*u2+1.0f) * P0T0P1T1BN.col(0)
//                    /*************/ + g*(      u3-2.0f*u2+u1)   * P0T0P1T1BN.col(1)
//                    /*************/ +   (-2.0f*u3+3.0f*u2)      * P0T0P1T1BN.col(2)
//                    /*************/ + g*(      u3-u2)           * P0T0P1T1BN.col(3);
//                    
//                    points->InsertNextPoint(P.data());
//                    radii->InsertNextValue(burgersNorm*tubeRadius);
//                    line->GetPointIds()->SetId(k,ptID);
//                    
//                    ptID++;
//                }
//                
//                
//                if(burgers.squaredNorm()>FLT_EPSILON)
//                {
//                    
//                    Eigen::Matrix<int,3,1> colorVector=computeColor(burgers,chord,planeNormal);
//                    
//                    //                    unsigned char lineClr[3]={51,153,255};
//                    unsigned char lineClr[3]={(unsigned char) colorVector(0),(unsigned char) colorVector(1),(unsigned char) colorVector(2)};
//                    
//                    if(meshLocation>=1)
//                    {
//                        cellsBnd->InsertNextCell(line);
//                        colorsBnd->InsertNextTypedTuple(lineClr);
//                        
//                    }
//                    else
//                    {
//                        cells->InsertNextCell(line);
//                        
//                        if(meshLocation==2 && blackGrainBoundarySegments)
//                        {
//                            unsigned char lineClr1[3]={1,1,1};
//                            colors->InsertNextTypedTuple(lineClr1);
//                        }
//                        else
//                        {
//                            colors->InsertNextTypedTuple(lineClr);
//                        }
//                    }
//                    
//                    //                        if(meshLocation==1)
//                    //                        {
//                    //                            cellsBnd->InsertNextCell(line);
//                    //                            colorsBnd->InsertNextTypedTuple(lineClr);
//                    //
//                    //                        }
//                    //                        else
//                    //                        {
//                    //                            cells->InsertNextCell(line);
//                    //
//                    //                            if(meshLocation==2 && blackGrainBoundarySegments)
//                    //                            {
//                    //                                unsigned char lineClr1[3]={1,1,1};
//                    //                                colors->InsertNextTypedTuple(lineClr1);
//                    //                            }
//                    //                            else
//                    //                            {
//                    //                                colors->InsertNextTypedTuple(lineClr);
//                    //                            }
//                    //                        }
//                }
//                else
//                {
//                    cells0->InsertNextCell(line);
//                }
//                
//                //                }
//            }
//            
//            polyData->SetPoints(points);
//            polyData->SetLines(cells);
//            //            polyData->GetCellData()->SetScalars(colors);
//            
//            //            polyData->GetPointData()->SetScalars(radii);
//            //            polyData->GetCellData()->SetScalars(radii);
//            polyData->GetCellData()->AddArray(colors);
//            polyData->GetPointData()->AddArray(radii);
//            polyData->GetPointData()->SetActiveScalars("TubeRadius");
//            
//            
//            polyDataBnd->SetPoints(points);
//            polyDataBnd->SetLines(cellsBnd);
//            polyDataBnd->GetCellData()->SetScalars(colorsBnd);
//            
//            
//            polyData0->SetPoints(points);
//            polyData0->SetLines(cells0);
//            //            polyData0->GetCellData()->SetScalars(colors);
//            
//        }
        
        /**********************************************************************/
        void createLoopLinks()
        {
            
            std::map<size_t,std::map<size_t,size_t>> loopEdgeMap;
            
            for(const auto& looplink : this->links())
            {
                loopEdgeMap[looplink.loopID].emplace(looplink.sourceID,looplink.sinkID);
                edgeMap()[std::make_pair(std::min(looplink.sourceID,looplink.sinkID),std::max(looplink.sourceID,looplink.sinkID))].insert(&looplink);
            }
            
            assert(loopEdgeMap.size()==this->links().size());
            
            
            size_t loopLumber=1;
            for(const auto& loop : this->loops())
            {
                
                //                loopMap().emplace(loop.sID,&loop);
                
                const auto loopFound=loopEdgeMap.find(loop.sID);
                assert(loopFound!=loopEdgeMap.end());
                
                std::vector<size_t> nodeIDs;
                nodeIDs.push_back(loopFound->second.begin()->first);
                std::deque<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d>> nodePositions;
                auto vIter0(this->nodeMap().find(loopFound->second.begin()->first)); //source
                assert(vIter0!=this->nodeMap().end() && "VERTEX NOT FOUND IN V-FILE");
                nodePositions.push_back(vIter0->second->P);
                for(size_t k=0;k<loopFound->second.size();++k)
                {
                    const auto nodeFound=loopFound->second.find(*nodeIDs.rbegin());
                    if(k<loopFound->second.size()-1)
                    {
                        nodeIDs.push_back(nodeFound->second);
                        
                        auto vIter(this->nodeMap().find(nodeFound->second)); //source
                        assert(vIter!=this->nodeMap().end() && "VERTEX NOT FOUND IN V-FILE");
                        //Eigen::Map<const Eigen::Matrix<double,1,10>> vertexRow(vIter->second.data());
                        nodePositions.push_back(vIter->second->P);
                    }
                    else
                    {
                        assert(nodeFound->second==nodeIDs[0]);
                    }
                }
                
                if(fabs(loop.B.dot(loop.N))<FLT_EPSILON)
                {
                    PlanarPolygon pp(loop.B,loop.N);
                    pp.assignPoints(nodePositions);
                    std::deque<std::array<size_t, 3>> tri=pp.triangulate();
                    
                    for(const auto& triID : tri)
                    {
                        
                        const size_t& nodeID0(nodeIDs[std::get<0>(triID)]);
                        const size_t& nodeID1(nodeIDs[std::get<1>(triID)]);
                        const size_t& nodeID2(nodeIDs[std::get<2>(triID)]);
                        
                        const size_t ptID0=std::distance(this->nodeMap().begin(),this->nodeMap().find(nodeID0));
                        const size_t ptID1=std::distance(this->nodeMap().begin(),this->nodeMap().find(nodeID1));
                        const size_t ptID2=std::distance(this->nodeMap().begin(),this->nodeMap().find(nodeID2));
                        
                        vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
                        triangle->GetPointIds()->SetId ( 0, ptID0 );
                        triangle->GetPointIds()->SetId ( 1, ptID1 );
                        triangle->GetPointIds()->SetId ( 2, ptID2 );
                        triangles->InsertNextCell ( triangle );
                    }
                }

            }
            
            trianglePolyData->SetPoints ( nodePoints );
            trianglePolyData->SetPolys ( triangles );
            trianglePolyData->Modified();
        }
        
        /**********************************************************************/
        DislocationSegmentActor(const size_t& frameID,
                                /*const size_t& lastFrameID,
                                 const float& anglePerStep,
                                 Eigen::Matrix<float,3,1> spinAxis,*/
                                vtkRenderer* const ren) :
        /* init */ renderer(ren),
        /* init */ lineActor(vtkSmartPointer<vtkActor>::New()),
        /* init */ tubeActor(vtkSmartPointer<vtkActor>::New()),
        /* init */ tubeActorBnd(vtkSmartPointer<vtkActor>::New()),
        /* init */ lineActor0(vtkSmartPointer<vtkActor>::New()),
        /* init */ tubeActor0(vtkSmartPointer<vtkActor>::New()),
        /* init */ points(vtkSmartPointer<vtkPoints>::New()),
        /* init */ cells(vtkSmartPointer<vtkCellArray>::New()),
        /* init */ cellsBnd(vtkSmartPointer<vtkCellArray>::New()),
        /* init */ polyData(vtkSmartPointer<vtkPolyData>::New()),
        /* init */ polyDataBnd(vtkSmartPointer<vtkPolyData>::New()),
        /* init */ cells0(vtkSmartPointer<vtkCellArray>::New()),
        /* init */ polyData0(vtkSmartPointer<vtkPolyData>::New()),
        /* init */ radii(vtkSmartPointer<vtkFloatArray>::New()),
        /* init */ colors(vtkSmartPointer<vtkUnsignedCharArray>::New()),
        /* init */ colorsBnd(vtkSmartPointer<vtkUnsignedCharArray>::New()),
        /* init */ lineMapper(vtkSmartPointer<vtkPolyDataMapper>::New()),
        /* init */ tubeFilter(vtkSmartPointer<vtkTubeFilter>::New()),
        /* init */ tubeFilterBnd(vtkSmartPointer<vtkTubeFilter>::New()),
        /* init */ tubeMapper(vtkSmartPointer<vtkPolyDataMapper>::New()),
        /* init */ tubeMapperBnd(vtkSmartPointer<vtkPolyDataMapper>::New()),
        /* init */ lineMapper0(vtkSmartPointer<vtkPolyDataMapper>::New()),
        /* init */ tubeFilter0(vtkSmartPointer<vtkTubeFilter>::New()),
        /* init */ tubeMapper0(vtkSmartPointer<vtkPolyDataMapper>::New()),
        /* init */ nodePoints(vtkSmartPointer<vtkPoints>::New()),
        //        /* init */ nodeLabels(vtkSmartPointer<vtkStringArray>::New()),
        /* init */ sphereSource(vtkSmartPointer<vtkSphereSource>::New()),
        /* init */ nodeData(vtkSmartPointer<vtkPolyData>::New()),
        /* init */ nodeColors(vtkSmartPointer<vtkUnsignedCharArray>::New()),
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
        /* init */ labelActor(vtkSmartPointer<vtkActor2D>::New()),
        /* init */ singleNodePoint(vtkSmartPointer<vtkPoints>::New()),
        /* init */ singleNodeLabelPolyData(vtkSmartPointer<vtkPolyData>::New()),
        /* init */ singleNodeLabelScalars(vtkSmartPointer<vtkDoubleArray>::New()),
        /* init */ singleNodeLabelMapper(vtkSmartPointer<vtkLabeledDataMapper>::New()),
        /* init */ singleNodeLabelActor(vtkSmartPointer<vtkActor2D>::New()),
        /* init */ triangles(vtkSmartPointer<vtkCellArray>::New()),
        /* init */ trianglePolyData(vtkSmartPointer<vtkPolyData>::New()),
        /* init */ triangleMapper(vtkSmartPointer<vtkPolyDataMapper>::New()),
        /* init */ triangleActor(vtkSmartPointer<vtkActor>::New())
        {
            
            radii->SetName("TubeRadius");
            
            colors->SetName("Colors");
            colors->SetNumberOfComponents(3);
            colorsBnd->SetNumberOfComponents(3);
            nodeColors->SetNumberOfComponents(3);
            velocityVectors->SetNumberOfComponents(3);
            velocityVectors->SetName("nodeVelocity");
            velocityColors->SetNumberOfComponents(3);
            
            labelScalars->SetNumberOfComponents(1);
            labelScalars->SetName("node IDs");
            
            
            //            nodeLabels->SetName("node IDs");
            
            //            this->readBin(frameID);
            if(EVLio<3>::isBinGood(frameID))
            {
                this->readBin(frameID);
            }
            else
            {
                //                assert(EVLio<3>::isTxtGood(frameID));
                this->readTxt(frameID);
            }
            
            createNodes();
            createLoopLinks();
            createSegments();
            //            readSegments(frameID);
            //            //            if(showSlippedArea)
            //            //            {
            //            readLoopLinks(frameID);
            //            }
            
            // Populate polyData
            
            // Segments
            tubeFilter->SetInputData(polyData);
            tubeFilter->SetRadius(tubeRadius); // this must be a function similar to setColor
            //            tubeFilter->SetRadiusFactor(tubeRadius);
            //            tubeFilter->SetVaryRadiusToVaryRadiusByScalar();
            if(scaleRadiusByBurgers)
            {
                tubeFilter->SetVaryRadiusToVaryRadiusByAbsoluteScalar();
            }
            //            tubeFilter->SetColorModeToColorByVector();
            tubeFilter->SetNumberOfSides(10);
            tubeFilter->Update();
            tubeMapper->SetInputConnection(tubeFilter->GetOutputPort());
            tubeMapper->ScalarVisibilityOn();
            
            tubeMapper->SetScalarModeToUseCellFieldData();
            tubeMapper->SelectColorArray("Colors");
            
            
            lineMapper->SetInputData(polyData);
            tubeActor->SetMapper(tubeMapper);
            //tube->GetProperty()->SetColor(colorVector(0),colorVector(1),colorVector(2)); // Give some color to the tube
            lineActor->SetMapper(lineMapper);
            
            // Boundary segments
            tubeFilterBnd->SetInputData(polyDataBnd);
            tubeFilterBnd->SetRadius(tubeRadius); // this must be a function similar to setColor
            tubeFilterBnd->SetNumberOfSides(10);
            tubeFilterBnd->Update();
            tubeMapperBnd->SetInputConnection(tubeFilterBnd->GetOutputPort());
            tubeMapperBnd->ScalarVisibilityOn();
            tubeActorBnd->SetMapper(tubeMapperBnd);
            //            tubeActorBnd->GetProperty()->SetColor(0.5, 0.0, 0.5); //(R,G,B)
            tubeActorBnd->GetProperty()->SetOpacity(0.3); //(R,G,B)
            
            // Zero-Burgers Segments
            tubeFilter0->SetInputData(polyData0);
            tubeFilter0->SetRadius(tubeRadius); // this must be a function similar to setColor
            tubeFilter0->SetNumberOfSides(10);
            tubeFilter0->Update();
            tubeMapper0->SetInputConnection(tubeFilter0->GetOutputPort());
            tubeMapper0->ScalarVisibilityOn();
            lineMapper0->SetInputData(polyData0);
            tubeActor0->SetMapper(tubeMapper0);
            tubeActor0->GetProperty()->SetColor(0.5, 0.5, 0.5); //(R,G,B)
            tubeActor0->GetProperty()->SetOpacity(0.3); //(R,G,B)
            lineActor0->SetMapper(lineMapper0);
            
            
            // Nodes
            nodeGlyphs->SetSourceConnection(sphereSource->GetOutputPort());
            nodeGlyphs->SetInputData(nodeData);
            nodeGlyphs->ScalingOn();
            nodeGlyphs->SetScaleModeToScaleByVector();
            nodeGlyphs->SetScaleFactor(2.0*tubeRadius*1.2);
            //            nodeGlyphs->OrientOn();
            //            nodeGlyphs->ClampingOff();
            //            nodeGlyphs->SetVectorModeToUseVector();
            //            nodeGlyphs->SetIndexModeToOff();
            nodeGlyphs->SetColorModeToColorByScalar();
            nodeGlyphs->Update();
            //            nodeGlyphs->SetColorModeToColorByVector();
            nodeMapper->SetInputConnection(nodeGlyphs->GetOutputPort());
            //            nodeMapper->ScalarVisibilityOff();
            nodeActor->SetMapper(nodeMapper);
            
            
            // Velocities
            velocityGlyphs->SetSourceConnection(velocityArrowSource->GetOutputPort());
            velocityGlyphs->SetInputData(velocityPolyData);
            velocityGlyphs->ScalingOn();
            velocityGlyphs->SetScaleModeToScaleByVector();
            //            velocityGlyphs->SetColorModeToColorByScalar();
            //            velocityGlyphs->SetColorModeToColorByScalar();
            //velocityGlyphs->SetColorModeToColorByVector();
            //            velocityGlyphs->SetScaleFactor(velocityFactor);
            velocityGlyphs->OrientOn();
            velocityGlyphs->ClampingOff();
            velocityGlyphs->SetVectorModeToUseVector();
            velocityGlyphs->SetIndexModeToOff();
            velocityMapper->SetInputConnection(velocityGlyphs->GetOutputPort());
            velocityMapper->ScalarVisibilityOff();
            velocityActor->SetMapper(velocityMapper);
            velocityActor->GetProperty()->SetColor(1.0, 0.0, 1.0); //(R,G,B)
            
            
            
            // Labels
            labelMapper->SetInputData(labelPolyData);
            labelMapper->SetLabelModeToLabelScalars();
            labelMapper->SetLabelFormat("%1.0f");
            labelActor->SetMapper(labelMapper);
            labelActor->GetProperty()->SetColor(0.0, 0.0, 0.0); //(R,G,B)
            
            
            
            // Single node Label
            singleNodeLabelMapper->SetInputData(singleNodeLabelPolyData);
            singleNodeLabelMapper->SetLabelModeToLabelScalars();
            singleNodeLabelMapper->SetLabelFormat("%1.0f");
            singleNodeLabelActor->SetMapper(singleNodeLabelMapper);
            singleNodeLabelActor->GetProperty()->SetColor(1.0, 0.0, 0.0); //(R,G,B)
            singleNodeLabelActor->VisibilityOff();
            
            
            
            
            // SlippedArea
            triangleMapper->SetInputData(trianglePolyData);
            triangleActor->SetMapper(triangleMapper);
            
            
            
            //            const double axisNorm(spinAxis.norm());
            //            if(anglePerStep && axisNorm>0.0)
            //            {
            //                spinAxis/=axisNorm;
            //                //                std::cout<<"rotating"<< frameID*anglePerStep<<std::endl;
            //                const double splinAngle((frameID-lastFrameID)*anglePerStep);
            //                std::cout<<"ddSegments "<<splinAngle<<std::endl;
            //
            //                lineActor->RotateWXYZ(splinAngle,spinAxis(0),spinAxis(1),spinAxis(2));
            //                tubeActor->RotateWXYZ(splinAngle,spinAxis(0),spinAxis(1),spinAxis(2));
            //                tubeActorBnd->RotateWXYZ(splinAngle,spinAxis(0),spinAxis(1),spinAxis(2));
            //                tubeActor0->RotateWXYZ(splinAngle,spinAxis(0),spinAxis(1),spinAxis(2));
            //                lineActor0->RotateWXYZ(splinAngle,spinAxis(0),spinAxis(1),spinAxis(2));
            //                nodeActor->RotateWXYZ(splinAngle,spinAxis(0),spinAxis(1),spinAxis(2));
            //                velocityActor->RotateWXYZ(splinAngle,spinAxis(0),spinAxis(1),spinAxis(2));
            //                //labelActor->RotateWXYZ(splinAngle,spinAxis(0),spinAxis(1),spinAxis(2));
            //                //singleNodeLabelActor->RotateWXYZ(splinAngle,spinAxis(0),spinAxis(1),spinAxis(2));
            //                triangleActor->RotateWXYZ(splinAngle,spinAxis(0),spinAxis(1),spinAxis(2));
            //            }
            renderer->AddActor(tubeActor);
            renderer->AddActor(tubeActorBnd);
            renderer->AddActor(tubeActor0);
            renderer->AddActor(nodeActor);
            renderer->AddActor(velocityActor);
            renderer->AddActor(labelActor);
            renderer->AddActor(singleNodeLabelActor);
            renderer->AddActor(triangleActor);
            
            
            modify();
        }
        
        /**********************************************************************/
        ~DislocationSegmentActor()
        {
            renderer->RemoveActor(tubeActor);
            renderer->RemoveActor(tubeActorBnd);
            renderer->RemoveActor(tubeActor0);
            renderer->RemoveActor(nodeActor);
            renderer->RemoveActor(velocityActor);
            renderer->RemoveActor(labelActor);
            renderer->RemoveActor(singleNodeLabelActor);
            renderer->RemoveActor(triangleActor);
            
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
            
            triangleActor->GetProperty()->SetColor(0.5,0.0,0.5);
            triangleActor->GetProperty()->SetOpacity(slippedAreaOpacity);
            
            if(showSlippedArea)
            {
                triangleActor->VisibilityOn();
            }
            else
            {
                triangleActor->VisibilityOff();
            }
            
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
            
            
            if(showSingleNode)
            {
                // HERE WE SHOULD CHANGE THE NODE POSITION BASED ON NODE ID
                // OTHERWISE THE SELECTED NODE WILL BE VISIBLE ONLY UPON LOADING A NEW FRAME
                std::cout<<"RELOAD FRAME TO SHOW SELECTED NODE"<<std::endl;
                singleNodeLabelActor->VisibilityOn();
                
            }
            else
            {
                singleNodeLabelActor->VisibilityOff();
            }
            
            if(showNodes)
            {
                // HERE WE SHOULD CHANGE THE NODE POSITION BASED ON NODE ID
                // OTHERWISE THE SELECTED NODE WILL BE VISIBLE ONLY UPON LOADING A NEW FRAME
                std::cout<<"RELOAD FRAME TO SHOW SELECTED NODE"<<std::endl;
                nodeActor->VisibilityOn();
                
            }
            else
            {
                nodeActor->VisibilityOff();
            }
            
        }
        
    };
    
    // Static data members
    float DislocationSegmentActor::tubeRadius=5.0;
    bool DislocationSegmentActor::scaleRadiusByBurgers=false;
    float DislocationSegmentActor::alpha=0.5;
    DislocationSegmentActor::ColorScheme DislocationSegmentActor::clr=DislocationSegmentActor::colorBurgers;
    size_t DislocationSegmentActor::Np=2;
    bool DislocationSegmentActor::showBoundarySegments=false;
    bool DislocationSegmentActor::blackGrainBoundarySegments=false;
    bool DislocationSegmentActor::showVelocities=false;
    bool DislocationSegmentActor::showNodeIDs=false;
    float DislocationSegmentActor::velocityFactor=100.0;
    bool DislocationSegmentActor::showZeroBuergers=false;
    bool DislocationSegmentActor::showSingleNode=false;
    size_t DislocationSegmentActor::singleNodeID=0;
    bool DislocationSegmentActor::showNodes=false;
    bool DislocationSegmentActor::showSlippedArea=false;
    float DislocationSegmentActor::slippedAreaOpacity=0.1;
    
} // namespace model
#endif







