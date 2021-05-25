/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DDauxVtk_H_
#define model_DDauxVtk_H_

#include <Eigen/Dense>

#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkPolyData.h>
#include <vtkArrowSource.h>
#include <vtkGlyph3D.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkTriangle.h>
#include <DDauxIO.h>


namespace model
{
    
    /**************************************************************************/
    /**************************************************************************/
    struct DDauxVtk : public DDauxIO<3>
    {
        
        vtkRenderer* const renderer;
        
        vtkSmartPointer<vtkPoints> glidePlanePoints;
        vtkSmartPointer<vtkPolyData> glidePlanePolydata;
        vtkSmartPointer<vtkPolyDataMapper> glidePlaneMapper;
        vtkSmartPointer<vtkActor> glidePlaneActor;

        vtkSmartPointer<vtkPoints> periodicGlidePlanePoints;
        vtkSmartPointer<vtkPolyData> periodicGlidePlanePolydata;
        vtkSmartPointer<vtkPolyDataMapper> periodicGlidePlaneMapper;
        vtkSmartPointer<vtkActor> periodicGlidePlaneActor;

        vtkSmartPointer<vtkPoints> periodicLoopPoints;
        vtkSmartPointer<vtkPolyData> periodicLoopPolydata;
        vtkSmartPointer<vtkPolyDataMapper> periodicLoopMapper;
        vtkSmartPointer<vtkActor> periodicLoopActor;

        
        vtkSmartPointer<vtkPoints> quadraturePositions;
        vtkSmartPointer<vtkArrowSource> arrowSource;

        vtkSmartPointer<vtkDoubleArray> quadraturePk;
        vtkSmartPointer<vtkPolyData> quadraturePkPolyData;
        vtkSmartPointer<vtkGlyph3D> quadraturePkGlyph;
        vtkSmartPointer<vtkPolyDataMapper> quadraturePkMapper;
        vtkSmartPointer<vtkActor> quadraturePkActor;
        
        vtkSmartPointer<vtkDoubleArray> quadratureGlideVelocities;
        vtkSmartPointer<vtkPolyData> quadratureGlideVelocitiesPolyData;
        vtkSmartPointer<vtkGlyph3D> quadratureGlideVelocitiesGlyph;
        vtkSmartPointer<vtkPolyDataMapper> quadratureGlideVelocitiesMapper;
        vtkSmartPointer<vtkActor> quadratureGlideVelocitiesActor;

        std::map<size_t,std::vector<const GlidePlaneBoundaryIO<3>*>> glidePlaneBoundaryMap;
        
        static bool showGlidePlanes;
        static float glidePlaneOpacity;
        static bool showPeriodicGlidePlanes;
        static bool showPeriodicLoops;
        static bool showPkforces;
        static bool showGlideVelocities;
        static float pkFactor;
        static float glideVelocitiesFactor;


        /**********************************************************************/
        DDauxVtk(const size_t& frameID,
                 vtkRenderer* const ren) :
        /* init */ renderer(ren)
        /* init */,glidePlanePoints(vtkSmartPointer<vtkPoints>::New())
        /* init */,glidePlanePolydata(vtkSmartPointer<vtkPolyData>::New())
        /* init */,glidePlaneMapper(vtkSmartPointer<vtkPolyDataMapper>::New())
        /* init */,glidePlaneActor(vtkSmartPointer<vtkActor>::New())
        /* init */,periodicGlidePlanePoints(vtkSmartPointer<vtkPoints>::New())
        /* init */,periodicGlidePlanePolydata(vtkSmartPointer<vtkPolyData>::New())
        /* init */,periodicGlidePlaneMapper(vtkSmartPointer<vtkPolyDataMapper>::New())
        /* init */,periodicGlidePlaneActor(vtkSmartPointer<vtkActor>::New())
        /* init */,periodicLoopPoints(vtkSmartPointer<vtkPoints>::New())
        /* init */,periodicLoopPolydata(vtkSmartPointer<vtkPolyData>::New())
        /* init */,periodicLoopMapper(vtkSmartPointer<vtkPolyDataMapper>::New())
        /* init */,periodicLoopActor(vtkSmartPointer<vtkActor>::New())
        /* init */,quadraturePositions(vtkSmartPointer<vtkPoints>::New())
        /* init */,arrowSource(vtkSmartPointer<vtkArrowSource>::New())
        /* init */,quadraturePk(vtkSmartPointer<vtkDoubleArray>::New())
        /* init */,quadraturePkPolyData(vtkSmartPointer<vtkPolyData>::New())
        /* init */,quadraturePkGlyph(vtkSmartPointer<vtkGlyph3D>::New())
        /* init */,quadraturePkMapper(vtkSmartPointer<vtkPolyDataMapper>::New())
        /* init */,quadraturePkActor(vtkSmartPointer<vtkActor>::New())
        /* init */,quadratureGlideVelocities(vtkSmartPointer<vtkDoubleArray>::New())
        /* init */,quadratureGlideVelocitiesPolyData(vtkSmartPointer<vtkPolyData>::New())
        /* init */,quadratureGlideVelocitiesGlyph(vtkSmartPointer<vtkGlyph3D>::New())
        /* init */,quadratureGlideVelocitiesMapper(vtkSmartPointer<vtkPolyDataMapper>::New())
        /* init */,quadratureGlideVelocitiesActor(vtkSmartPointer<vtkActor>::New())
        {
            
            quadraturePk->SetNumberOfComponents(3);
            quadraturePk->SetName("PkForce");
            
            quadratureGlideVelocities->SetNumberOfComponents(3);
            quadratureGlideVelocities->SetName("glideVelocity");

            if(this->isBinGood(frameID))
            {
                this->readBin(frameID);
            }
            else
            {
                this->readTxt(frameID);
            }
            
            plotGlidePlaneBoundaries();
            plotPeriodicGlidePlaneBoundaries();
            plotPeriodicLoops();
            plotQuadraturePointData();

        }
        
        /**********************************************************************/
        ~DDauxVtk()
        {
            renderer->RemoveActor(glidePlaneActor);
            renderer->RemoveActor(periodicGlidePlaneActor);
            renderer->RemoveActor(periodicLoopActor);
            renderer->RemoveActor(quadraturePkActor);
            renderer->RemoveActor(quadratureGlideVelocitiesActor);
        }
        
        /**********************************************************************/
        void plotGlidePlaneBoundaries()
        {
            glidePlanePolydata->Allocate();
            size_t connectivityID=0;
            for(const auto& gpb : this->glidePlanesBoundaries())
            {
                glidePlaneBoundaryMap[gpb.glidePlaneID].push_back(&gpb);
                
                glidePlanePoints->InsertNextPoint(gpb.P0(0),
                                                  gpb.P0(1),
                                                  gpb.P0(2));
                
                glidePlanePoints->InsertNextPoint(gpb.P1(0),
                                                  gpb.P1(1),
                                                  gpb.P1(2));
                
                vtkIdType connectivity[2]; //points IDs of each line segment
                connectivity[0] = connectivityID;
                connectivity[1] = connectivityID+1;
                glidePlanePolydata->InsertNextCell(VTK_LINE,2,connectivity);
                connectivityID+=2;
            }
            glidePlanePolydata->SetPoints(glidePlanePoints);
            glidePlaneMapper->SetInputData(glidePlanePolydata);
            glidePlaneActor->SetMapper(glidePlaneMapper);
            glidePlaneActor->GetProperty()->SetLineWidth(1.5);
            glidePlaneActor->GetProperty()->SetColor(0.5,0,0.7);
            glidePlaneActor->GetProperty()->SetOpacity(glidePlaneOpacity);
            glidePlaneActor->SetVisibility(showGlidePlanes);
            renderer->AddActor(glidePlaneActor);
        }
        
        /**********************************************************************/
        void plotPeriodicLoops()
        {
            
            std::map<size_t,size_t> periodicLoopMap; // nodeID,nodePositionInDDauxIO
            for(size_t k=0;k<periodicLoopNodes().size();++k)
            {
                const auto& node(periodicLoopNodes()[k]);
                periodicLoopMap.emplace(node.sID,k);

                periodicLoopPoints->InsertNextPoint(node.Pg(0),
                                                    node.Pg(1),
                                                    node.Pg(2));
            }
            
            periodicLoopPolydata->Allocate();
            size_t connectivityID=0;
            for(const auto& link : this->periodicLoopLinks())
            {
//                glidePlaneBoundaryMap[gpb.glidePlaneID].push_back(&gpb);
                
//                glidePlanePoints->InsertNextPoint(gpb.P0(0),
//                                                  gpb.P0(1),
//                                                  gpb.P0(2));
//
//                glidePlanePoints->InsertNextPoint(gpb.P1(0),
//                                                  gpb.P1(1),
//                                                  gpb.P1(2));
                
                const auto sourceIter(periodicLoopMap.find(link.sourceID));
                assert(sourceIter!=periodicLoopMap.end());
                const auto sinkIter(periodicLoopMap.find(link.sinkID));
                assert(sinkIter!=periodicLoopMap.end());

                vtkIdType connectivity[2]; //points IDs of each line segment
                connectivity[0] = sourceIter->second;
                connectivity[1] = sinkIter->second;
                periodicLoopPolydata->InsertNextCell(VTK_LINE,2,connectivity);
                connectivityID+=2;
            }
            periodicLoopPolydata->SetPoints(periodicLoopPoints);
            periodicLoopMapper->SetInputData(periodicLoopPolydata);
            periodicLoopActor->SetMapper(periodicLoopMapper);
            periodicLoopActor->GetProperty()->SetLineWidth(2.0);
            periodicLoopActor->GetProperty()->SetColor(0.0,0,1.0);
            periodicLoopActor->GetProperty()->SetOpacity(glidePlaneOpacity);
            periodicLoopActor->SetVisibility(showPeriodicLoops);
            renderer->AddActor(periodicLoopActor);
        }

        /**********************************************************************/
        void plotPeriodicGlidePlaneBoundaries()
        {
            periodicGlidePlanePolydata->Allocate();
            size_t connectivityID=0;
            for(const auto& patch : this->periodicGlidePlanePatches())
            {
                for(const auto& gpb : glidePlaneBoundaryMap[patch.glidePlaneID])
                {
                    
                    periodicGlidePlanePoints->InsertNextPoint(gpb->P0(0)-patch.shift(0),
                                                              gpb->P0(1)-patch.shift(1),
                                                              gpb->P0(2)-patch.shift(2));
                    
                    periodicGlidePlanePoints->InsertNextPoint(gpb->P1(0)-patch.shift(0),
                                                              gpb->P1(1)-patch.shift(1),
                                                              gpb->P1(2)-patch.shift(2));
                    
                    vtkIdType connectivity[2]; //points IDs of each line segment
                    connectivity[0] = connectivityID;
                    connectivity[1] = connectivityID+1;
                    periodicGlidePlanePolydata->InsertNextCell(VTK_LINE,2,connectivity);
                    connectivityID+=2;
                }
            }
            periodicGlidePlanePolydata->SetPoints(periodicGlidePlanePoints);
            periodicGlidePlaneMapper->SetInputData(periodicGlidePlanePolydata);
            periodicGlidePlaneActor->SetMapper(periodicGlidePlaneMapper);
            periodicGlidePlaneActor->GetProperty()->SetLineWidth(1.5);
            periodicGlidePlaneActor->GetProperty()->SetColor(0.0,0.5,0.7);
            periodicGlidePlaneActor->GetProperty()->SetOpacity(glidePlaneOpacity);
            periodicGlidePlaneActor->SetVisibility(showPeriodicGlidePlanes);
            renderer->AddActor(periodicGlidePlaneActor);
        }
        
        /**********************************************************************/
        void plotQuadraturePointData()
        {
            for(const auto& q : this->quadraturePoints())
            {
                quadraturePositions->InsertNextPoint(q.r.data());  // origin of arrow
                quadraturePk->InsertNextTuple(q.pkForce.data()); // arrow vactor
                quadratureGlideVelocities->InsertNextTuple(q.glideVelocity.data()); // arrow vactor
            }
            
            quadraturePkPolyData->SetPoints(quadraturePositions);
            quadraturePkPolyData->GetPointData()->SetVectors(quadraturePk);
            quadraturePkPolyData->Modified();
            quadraturePkGlyph->SetSourceConnection(arrowSource->GetOutputPort());
            quadraturePkGlyph->SetInputData(quadraturePkPolyData);
            quadraturePkGlyph->ScalingOn();
            quadraturePkGlyph->SetScaleModeToScaleByVector();
            quadraturePkGlyph->OrientOn();
            quadraturePkGlyph->ClampingOff();
            quadraturePkGlyph->SetVectorModeToUseVector();
            quadraturePkGlyph->SetIndexModeToOff();
            quadraturePkGlyph->SetScaleFactor(pkFactor);
            quadraturePkMapper->SetInputConnection(quadraturePkGlyph->GetOutputPort());
            quadraturePkMapper->ScalarVisibilityOff();
            quadraturePkActor->SetMapper(quadraturePkMapper);
            quadraturePkActor->GetProperty()->SetColor(0.0,0.0,1.0);
            quadraturePkActor->SetVisibility(showPkforces);
            renderer->AddActor(quadraturePkActor);
            
            quadratureGlideVelocitiesPolyData->SetPoints(quadraturePositions);
            quadratureGlideVelocitiesPolyData->GetPointData()->SetVectors(quadratureGlideVelocities);
            quadratureGlideVelocitiesPolyData->Modified();
            quadratureGlideVelocitiesGlyph->SetSourceConnection(arrowSource->GetOutputPort());
            quadratureGlideVelocitiesGlyph->SetInputData(quadratureGlideVelocitiesPolyData);
            quadratureGlideVelocitiesGlyph->ScalingOn();
            quadratureGlideVelocitiesGlyph->SetScaleModeToScaleByVector();
            quadratureGlideVelocitiesGlyph->OrientOn();
            quadratureGlideVelocitiesGlyph->ClampingOff();
            quadratureGlideVelocitiesGlyph->SetVectorModeToUseVector();
            quadratureGlideVelocitiesGlyph->SetIndexModeToOff();
            quadratureGlideVelocitiesGlyph->SetScaleFactor(glideVelocitiesFactor);
            quadratureGlideVelocitiesMapper->SetInputConnection(quadratureGlideVelocitiesGlyph->GetOutputPort());
            quadratureGlideVelocitiesMapper->ScalarVisibilityOff();
            quadratureGlideVelocitiesActor->SetMapper(quadratureGlideVelocitiesMapper);
            quadratureGlideVelocitiesActor->GetProperty()->SetColor(0.0,1.0,0.0);
            quadratureGlideVelocitiesActor->SetVisibility(showGlideVelocities);
            renderer->AddActor(quadratureGlideVelocitiesActor);

        }
        
        /**********************************************************************/
        void modify()
        {
            glidePlaneActor->SetVisibility(showGlidePlanes);
            glidePlaneActor->GetProperty()->SetOpacity(glidePlaneOpacity);

            quadraturePkGlyph->SetScaleFactor(pkFactor);
            quadraturePkActor->SetVisibility(showPkforces);
            
            quadratureGlideVelocitiesGlyph->SetScaleFactor(glideVelocitiesFactor);
            quadratureGlideVelocitiesActor->SetVisibility(showGlideVelocities);
            
            periodicGlidePlaneActor->SetVisibility(showPeriodicGlidePlanes);
            periodicGlidePlaneActor->GetProperty()->SetOpacity(glidePlaneOpacity);

            periodicLoopActor->SetVisibility(showPeriodicLoops);
            periodicLoopActor->GetProperty()->SetOpacity(glidePlaneOpacity);

            
        }

    };
    
    bool DDauxVtk::showGlidePlanes=false;
    float DDauxVtk::glidePlaneOpacity=0.5;
    bool DDauxVtk::showPeriodicGlidePlanes=false;
    bool DDauxVtk::showPeriodicLoops=false;
    bool DDauxVtk::showPkforces=false;
    bool DDauxVtk::showGlideVelocities=false;
    float DDauxVtk::pkFactor=0.5;
    float DDauxVtk::glideVelocitiesFactor=0.5;

}
#endif







