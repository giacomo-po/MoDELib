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

        std::map<size_t,std::vector<const GlidePlaneBoundaryIO<3>*>> glidePlaneBoundaryMap;
        static bool showGlidePlanes;
        static bool showPeriodicGlidePlanes;
        static float glidePlaneOpacity;

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
        {
            if(this->isBinGood(frameID))
            {
                this->readBin(frameID);
            }
            else
            {
                this->readTxt(frameID);
            }
            
            makeGlidePlaneBoundaries();
            makePeriodicGlidePlaneBoundaries();

        }
        
        /**********************************************************************/
        ~DDauxVtk()
        {
            renderer->RemoveActor(glidePlaneActor);
            renderer->RemoveActor(periodicGlidePlaneActor);
        }
        
        /**********************************************************************/
        void makeGlidePlaneBoundaries()
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
        void makePeriodicGlidePlaneBoundaries()
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
        void modify()
        {
            glidePlaneActor->SetVisibility(showGlidePlanes);
            glidePlaneActor->GetProperty()->SetOpacity(glidePlaneOpacity);

            periodicGlidePlaneActor->SetVisibility(showPeriodicGlidePlanes);
            periodicGlidePlaneActor->GetProperty()->SetOpacity(glidePlaneOpacity);

        }

    };
    
    bool DDauxVtk::showGlidePlanes=false;
    bool DDauxVtk::showPeriodicGlidePlanes=false;
    float DDauxVtk::glidePlaneOpacity=0.5;

}
#endif







