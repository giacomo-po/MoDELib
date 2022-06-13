/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_GlidePlaneActor_cpp_
#define model_GlidePlaneActor_cpp_

#include <vtkLine.h>
#include <GlidePlaneActor.h>

namespace model
{

    GlidePlaneActor::GlidePlaneActor(vtkGenericOpenGLRenderWindow* const renWin,vtkRenderer* const renderer,const Polycrystal<3>& poly,const DDtraitsIO& traitsIO) :
    /* init */ renderWindow(renWin)
    /* init */,mainLayout(new QGridLayout(this))
    /* init */,showGlidePlanes(new QCheckBox(this))
    /* init */,showCompositionNoise(new QCheckBox(this))
    /* init */,glidePlanePolydata(vtkSmartPointer<vtkPolyData>::New())
    /* init */,glidePlaneMapper(vtkSmartPointer<vtkPolyDataMapper>::New())
    /* init */,glidePlaneActor(vtkSmartPointer<vtkActor>::New())
    /* init */,glidePlaneFactory(poly)
    /* init */,planeNoise(traitsIO)
    {
        
        showGlidePlanes->setChecked(false);
        showGlidePlanes->setText("show GlidePlanes");
        glidePlaneActor->SetVisibility(false);

        showCompositionNoise->setChecked(false);
        showCompositionNoise->setText("show CompositionNoise");

        
        mainLayout->addWidget(showGlidePlanes,0,0,1,1);
        mainLayout->addWidget(showCompositionNoise,1,0,1,1);

        
        
        this->setLayout(mainLayout);
        connect(showGlidePlanes,SIGNAL(stateChanged(int)), this, SLOT(modify()));
        connect(showCompositionNoise,SIGNAL(stateChanged(int)), this, SLOT(modify()));

        glidePlaneMapper->SetInputData(glidePlanePolydata);
        glidePlaneActor->SetMapper(glidePlaneMapper);
        glidePlaneActor->GetProperty()->SetColor(1.0, 0.0, 1.0); //(R,G,B)

        renderer->AddActor(glidePlaneActor);

        
    }

    void GlidePlaneActor::updateConfiguration(const DDauxIO<3>& auxIO)
    {
        vtkSmartPointer<vtkPoints> glidePlanePoints = vtkSmartPointer<vtkPoints>::New();
        vtkSmartPointer<vtkCellArray> glidePlaneCells(vtkSmartPointer<vtkCellArray>::New());

        size_t pointIDs(0);
        for(const auto& gpIO : auxIO.glidePlanes())
        {
            GlidePlaneKey<3> planeKey(gpIO.r,gpIO.h,gpIO.latticeID);
            const auto glidePlane(glidePlaneFactory.getFromKey(planeKey));
            for(const auto& mshInt : glidePlane->meshIntersections)
            {
                glidePlanePoints->InsertNextPoint(mshInt->P0(0),
                                                  mshInt->P0(1),
                                                  mshInt->P0(2));
                
                glidePlanePoints->InsertNextPoint(mshInt->P1(0),
                                                  mshInt->P1(1),
                                                  mshInt->P1(2));
                
                vtkSmartPointer<vtkLine> line(vtkSmartPointer<vtkLine>::New());
                line->GetPointIds()->SetId(0, pointIDs); // the second 0 is the index of the Origin in linesPolyData's points
                line->GetPointIds()->SetId(1, pointIDs+1);
                glidePlaneCells->InsertNextCell(line);
                pointIDs+=2;
            }
        }
        
        
        glidePlanePolydata->SetPoints(glidePlanePoints);
        glidePlanePolydata->SetLines(glidePlaneCells);
        glidePlanePolydata->Modified();

        
    }

    void GlidePlaneActor::modify()
    {
        glidePlaneActor->SetVisibility(showGlidePlanes->isChecked());

        renderWindow->Render();
    }


} // namespace model
#endif
