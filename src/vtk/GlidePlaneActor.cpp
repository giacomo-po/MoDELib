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
#include <vtkStructuredGrid.h>
#include <vtkLookupTable.h>
#include <vtkStructuredGridAppend.h>

#include <GlidePlaneActor.h>

namespace model
{

    GlidePlaneActor::GlidePlaneActor(vtkGenericOpenGLRenderWindow* const renWin,vtkRenderer* const ren,const Polycrystal<3>& poly,const DDtraitsIO& traitsIO) :
    /* init */ renderWindow(renWin)
    /* init */,renderer(ren)
    /* init */,mainLayout(new QGridLayout(this))
    /* init */,showGlidePlanes(new QCheckBox(this))
    /* init */,showGlidePlanesNoise(new QCheckBox(this))
    /* init */,glidePlanesNoiseBox(new QComboBox(this))
    /* init */,glidePlanePolydata(vtkSmartPointer<vtkPolyData>::New())
    /* init */,glidePlaneMapper(vtkSmartPointer<vtkPolyDataMapper>::New())
    /* init */,glidePlaneActor(vtkSmartPointer<vtkActor>::New())
    /* init */,glidePlaneFactory(poly)
    /* init */,planeNoise(traitsIO,poly)
    {
        
        showGlidePlanes->setChecked(false);
        showGlidePlanes->setText("show GlidePlanes");
        glidePlaneActor->SetVisibility(false);
        
        showGlidePlanesNoise->setChecked(false);
        showGlidePlanesNoise->setText("show GlidePlanesNoise");
        glidePlanesNoiseBox->setEnabled(false);
        if(planeNoise.solidSolution)
        {
            glidePlanesNoiseBox->addItem("solidSolution stress_xz");
            glidePlanesNoiseBox->addItem("solidSolution stress_yz");
        }
        if(planeNoise.stackingFault)
        {
            glidePlanesNoiseBox->addItem("ISF energy");
        }
        
        mainLayout->addWidget(showGlidePlanes,0,0,1,1);
        mainLayout->addWidget(showGlidePlanesNoise,1,0,1,1);
        mainLayout->addWidget(glidePlanesNoiseBox,1,1,1,1);
        
        this->setLayout(mainLayout);
        connect(showGlidePlanes,SIGNAL(stateChanged(int)), this, SLOT(modify()));
        connect(showGlidePlanesNoise,SIGNAL(stateChanged(int)), this, SLOT(modify()));
        connect(glidePlanesNoiseBox,SIGNAL(currentIndexChanged(int)), this, SLOT(modify()));
        
        
        
        glidePlaneMapper->SetInputData(glidePlanePolydata);
        glidePlaneActor->SetMapper(glidePlaneMapper);
        glidePlaneActor->GetProperty()->SetColor(1.0, 0.0, 1.0); //(R,G,B)
        
        renderer->AddActor(glidePlaneActor);
        //        renderer->AddActor(solidSolutionNoiseActorXZ);
        
        
    }

    void GlidePlaneActor::updateConfiguration(const DDauxIO<3>& auxIO)
    {
        
        for(const auto& actor : noiseActors)
        {
            renderer->RemoveActor(actor);
        }
        noiseMappers.clear();
        noiseActors.clear();
        
        vtkSmartPointer<vtkPoints> glidePlanePoints = vtkSmartPointer<vtkPoints>::New();
        vtkSmartPointer<vtkCellArray> glidePlaneCells(vtkSmartPointer<vtkCellArray>::New());
        //        vtkNew<vtkStructuredGridAppend> gridsAppend;
        
        
        size_t pointIDs(0);
        double noiseValMin(std::numeric_limits<double>::max());
        double noiseValMax(std::numeric_limits<double>::min());
        for(const auto& gpIO : auxIO.glidePlanes())
        {
            // Plot GlidePlaneBoundary
            GlidePlaneKey<3> planeKey(gpIO.r,gpIO.h,gpIO.latticeID);
            const auto glidePlane(glidePlaneFactory.getFromKey(planeKey));
            if(glidePlane->slipSystems().size())
            {
                
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
            

                // Plot solidSolution-noise
                vtkNew<vtkPoints> glidePlaneNoisePoints;
                
                vtkNew<vtkDoubleArray> glidePlaneNoisePointsValuesXZ;
                vtkNew<vtkDoubleArray> glidePlaneNoisePointsValuesYZ;
                if(planeNoise.solidSolution)
                {
                    glidePlaneNoisePointsValuesXZ->SetNumberOfValues(planeNoise.gridSize.array().prod());
                    glidePlaneNoisePointsValuesXZ->SetName(glidePlanesNoiseBox->itemText(0).toStdString().c_str());
                    
                    glidePlaneNoisePointsValuesYZ->SetNumberOfValues(planeNoise.gridSize.array().prod());
                    glidePlaneNoisePointsValuesYZ->SetName(glidePlanesNoiseBox->itemText(1).toStdString().c_str());
                }
                
                vtkNew<vtkDoubleArray> sfNoise;
                if(planeNoise.stackingFault)
                {
                    sfNoise->SetNumberOfValues(planeNoise.gridSize.array().prod());
                    sfNoise->SetName(glidePlanesNoiseBox->itemText(2).toStdString().c_str());
                }
                


                
                for (int k = 0; k < planeNoise.gridSize(2); k++)
                {
                    for (int j = 0; j < planeNoise.gridSize(1); j++)
                    {
                        for (int i = 0; i < planeNoise.gridSize(0); i++)
                        {
                            
                            const double x(i*planeNoise.gridSpacing(0));
                            const double y(j*planeNoise.gridSpacing(1));
                            
                            const Eigen::Matrix<double,2,1> localPos((Eigen::Matrix<double,2,1>()<<x,y).finished());
                            const Eigen::Matrix<double,3,1> globalPos(glidePlane->globalPosition(localPos));
                            
                            glidePlaneNoisePoints->InsertNextPoint(globalPos(0), globalPos(1), globalPos(2));
                            
                            const auto ind = planeNoise.gridSize(1)*planeNoise.gridSize(2)*i + j*planeNoise.gridSize(2) + k;
                            
                            if(planeNoise.solidSolution)
                            {
                                const auto& noiseVal(planeNoise.solidSolution->operator[](ind));
                                glidePlaneNoisePointsValuesXZ->SetValue(ind, noiseVal(0));
                                glidePlaneNoisePointsValuesYZ->SetValue(ind, noiseVal(1));
                                noiseValMin=std::min(noiseValMin,std::min(noiseVal(0),noiseVal(1)));
                                noiseValMax=std::max(noiseValMax,std::max(noiseVal(0),noiseVal(1)));
                            }
                            
                            if(planeNoise.stackingFault)
                            {
                                const auto& noiseVal(planeNoise.stackingFault->operator[](ind));
                                sfNoise->SetValue(ind, noiseVal);
                                                            noiseValMin=std::min(noiseValMin,noiseVal);
                                                            noiseValMax=std::max(noiseValMax,noiseVal);
                            }

                            
                        }
                    }
                }
                vtkNew<vtkStructuredGrid> glidePlaneNoiseGrid;
                glidePlaneNoiseGrid->SetDimensions(planeNoise.gridSize(0),planeNoise.gridSize(1),planeNoise.gridSize(2));
                glidePlaneNoiseGrid->SetPoints(glidePlaneNoisePoints);
                glidePlaneNoiseGrid->GetPointData()->AddArray(glidePlaneNoisePointsValuesXZ);
                glidePlaneNoiseGrid->GetPointData()->AddArray(glidePlaneNoisePointsValuesYZ);
                glidePlaneNoiseGrid->GetPointData()->AddArray(sfNoise);

                
                
                vtkNew<vtkLookupTable> glidePlaneNoiseLut;
                glidePlaneNoiseLut->SetNumberOfTableValues(planeNoise.gridSize.array().prod());
                glidePlaneNoiseLut->Build();
                
                noiseMappers.push_back(vtkSmartPointer<vtkDataSetMapper>::New());
                noiseMappers.back()->SetInputData(glidePlaneNoiseGrid);
                noiseMappers.back()->SetScalarModeToUsePointFieldData();
                noiseMappers.back()->SelectColorArray(1);
                noiseMappers.back()->SetLookupTable(glidePlaneNoiseLut);
                noiseMappers.back()->SetScalarRange(noiseValMin, noiseValMax);
                noiseMappers.back()->ScalarVisibilityOn();
                noiseActors.push_back(vtkSmartPointer<vtkActor>::New());
                noiseActors.back()->SetMapper(noiseMappers.back());
                noiseActors.back()->SetVisibility(showGlidePlanesNoise->isChecked());
                renderer->AddActor( noiseActors.back());
                
            }
        }
        
        glidePlanePolydata->SetPoints(glidePlanePoints);
        glidePlanePolydata->SetLines(glidePlaneCells);
        glidePlanePolydata->Modified();
        
        
    }

    void GlidePlaneActor::modify()
    {
        glidePlaneActor->SetVisibility(showGlidePlanes->isChecked());
        
        glidePlanesNoiseBox->setEnabled(showGlidePlanesNoise->isChecked());
        for(const auto& mapper : noiseMappers)
        {
            mapper->SelectColorArray(glidePlanesNoiseBox->currentIndex());
        }
        
        
        for(const auto& actor : noiseActors)
        {
            actor->SetVisibility(showGlidePlanesNoise->isChecked());
        }
        
        
        renderWindow->Render();
    }


} // namespace model
#endif
