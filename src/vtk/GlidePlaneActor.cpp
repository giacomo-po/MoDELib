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

    GlidePlaneActor::GlidePlaneActor(vtkGenericOpenGLRenderWindow* const renWin,vtkRenderer* const ren,const Polycrystal<3>& poly_in,const DDtraitsIO& traitsIO) :
    /* init */ renderWindow(renWin)
    /* init */,renderer(ren)
    /* init */,poly(poly_in)
    /* init */,mainLayout(new QGridLayout(this))
    /* init */,showGlidePlanes(new QCheckBox(this))
    /* init */,showGlidePlanesNoise(new QCheckBox(this))
    /* init */,glidePlanesNoiseBox(new QComboBox(this))
    /* init */,slipSystemNoiseBox(new QComboBox(this))
    /* init */,glidePlanePolydata(vtkSmartPointer<vtkPolyData>::New())
    /* init */,glidePlaneMapper(vtkSmartPointer<vtkPolyDataMapper>::New())
    /* init */,glidePlaneActor(vtkSmartPointer<vtkActor>::New())
    /* init */,glidePlaneFactory(poly)
//    /* init */,planeNoise(traitsIO.polyFile,poly)
    {
        
        showGlidePlanes->setChecked(false);
        showGlidePlanes->setText("show GlidePlanes");
        glidePlaneActor->SetVisibility(false);
        
        showGlidePlanesNoise->setChecked(false);
        showGlidePlanesNoise->setText("show GlidePlanesNoise");
        glidePlanesNoiseBox->setEnabled(false);
//        if(planeNoise->solidSolution)
//        {
//            glidePlanesNoiseBox->addItem("solidSolution stress_xz");
//            glidePlanesNoiseBox->addItem("solidSolution stress_yz");
            glidePlanesNoiseBox->addItem("RSS");

//        }
//        if(planeNoise->stackingFault)
//        {
            glidePlanesNoiseBox->addItem("ISF energy");
//        }
//        if(planeNoise->solidSolution || planeNoise->stackingFault)
//        {
            const auto& grain(poly.grains.begin()->second);
            
            int k=0;
            for(const auto& ss : grain.singleCrystal->slipSystems())
            {
                slipSystemNoiseBox->addItem(QString::fromStdString("slipSystem "+ std::to_string(k)));
                k++;
            }
//        }
        
        mainLayout->addWidget(showGlidePlanes,0,0,1,1);
        mainLayout->addWidget(showGlidePlanesNoise,1,0,1,1);
        mainLayout->addWidget(glidePlanesNoiseBox,1,1,1,1);
        mainLayout->addWidget(slipSystemNoiseBox,2,1,1,1);
        
        
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
        
        std::cout<<"Updating GlidePlanes..."<<std::flush;
        const auto t0= std::chrono::system_clock::now();
        
        
        for(const auto& actorPair : noiseActors)
        {
            for(const auto& actor : actorPair.second)
            {
                renderer->RemoveActor(actor);
            }
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
                for(size_t slipSystemID=0; slipSystemID<poly.grain(gpIO.latticeID).singleCrystal->slipSystems().size();++slipSystemID)
                {
                    const auto& slipSystem(poly.grain(gpIO.latticeID).singleCrystal->slipSystems()[slipSystemID]);
                    if(slipSystem->planeNoise)
                    {
                        if(slipSystem->planeNoise->solidSolution || slipSystem->planeNoise->stackingFault)
                        {
                            if(glidePlane->n.dot(slipSystem->s.dir)==0)
                            {// glide plane includes slip system
//                                vtkNew<vtkPoints> glidePlaneNoisePoints;
//
//                                vtkNew<vtkDoubleArray> glidePlaneNoisePointsValuesXZ;
//                                vtkNew<vtkDoubleArray> glidePlaneNoisePointsValuesYZ;
//                                if(slipSystem->planeNoise->solidSolution)
//                                {
//                                    glidePlaneNoisePointsValuesXZ->SetNumberOfValues(slipSystem->planeNoise->gridSize.array().prod());
//                                    //                    glidePlaneNoisePointsValuesXZ->SetName(glidePlanesNoiseBox->itemText(0).toStdString().c_str());
//
//                                    glidePlaneNoisePointsValuesYZ->SetNumberOfValues(slipSystem->planeNoise->gridSize.array().prod());
//                                    //                    glidePlaneNoisePointsValuesYZ->SetName(glidePlanesNoiseBox->itemText(1).toStdString().c_str());
//                                }
//
                                vtkNew<vtkDoubleArray> ssNoise; // solid solution
                                if(slipSystem->planeNoise->solidSolution)
                                {
                                    ssNoise->SetNumberOfValues(slipSystem->planeNoise->gridSize.array().prod());
                                    ssNoise->SetName(glidePlanesNoiseBox->itemText(0).toStdString().c_str());
                                }

                                
                                vtkNew<vtkDoubleArray> sfNoise; // stacking-fault noise
                                if(slipSystem->planeNoise->stackingFault)
                                {
                                    sfNoise->SetNumberOfValues(slipSystem->planeNoise->gridSize.array().prod());
                                    sfNoise->SetName(glidePlanesNoiseBox->itemText(1).toStdString().c_str());
                                }
//
//
//
//
//                                for (int k = 0; k < slipSystem->planeNoise->gridSize(2); k++)
//                                {
//                                    for (int j = 0; j < slipSystem->planeNoise->gridSize(1); j++)
//                                    {
//                                        for (int i = 0; i < slipSystem->planeNoise->gridSize(0); i++)
//                                        {
//
//                                            const double x(i*slipSystem->planeNoise->gridSpacing(0));
//                                            const double y(j*slipSystem->planeNoise->gridSpacing(1));
//
//                                            const Eigen::Matrix<double,2,1> localPos((Eigen::Matrix<double,2,1>()<<x,y).finished());
//                                            const Eigen::Matrix<double,3,1> globalPos(slipSystem->localToGlobal(localPos)+glidePlane->P);
//
//                                            glidePlaneNoisePoints->InsertNextPoint(globalPos(0), globalPos(1), globalPos(2));
//
////                                            const auto ind = slipSystem->planeNoise->gridSize(1)*slipSystem->planeNoise->gridSize(2)*i + j*slipSystem->planeNoise->gridSize(2) + k;
//
//                                            const auto noiseVal(slipSystem->gridInterp(localPos));
//                                            
//                                            if(slipSystem->planeNoise->solidSolution)
//                                            {
////                                                const auto& noiseVal(slipSystem->planeNoise->solidSolution->operator[](ind));
//                                                glidePlaneNoisePointsValuesXZ->SetValue(ind, noiseVal(0));
//                                                glidePlaneNoisePointsValuesYZ->SetValue(ind, noiseVal(1));
//                                                noiseValMin=std::min(noiseValMin,std::min(noiseVal(0),noiseVal(1)));
//                                                noiseValMax=std::max(noiseValMax,std::max(noiseVal(0),noiseVal(1)));
//                                            }
//
//                                            if(slipSystem->planeNoise->stackingFault)
//                                            {
//                                                const auto& noiseVal(slipSystem->planeNoise->stackingFault->operator[](ind));
//                                                sfNoise->SetValue(ind, noiseVal);
//                                                noiseValMin=std::min(noiseValMin,noiseVal);
//                                                noiseValMax=std::max(noiseValMax,noiseVal);
//                                            }
//
//
//                                        }
//                                    }
//                                }
//                                vtkNew<vtkStructuredGrid> glidePlaneNoiseGrid;
//                                glidePlaneNoiseGrid->SetDimensions(slipSystem->planeNoise->gridSize(0),slipSystem->planeNoise->gridSize(1),slipSystem->planeNoise->gridSize(2));
//                                glidePlaneNoiseGrid->SetPoints(glidePlaneNoisePoints);
//                                if(slipSystem->planeNoise->solidSolution)
//                                {
//                                    glidePlaneNoiseGrid->GetPointData()->AddArray(glidePlaneNoisePointsValuesXZ);
//                                    glidePlaneNoiseGrid->GetPointData()->AddArray(glidePlaneNoisePointsValuesYZ);
//                                }
//                                if(slipSystem->planeNoise->stackingFault)
//                                {
//                                    glidePlaneNoiseGrid->GetPointData()->AddArray(sfNoise);
//                                }
//
//
//
//
//                                vtkNew<vtkLookupTable> glidePlaneNoiseLut;
//                                glidePlaneNoiseLut->SetNumberOfTableValues(slipSystem->planeNoise->gridSize.array().prod());
//                                glidePlaneNoiseLut->Build();
//
//                                noiseMappers.push_back(vtkSmartPointer<vtkDataSetMapper>::New());
//                                noiseMappers.back()->SetInputData(glidePlaneNoiseGrid);
//                                noiseMappers.back()->SetScalarModeToUsePointFieldData();
//                                noiseMappers.back()->SelectColorArray(glidePlanesNoiseBox->currentIndex());
//                                noiseMappers.back()->SetLookupTable(glidePlaneNoiseLut);
//                                noiseMappers.back()->SetScalarRange(noiseValMin, noiseValMax);
//                                noiseMappers.back()->ScalarVisibilityOn();
//                                noiseActors.push_back(vtkSmartPointer<vtkActor>::New());
//                                noiseActors.back()->SetMapper(noiseMappers.back());
//                                noiseActors.back()->SetVisibility(showGlidePlanesNoise->isChecked());
//                                renderer->AddActor( noiseActors.back());
                            }
                        }
                    }
                }
                
                
            }
        }
        
        glidePlanePolydata->SetPoints(glidePlanePoints);
        glidePlanePolydata->SetLines(glidePlaneCells);
        glidePlanePolydata->Modified();
        
        std::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
        
    }

    void GlidePlaneActor::modify()
    {
        glidePlaneActor->SetVisibility(showGlidePlanes->isChecked());
        
        const size_t slipSystemID(slipSystemNoiseBox->currentIndex());
        glidePlanesNoiseBox->setEnabled(showGlidePlanesNoise->isChecked());
        for(const auto& mapper : noiseMappers[slipSystemID])
        {
            mapper->SelectColorArray(glidePlanesNoiseBox->currentIndex());
        }
        
        
        for(const auto& actorPair : noiseActors)
        {
            for(const auto& actor : actorPair.second)
            {
                actor->SetVisibility(showGlidePlanesNoise->isChecked());
            }
        }
        
        
        renderWindow->Render();
    }


} // namespace model
#endif
