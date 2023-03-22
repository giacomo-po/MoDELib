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
#include <vtkPlanes.h>
#include <vtkStructuredGrid.h>
#include <vtkUnstructuredGrid.h>
#include <vtkLookupTable.h>
#include <vtkStructuredGridAppend.h>
#include <vtkTableBasedClipDataSet.h>

#include <GlidePlaneActor.h>
#include <Polygon2D.h>

namespace model
{

    GlidePlaneActor::GlidePlaneActor(vtkGenericOpenGLRenderWindow* const renWin,vtkRenderer* const ren,const Polycrystal<3>& poly_in,const DDtraitsIO& ) :
    /* init */ renderWindow(renWin)
    /* init */,renderer(ren)
    /* init */,poly(poly_in)
    /* init */,mainLayout(new QGridLayout(this))
    /* init */,showGlidePlanes(new QCheckBox(this))
    /* init */,showGlidePlanesNoise(new QCheckBox(this))
    /* init */,glidePlanesNoiseBox(new QComboBox(this))
    /* init */,slipSystemNoiseBox(new QComboBox(this))
    /* init */,ssNoiseMin(new QLineEdit("0.0"))
    /* init */,ssNoiseMax(new QLineEdit("0.0"))
    /* init */,sfNoiseMin(new QLineEdit("0.0"))
    /* init */,sfNoiseMax(new QLineEdit("0.0"))
    /* init */,glidePlanePolydata(vtkSmartPointer<vtkPolyData>::New())
    /* init */,glidePlaneMapper(vtkSmartPointer<vtkPolyDataMapper>::New())
    /* init */,glidePlaneActor(vtkSmartPointer<vtkActor>::New())
    /* init */,glidePlaneFactory(poly)
    /* init */,noiseLimits(Eigen::Array<double,2,2>::Zero())
//    /* init */,planeNoise(traitsIO.polyFile,poly)
    {
        
        showGlidePlanes->setChecked(false);
        showGlidePlanes->setText("show GlidePlanes");
        glidePlaneActor->SetVisibility(false);
        
        showGlidePlanesNoise->setChecked(false);
        showGlidePlanesNoise->setText("show GlidePlanesNoise");
        glidePlanesNoiseBox->setEnabled(false);
        slipSystemNoiseBox->setEnabled(false);

            glidePlanesNoiseBox->addItem("RSS");
            glidePlanesNoiseBox->addItem("ISF energy");
            const auto& grain(poly.grains.begin()->second);
                    
        for(size_t k=0; k<grain.singleCrystal->slipSystems().size();++k)
        {
            slipSystemNoiseBox->addItem(QString::fromStdString("slipSystem "+ std::to_string(k)));
        }
        
        mainLayout->addWidget(showGlidePlanes,0,0,1,1);
        mainLayout->addWidget(showGlidePlanesNoise,1,0,1,1);
        mainLayout->addWidget(glidePlanesNoiseBox,1,1,1,1);
        mainLayout->addWidget(slipSystemNoiseBox,2,1,1,1);
        mainLayout->addWidget(ssNoiseMin,3,0,1,1);
        mainLayout->addWidget(ssNoiseMax,3,1,1,1);
        mainLayout->addWidget(sfNoiseMin,4,0,1,1);
        mainLayout->addWidget(sfNoiseMax,4,1,1,1);

        
        this->setLayout(mainLayout);
        connect(showGlidePlanes,SIGNAL(stateChanged(int)), this, SLOT(modify()));
        connect(showGlidePlanesNoise,SIGNAL(stateChanged(int)), this, SLOT(modify()));
        connect(glidePlanesNoiseBox,SIGNAL(currentIndexChanged(int)), this, SLOT(modify()));
        connect(slipSystemNoiseBox,SIGNAL(currentIndexChanged(int)), this, SLOT(modify()));
        connect(ssNoiseMin,SIGNAL(returnPressed()), this, SLOT(modify()));
        connect(ssNoiseMax,SIGNAL(returnPressed()), this, SLOT(modify()));
        connect(sfNoiseMin,SIGNAL(returnPressed()), this, SLOT(modify()));
        connect(sfNoiseMax,SIGNAL(returnPressed()), this, SLOT(modify()));

        
        
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
//        double noiseSSMin(std::numeric_limits<double>::max());
//        double noiseSSMax(std::numeric_limits<double>::min());
//        double noiseSSMin(std::numeric_limits<double>::max());
//        double noiseSSMax(std::numeric_limits<double>::min());

        noiseLimits<<std::numeric_limits<double>::max(),std::numeric_limits<double>::min(),
        /*         */std::numeric_limits<double>::max(),std::numeric_limits<double>::min();
        
        for(const auto& gpIO : auxIO.glidePlanes())
        {
            // Plot GlidePlaneBoundary
            GlidePlaneKey<3> planeKey(gpIO.r,gpIO.h,gpIO.latticeID);
            const auto glidePlane(glidePlaneFactory.getFromKey(planeKey));
            if(glidePlane->slipSystems().size())
            {
                
                for(const auto& mshInt : glidePlane->meshIntersections)
                {
                    glidePlanePoints->InsertNextPoint(mshInt->P0.data());
                    
                    glidePlanePoints->InsertNextPoint(mshInt->P1.data());
                    
                    vtkSmartPointer<vtkLine> line(vtkSmartPointer<vtkLine>::New());
                    line->GetPointIds()->SetId(0, pointIDs); // the second 0 is the index of the Origin in linesPolyData's points
                    line->GetPointIds()->SetId(1, pointIDs+1);
                    glidePlaneCells->InsertNextCell(line);
                    pointIDs+=2;
                }
                
                
                auto clipNormals = vtkSmartPointer<vtkDoubleArray>::New();
                clipNormals->SetNumberOfComponents(3);
                vtkNew<vtkPoints> clipPts;
                for(const auto& meshInt : glidePlane->meshIntersections)
                {
                    clipNormals->InsertNextTuple(meshInt->outNormal.data());
                    clipPts->InsertNextPoint(meshInt->P0.data());
                }
                vtkNew<vtkPlanes> clipPlanes;
                clipPlanes->SetNormals(clipNormals);
                clipPlanes->SetPoints(clipPts);
                
                // Plot solidSolution-noise
                for(size_t slipSystemID=0; slipSystemID<poly.grain(gpIO.latticeID).singleCrystal->slipSystems().size();++slipSystemID)
                {
                    const auto& slipSystem(poly.grain(gpIO.latticeID).singleCrystal->slipSystems()[slipSystemID]);
                    if(slipSystem->planeNoise)
                    {
                        if(slipSystem->planeNoise->solidSolution || slipSystem->planeNoise->stackingFault)
                        {
                            if(glidePlane->n.cross(slipSystem->n).squaredNorm()==0)
                            {// glide plane includes slip system
                                
                                std::set<double> xPos;
                                std::set<double> yPos;
//                                std::vector<Eigen::Array<double,2,1>> bndPoly;
                                
//                                auto clipNormals = vtkSmartPointer<vtkDoubleArray>::New();
//                                clipNormals->SetNumberOfComponents(3);
//                                vtkNew<vtkPoints> clipPts;
                                for(const auto& meshInt : glidePlane->meshIntersections)
                                {
//                                    clipNormals->InsertNextTuple(meshInt->outNormal.data());
//                                    clipPts->InsertNextPoint(meshInt->P0.data());

                                    const auto localPos(slipSystem->globalToLocal(meshInt->P0-glidePlane->P));
 //                                   bndPoly.push_back(localPos);
                                    xPos.insert(localPos(0));
                                    yPos.insert(localPos(1));
                                }

//                                vtkNew<vtkPlanes> clipPlanes;
//                                clipPlanes->SetNormals(clipNormals);
//                                clipPlanes->SetPoints(clipPts);
                                
                                const Eigen::Array<double,2,1> lCorner(*xPos. begin(),*yPos. begin());
                                const Eigen::Array<double,2,1> hCorner(*xPos.rbegin(),*yPos.rbegin());

                                const Eigen::Array<int,2,1> lIdx(slipSystem->planeNoise->posToIdx(lCorner).first);
                                const Eigen::Array<int,2,1> hIdx(slipSystem->planeNoise->posToIdx(hCorner).second);
                                
                                const int Nx(hIdx(0)-lIdx(0)+1);
                                const int Ny(hIdx(1)-lIdx(1)+1);
                                
                                vtkNew<vtkPoints> glidePlaneNoisePoints;
                                glidePlaneNoisePoints->SetNumberOfPoints(Nx*Ny);
                                
                                vtkNew<vtkDoubleArray> glidePlaneSSNoise;
                                glidePlaneSSNoise->SetNumberOfValues(Nx*Ny);

                                vtkNew<vtkDoubleArray> glidePlaneSFNoise;
                                glidePlaneSFNoise->SetNumberOfValues(Nx*Ny);

                                
                                for(int i=0;i<Nx;++i)
                                {
                                    const int gridi(i+lIdx(0));
                                    
                                    for(int j=0;j<Ny;++j)
                                    {
                                        const int gridj(j+lIdx(1));
//                                        const int storageIdx = Ny*i + j;
                                        const int storageIdx = Nx*j + i;

                                        const Eigen::Array<int,2,1> idx(gridi,gridj);
                                        const Eigen::Array<double,2,1> localPos(slipSystem->planeNoise->idxToPos(idx));

                                        const Eigen::Array<double,3,1> globalPos(slipSystem->localToGlobal(localPos.matrix())+glidePlane->P);
                                        glidePlaneNoisePoints->SetPoint(storageIdx,globalPos.data());
                                        
                                        const auto noiseVal(slipSystem->gridVal(idx));
                                        glidePlaneSSNoise->SetValue(storageIdx,std::get<1>(noiseVal));
                                        glidePlaneSFNoise->SetValue(storageIdx,std::get<2>(noiseVal));
                                        noiseLimits<<std::min(noiseLimits(0,0),std::get<1>(noiseVal)),std::max(noiseLimits(0,1),std::get<1>(noiseVal)),
                                        /*         */std::min(noiseLimits(1,0),std::get<2>(noiseVal)),std::max(noiseLimits(1,1),std::get<2>(noiseVal));
                                        
                                    }
                                }
                                

                                ssNoiseMin->setText(QString::number(noiseLimits(0,0)));
                                ssNoiseMax->setText(QString::number(noiseLimits(0,1)));
                                sfNoiseMin->setText(QString::number(noiseLimits(1,0)));
                                sfNoiseMax->setText(QString::number(noiseLimits(1,1)));

                                
                                vtkNew<vtkStructuredGrid> glidePlaneNoiseGrid;
                                
                                glidePlaneNoiseGrid->SetDimensions(Nx,Ny,1);
                                glidePlaneNoiseGrid->SetPoints(glidePlaneNoisePoints);
                                glidePlaneNoiseGrid->GetPointData()->AddArray(glidePlaneSSNoise);
                                glidePlaneNoiseGrid->GetPointData()->AddArray(glidePlaneSFNoise);
                                
                                
                                vtkNew<vtkTableBasedClipDataSet> clipper;
                                clipper->SetInputData(glidePlaneNoiseGrid);
                                clipper->SetClipFunction(clipPlanes);
                                clipper->InsideOutOn();

//                                vtkNew<vtkLookupTable> glidePlaneNoiseLut;
//                                glidePlaneNoiseLut->SetNumberOfTableValues(Nx*Ny);
//                                glidePlaneNoiseLut->Build();

                                noiseMappers[slipSystemID].push_back(vtkSmartPointer<vtkDataSetMapper>::New());
//                                noiseMappers[slipSystemID].back()->SetInputData(glidePlaneNoiseGrid);
                                noiseMappers[slipSystemID].back()->SetInputConnection(clipper->GetOutputPort());
                                noiseMappers[slipSystemID].back()->SetScalarModeToUsePointFieldData();
                                noiseMappers[slipSystemID].back()->SelectColorArray(glidePlanesNoiseBox->currentIndex());
//                                noiseMappers[slipSystemID].back()->SetLookupTable(glidePlaneNoiseLut);
                                noiseMappers[slipSystemID].back()->SetScalarRange(noiseLimits(glidePlanesNoiseBox->currentIndex(),0), noiseLimits(glidePlanesNoiseBox->currentIndex(),1));
                                noiseMappers[slipSystemID].back()->ScalarVisibilityOn();
                                noiseActors[slipSystemID].push_back(vtkSmartPointer<vtkActor>::New());
                                noiseActors[slipSystemID].back()->SetMapper(noiseMappers[slipSystemID].back());
                                noiseActors[slipSystemID].back()->SetVisibility(showGlidePlanesNoise->isChecked() && slipSystemNoiseBox->currentIndex()==int(slipSystemID));
                                renderer->AddActor( noiseActors[slipSystemID].back());

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
        
        glidePlanesNoiseBox->setEnabled(showGlidePlanesNoise->isChecked());
        slipSystemNoiseBox->setEnabled(showGlidePlanesNoise->isChecked());

        // Select solidSolution or stacking-fault noise
        for(const auto& mapperPair : noiseMappers)
        {
            for(const auto& mapper : mapperPair.second)
            {
                mapper->SelectColorArray(glidePlanesNoiseBox->currentIndex());
                
                switch (glidePlanesNoiseBox->currentIndex())
                {
                    case 0:
                        mapper->SetScalarRange(ssNoiseMin->text().toDouble(), ssNoiseMax->text().toDouble());
                        break;
                        
                    case 1:
                        mapper->SetScalarRange(sfNoiseMin->text().toDouble(), sfNoiseMax->text().toDouble());
                        break;
                        
                    default:
                        break;
                }
                

            }
        }
        
        const size_t slipSystemID(slipSystemNoiseBox->currentIndex());
        for(const auto& actorPair : noiseActors)
        {
            for(const auto& actor : actorPair.second)
            {
                actor->SetVisibility(actorPair.first==slipSystemID && showGlidePlanesNoise->isChecked());
            }
        }
        
        
        renderWindow->Render();
    }


} // namespace model
#endif
