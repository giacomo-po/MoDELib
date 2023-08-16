/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DDFieldWidget_cpp_
#define model_DDFieldWidget_cpp_



#include <QVBoxLayout>
#include <vtkProperty.h>
#include <vtkTriangle.h>
#include <vtkCellData.h>

#include <DDFieldWidget.h>

namespace model
{

FieldDataPnt::FieldDataPnt(const Eigen::Matrix<double,3,1>& Pin):
/* init */ P(Pin)
/* init */,solidAngle(0.0)
/* init */,stressDD(Eigen::Matrix<double,3,3>::Zero())
/* init */,stressIN(Eigen::Matrix<double,3,3>::Zero())
{
    
}

double FieldDataPnt::value(const int& valID,const bool& useDD,const bool& useIN) const
{
    switch (valID)
    {
        case 0:
            return stressDD(0,0)*useDD+stressIN(0,0)*useIN;
            break;
        case 1:
            return stressDD(0,1)*useDD+stressIN(0,1)*useIN;
            break;
        case 2:
            return stressDD(0,2)*useDD+stressIN(0,2)*useIN;
            break;
        case 3:
            return stressDD(1,1)*useDD+stressIN(1,1)*useIN;
            break;
        case 4:
            return stressDD(1,2)*useDD+stressIN(1,2)*useIN;
            break;
        case 5:
            return stressDD(2,2)*useDD+stressIN(2,2)*useIN;
            break;
        case 6:
            return value(0,useDD,useIN)+value(3,useDD,useIN)+value(5,useDD,useIN);
            break;
        case 7:
        {
            const Eigen::Matrix<double,3,3> stress(stressDD*useDD+stressIN*useIN);
            const Eigen::Matrix<double,3,3> stressDev(stress-stress.trace()/3.0*Eigen::Matrix<double,3,3>::Identity());
            return std::sqrt((stressDev*stressDev).trace()*1.5);
            break;
        }
        case 8:
            return solidAngle;
            break;

        default:
            return 0.0;
            break;
    }
}


DDFieldWidget::DDFieldWidget(vtkGenericOpenGLRenderWindow* const renWin_in,
                             vtkRenderer* const renderer_in,
                             const  DDconfigFields<3>& configFields_in):
/* init */ mainLayout(new QGridLayout(this))
/* init */,boxLabel(new QLabel(tr("# of planes")))
/* init */,spinBox(new QSpinBox(this))
/* init */,groupBox(new QGroupBox(tr("&Planes")))
/* init */,computeButton(new QPushButton(tr("Compute")))
/* init */,fieldComboBox(new QComboBox(this))
/* init */,dislocationsCheck(new QCheckBox("dislocations",this))
/* init */,inclusionsCheck(new QCheckBox("inclusions",this))
/* init */,customScaleBox(new QGroupBox(tr("&Custom scale")))
/* init */,minScale(new QLineEdit(tr("0")))
/* init */,maxScale(new QLineEdit(tr("1")))
/* init */,lut(vtkSmartPointer<vtkLookupTable>::New())
/* init */,scalarBar(vtkSmartPointer<vtkScalarBarActor>::New())
/* init */,renWin(renWin_in)
/* init */,renderer(renderer_in)
/* init */,configFields(configFields_in)
{
    
    QVBoxLayout* groupBoxLayout = new QVBoxLayout();
    groupBox->setLayout(groupBoxLayout);

    QGridLayout* autoscaleLayout = new QGridLayout();
    customScaleBox->setCheckable(true);
    customScaleBox->setChecked(false);

    autoscaleLayout->addWidget(minScale,0,0,1,1);
    autoscaleLayout->addWidget(maxScale,0,1,1,1);
    customScaleBox->setLayout(autoscaleLayout);


    fieldComboBox->insertItem(0,"stress_11");
    fieldComboBox->insertItem(1,"stress_12");
    fieldComboBox->insertItem(2,"stress_13");
//    fieldComboBox->insertItem(3,"stress_21");
    fieldComboBox->insertItem(3,"stress_22");
    fieldComboBox->insertItem(4,"stress_23");
//    fieldComboBox->insertItem(6,"stress_31");
//    fieldComboBox->insertItem(7,"stress_32");
    fieldComboBox->insertItem(5,"stress_33");
    fieldComboBox->insertItem(6,"tr(stress)");
    fieldComboBox->insertItem(7,"stress_VM");
    fieldComboBox->insertItem(8,"solid angle");

    dislocationsCheck->setChecked(true);
    inclusionsCheck->setChecked(true);

    
    connect(spinBox,SIGNAL(valueChanged(int)), this, SLOT(resetPlanes()));
    connect(computeButton,SIGNAL(released()), this, SLOT(compute()));
    connect(fieldComboBox,SIGNAL(currentIndexChanged(int)), this, SLOT(plotField()));
    connect(minScale,SIGNAL(returnPressed()), this, SLOT(plotField()));
    connect(maxScale,SIGNAL(returnPressed()), this, SLOT(plotField()));
    connect(dislocationsCheck,SIGNAL(stateChanged(int)), this, SLOT(plotField()));
    connect(inclusionsCheck,SIGNAL(stateChanged(int)), this, SLOT(plotField()));
    connect(customScaleBox,SIGNAL(toggled(bool)), this, SLOT(plotField()));

    
    
    mainLayout->addWidget(boxLabel,0,0,1,1);
    mainLayout->addWidget(spinBox,0,1,1,1);
    mainLayout->addWidget(groupBox,1,0,1,2);
    mainLayout->addWidget(computeButton,2,0,1,1);
    mainLayout->addWidget(fieldComboBox,3,0,1,1);
    mainLayout->addWidget(dislocationsCheck,2,1,1,1);
    mainLayout->addWidget(inclusionsCheck,3,1,1,1);
    mainLayout->addWidget(customScaleBox,4,0,1,2);
    this->setLayout(mainLayout);

    
    scalarBar->VisibilityOff();
    scalarBar->SetNumberOfLabels(4);
    scalarBar->GetLabelTextProperty()->SetColor(0,0,0);
    lut->SetHueRange(0.66667, 0.0);
    scalarBar->SetLookupTable( lut );
    lut->Build();


}

void DDFieldWidget::compute()
{
    if(groupBox->layout())
    {
        for(int k=0;k<groupBox->layout()->count();++k)
        {
            QLayoutItem *item = groupBox->layout()->itemAt(k);
            QWidget* widget = item->widget();
            if(widget)
            {
                auto* ddPlaneField = dynamic_cast<DDPlaneField*>(widget);
                if (ddPlaneField)
                {
                    ddPlaneField->compute(configFields);
                }
            }
        }
    }
    plotField();
}

void DDFieldWidget::plotField()
{
    if(groupBox->layout())
    {
        const int valID(fieldComboBox->currentIndex());

        // Compute lut range if customScaleBox is off
        if(!customScaleBox->isChecked())
        {// compute range
            double minValue=std::numeric_limits<double>::max();
            double maxValue=-std::numeric_limits<double>::max();
            for(int k=0;k<groupBox->layout()->count();++k)
            {
                QLayoutItem *item = groupBox->layout()->itemAt(k);
                QWidget* widget = item->widget();
                if(widget)
                {
                    auto* ddPlaneField = dynamic_cast<DDPlaneField*>(widget);
                    if (ddPlaneField)
                    {
                        if(ddPlaneField->groupBox->isChecked())
                        {
                            for(const auto& vtx : ddPlaneField->dataPnts())
                            {
                                const double value(vtx.value(valID,dislocationsCheck->isChecked(),inclusionsCheck->isChecked()));
                                minValue=std::min(minValue,value);
                                maxValue=std::max(maxValue,value);
                            }
                        }
                    }
                }
            }
            minScale->setText(QString::fromStdString(std::to_string(minValue)));
            maxScale->setText(QString::fromStdString(std::to_string(maxValue)));
        }
        
        // Set lut range
        try
        {
            double lutMin(std::stod(minScale->text().toStdString().c_str()));
//            minScale->setStyleSheet("color: black");
            try
            {
                double lutMax(std::stod(maxScale->text().toStdString().c_str()));
//                maxScale->setStyleSheet("color: black");
                lut->SetTableRange(lutMin, lutMax);
                
                for(int k=0;k<groupBox->layout()->count();++k)
                {
                    QLayoutItem *item = groupBox->layout()->itemAt(k);
                    QWidget* widget = item->widget();
                    if(widget)
                    {
                        auto* ddPlaneField = dynamic_cast<DDPlaneField*>(widget);
                        if (ddPlaneField)
                        {
                            if(ddPlaneField->groupBox->isChecked())
                            {
                                ddPlaneField->plotField(valID,dislocationsCheck->isChecked(),inclusionsCheck->isChecked(),lut);
                            }
                        }
                    }
                }
                renWin->Render();
            }
            catch (const std::exception& e)
            {
//                maxScale->setStyleSheet("color: red");
            }
        }
        catch (const std::exception& e)
        {
//            minScale->setStyleSheet("color: red");
        }
    }
}

void DDPlaneField::compute(const DDconfigFields<3>& configFields)
{
    std::cout<<"DDPlaneField computing "<<dataPnts().size()<<" points..."<<std::flush;
    const auto t0= std::chrono::system_clock::now();
    for(size_t vtkID=0;vtkID<dataPnts().size();++vtkID)
    {
        auto& vtx(dataPnts()[vtkID]);
        vtx.solidAngle=configFields.solidAngle(vtx.P);
        vtx.stressDD=configFields.dislocationStress(vtx.P);
        vtx.stressIN=configFields.inclusionStress(vtx.P);
    }
    std::cout<<magentaColor<<"["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
}

void DDFieldWidget::clearLayout(QLayout *layout)
{
    if(layout)
    {
        while(layout->count() > 0)
        {
            QLayoutItem *item = layout->takeAt(0);
            QWidget* widget = item->widget();
            if(widget)
            {
                delete widget;
            }
            delete item;
        }
    }
}

void DDFieldWidget::resetPlanes()
{
    clearLayout(groupBox->layout());
    std::cout<<"Resetting "<<spinBox->value()<<" planes"<<std::endl;
    for(int k=0;k<spinBox->value();++k)
    {
        DDPlaneField* planeField(new DDPlaneField(renWin,renderer,configFields.ddBase.poly));
        planeField->resetPlane();
        groupBox->layout()->addWidget(planeField);
    }
    renWin->Render();
}

void DDPlaneField::plotField(const int& valID,const bool& useDD,const bool& useIN,const vtkSmartPointer<vtkLookupTable>& lut)
{
    vtkSmartPointer<vtkUnsignedCharArray> fieldColors(vtkSmartPointer<vtkUnsignedCharArray>::New());
    fieldColors->SetNumberOfComponents(3);
    for(auto& vtx : dataPnts())
    {
        const double value(vtx.value(valID,useDD,useIN));
        double dclr[3];
        lut->GetColor(value, dclr);
        unsigned char cclr[3];
        for(unsigned int j = 0; j < 3; j++)
        {
            cclr[j] = static_cast<unsigned char>(255.0 * dclr[j]);
        }
        fieldColors->InsertNextTypedTuple(cclr);
    }
    meshPolydata->GetPointData()->SetScalars(fieldColors);
    meshPolydata->Modified();
    meshMapper->SetScalarModeToUsePointData();
}

const std::deque<FieldDataPnt>& DDPlaneField::dataPnts() const
{
    return *this;
}

std::deque<FieldDataPnt>& DDPlaneField::dataPnts()
{
    return *this;
}

DDPlaneField::DDPlaneField(vtkGenericOpenGLRenderWindow* const renWin_in,
                           vtkRenderer* const renderer_in,
                           const Polycrystal<3>& poly_in):
/* init */ mainLayout(new QGridLayout(this))
/* init */,groupBox(new QGroupBox(tr("&Plane")))
/* init */,posEdit(new QLineEdit("0 0 0"))
/* init */,normalEdit(new QLineEdit("1 1 1"))
/* init */,meshSizeEdit(new QLineEdit("1000"))
/* init */,meshPolydata(vtkSmartPointer<vtkPolyData>::New())
/* init */,meshMapper(vtkSmartPointer<vtkPolyDataMapper>::New())
/* init */,meshActor(vtkSmartPointer<vtkActor>::New())
/* init */,renWin(renWin_in)
/* init */,renderer(renderer_in)
/* init */,poly(poly_in)
{
    const VectorDim c(0.5*(poly.mesh.xMax()+poly.mesh.xMin()));
    posEdit->setText(QString::fromStdString(std::to_string(c(0))+" "+std::to_string(c(1))+" "+std::to_string(c(2))));
    
    const VectorDim n(Eigen::Matrix<double,3,1>::Random());
    normalEdit->setText(QString::fromStdString(std::to_string(n(0))+" "+std::to_string(n(1))+" "+std::to_string(n(2))));

    QGridLayout* groupBoxLayout = new QGridLayout();
    groupBoxLayout->addWidget(posEdit,0,0,1,1);
    groupBoxLayout->addWidget(normalEdit,1,0,1,1);
    groupBoxLayout->addWidget(meshSizeEdit,2,0,1,1);
    groupBox->setLayout(groupBoxLayout);

    groupBox->setCheckable(true);
    mainLayout->addWidget(groupBox,0,0,1,1);
    this->setLayout(mainLayout);

    connect(posEdit,SIGNAL(returnPressed()), this, SLOT(resetPlane()));
    connect(normalEdit,SIGNAL(returnPressed()), this, SLOT(resetPlane()));
    connect(meshSizeEdit,SIGNAL(returnPressed()), this, SLOT(resetPlane()));
    connect(groupBox,SIGNAL(toggled(bool)), this, SLOT(modify()));

    meshPolydata->Allocate();
    meshMapper->SetInputData(meshPolydata);
    meshActor->SetMapper ( meshMapper );
    meshActor->GetProperty()->SetOpacity(0.8); //Make the mesh have some transparency.
    renderer->AddActor(meshActor);
}

DDPlaneField::~DDPlaneField()
{
    renderer->RemoveActor(meshActor);
}

void DDPlaneField::modify()
{
    meshActor->SetVisibility(groupBox->isChecked());
    renWin->Render();
}

void DDPlaneField::resetPlane()
{
    
    double meshSize;
    std::stringstream ssM(meshSizeEdit->text().toStdString());
    if(ssM >> meshSize)
    {
        meshSizeEdit->setStyleSheet("background-color: white");

        VectorDim P;
            std::stringstream ssP(posEdit->text().toStdString());
            if(ssP >> P(0) && ssP >> P(1) && ssP >> P(2))
            {
                posEdit->setStyleSheet("background-color: white");
                std::cout<<"P="<<P.transpose()<<std::endl;
                
                VectorDim N;
                std::stringstream ssN(normalEdit->text().toStdString());
                if(ssN >> N(0) && ssN >> N(1) && ssN >> N(2))
                {
                    normalEdit->setStyleSheet("background-color: white");
                    std::cout<<"N="<<N.transpose()<<std::endl;
                    const double nNorm(N.norm());
                    if(nNorm>FLT_EPSILON)
                    {
                        plane.reset(new MeshPlane<3>(poly.mesh,P,N));
                        std::deque<Eigen::Matrix<double,2,1>> boundaryPts;
                        std::deque<Eigen::Matrix<double,2,1>> internalPts;
                        for(const auto& bndLine : plane->meshIntersections)
                        {
                            boundaryPts.push_back(plane->localPosition(bndLine->P0));

                        }
                        this->reMesh(boundaryPts,internalPts,meshSize);
                        
                        dataPnts().clear();
                        vtkSmartPointer<vtkPoints> meshPts(vtkSmartPointer<vtkPoints>::New());
                        for(const auto& point2d : this->vertices())
                        {
                            const auto point3d(plane->globalPosition(point2d));
                            meshPts->InsertNextPoint(point3d(0),point3d(1),point3d(2));
                            dataPnts().emplace_back(point3d);
                        }
                        
                        vtkSmartPointer<vtkCellArray> meshTriangles(vtkSmartPointer<vtkCellArray>::New());
                        vtkSmartPointer<vtkUnsignedCharArray> meshColors(vtkSmartPointer<vtkUnsignedCharArray>::New());
                        meshColors->SetNumberOfComponents(3);
                        for(const auto& tri : this->triangles())
                        {
                            vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
                            triangle->GetPointIds()->SetId (0,tri(0));
                            triangle->GetPointIds()->SetId (1,tri(1));
                            triangle->GetPointIds()->SetId (2,tri(2));
                            meshTriangles->InsertNextCell ( triangle );
                            const auto triColor(Eigen::Matrix<int,1,3>::Random()*255);
                            meshColors->InsertNextTuple3(triColor(0),triColor(1),triColor(2)); // use this to assig color to each vertex
                        }
                        
                        meshPolydata->SetPoints ( meshPts );
                        meshPolydata->SetPolys ( meshTriangles );
                        meshPolydata->GetCellData()->SetScalars(meshColors);
                        meshPolydata->Modified();
                        meshMapper->SetScalarModeToUseCellData();
                        renWin->Render();
                    }
                    else
                    {
                        normalEdit->setStyleSheet("background-color: red");
                    }
                }
                else
                {
                    normalEdit->setStyleSheet("background-color: red");
                }
            }
        else
        {
            posEdit->setStyleSheet("background-color: red");
        }
    }
    else
    {
        meshSizeEdit->setStyleSheet("background-color: red");
    }
}


} // namespace model
#endif
