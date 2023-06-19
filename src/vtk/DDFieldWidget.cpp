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

DDFieldWidget::DDFieldWidget(vtkGenericOpenGLRenderWindow* const renWin_in,
                             vtkRenderer* const renderer_in,
                             const SimplicialMesh<3>& mesh_in):
/* init */ mainLayout(new QGridLayout(this))
/* init */,boxLabel(new QLabel(tr("# of planes")))
/* init */,spinBox(new QSpinBox(this))
/* init */,computeButton(new QPushButton(tr("Compute")))
/* init */,groupBox(new QGroupBox(tr("&Planes")))
/* init */,renWin(renWin_in)
/* init */,renderer(renderer_in)
/* init */,mesh(mesh_in)
{
    
    QVBoxLayout* groupBoxLayout = new QVBoxLayout();
    groupBox->setLayout(groupBoxLayout);

    mainLayout->addWidget(boxLabel,0,0,1,1);
    mainLayout->addWidget(spinBox,0,1,1,1);
    mainLayout->addWidget(computeButton,1,0,1,2);
    mainLayout->addWidget(groupBox,2,0,1,2);

    this->setLayout(mainLayout);
    
    connect(spinBox,SIGNAL(valueChanged(int)), this, SLOT(resetPlanes()));
    connect(computeButton,SIGNAL(released()), this, SLOT(compute()));

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
                    ddPlaneField->compute();
                }
            }
        }
    }
}

void DDPlaneField::compute()
{
    std::cout<<"DDPlaneField computing "<<std::endl;
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
        DDPlaneField* planeField(new DDPlaneField(renWin,renderer,mesh));
        planeField->resetPlane();
        groupBox->layout()->addWidget(planeField);
    }
    renWin->Render();
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
                           const SimplicialMesh<3>& mesh_in):
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
/* init */,mesh(mesh_in)
{
    const VectorDim c(0.5*(mesh.xMax()+mesh.xMin()));
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
                        plane.reset(new MeshPlane<3>(mesh,P,N));
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

FieldDataPnt::FieldDataPnt(const Eigen::Matrix<double,3,1>& Pin) :
    /* init */ P(Pin)
    /* init */,solidAngle(0.0)
{

}

} // namespace model
#endif
