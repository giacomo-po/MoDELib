/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_QuadratureActor_cpp_
#define model_QuadratureActor_cpp_

#include <vtkLine.h>
#include <vtkStructuredGrid.h>
#include <vtkLookupTable.h>
#include <vtkStructuredGridAppend.h>

#include <QuadratureActor.h>

namespace model
{

    QuadratureActor::QuadratureActor(vtkGenericOpenGLRenderWindow* const renWin,vtkRenderer* const ren,const Polycrystal<3>& poly,const DDtraitsIO& traitsIO) :
    /* init */ renderWindow(renWin)
    /* init */,renderer(ren)
    /* init */,mainLayout(new QGridLayout(this))
    /* init */,pkBox(new QGroupBox(tr("&PK forces")))
    /* init */,pkScaleEdit(new QLineEdit("10000"))
    /* init */,velocityBox(new QGroupBox(tr("&velocities")))
    /* init */,velocityScaleEdit(new QLineEdit("100"))
    /* init */,pkPolydata(vtkSmartPointer<vtkPolyData>::New())
    /* init */,pkGlyphs(vtkSmartPointer<vtkGlyph3D>::New())
    /* init */,pkMapper(vtkSmartPointer<vtkPolyDataMapper>::New())
    /* init */,pkActor(vtkSmartPointer<vtkActor>::New())
    /* init */,velocityPolydata(vtkSmartPointer<vtkPolyData>::New())
    /* init */,velocityGlyphs(vtkSmartPointer<vtkGlyph3D>::New())
    /* init */,velocityMapper(vtkSmartPointer<vtkPolyDataMapper>::New())
    /* init */,velocityActor(vtkSmartPointer<vtkActor>::New())
{
        
        pkBox->setCheckable(true);
        pkBox->setChecked(false);
        pkActor->SetVisibility(false);
        QVBoxLayout *pkLayout = new QVBoxLayout();
        pkLayout->addWidget(pkScaleEdit);
        pkBox->setLayout(pkLayout);

        velocityBox->setCheckable(true);
        velocityBox->setChecked(false);
        velocityActor->SetVisibility(false);
        QVBoxLayout *velocityLayout = new QVBoxLayout();
        velocityLayout->addWidget(velocityScaleEdit);
        velocityBox->setLayout(velocityLayout);

        
        mainLayout->addWidget(pkBox,0,0,1,1);
        mainLayout->addWidget(velocityBox,1,0,1,1);

        
        this->setLayout(mainLayout);

        connect(pkBox,SIGNAL(toggled(bool)), this, SLOT(modify()));
        connect(pkScaleEdit,SIGNAL(returnPressed()), this, SLOT(modify()));
        connect(velocityBox,SIGNAL(toggled(bool)), this, SLOT(modify()));
        connect(velocityScaleEdit,SIGNAL(returnPressed()), this, SLOT(modify()));

        
        pkGlyphs->SetInputData(pkPolydata);
        pkGlyphs->SetSourceConnection(vtkSmartPointer<vtkArrowSource>::New()->GetOutputPort());
        pkGlyphs->ScalingOn();
        pkGlyphs->SetScaleModeToScaleByVector();
        pkGlyphs->OrientOn();
        pkGlyphs->ClampingOff();
        pkGlyphs->SetVectorModeToUseVector();
        pkGlyphs->SetIndexModeToOff();
        pkMapper->SetInputConnection(pkGlyphs->GetOutputPort());
        pkMapper->ScalarVisibilityOff();
        pkActor->SetMapper(pkMapper);
        pkActor->GetProperty()->SetColor(0.0, 0.0, 1.0); //(R,G,B)
        renderer->AddActor(pkActor);
    
        velocityGlyphs->SetInputData(velocityPolydata);
        velocityGlyphs->SetSourceConnection(vtkSmartPointer<vtkArrowSource>::New()->GetOutputPort());
        velocityGlyphs->ScalingOn();
        velocityGlyphs->SetScaleModeToScaleByVector();
        velocityGlyphs->OrientOn();
        velocityGlyphs->ClampingOff();
        velocityGlyphs->SetVectorModeToUseVector();
        velocityGlyphs->SetIndexModeToOff();
        velocityMapper->SetInputConnection(velocityGlyphs->GetOutputPort());
        velocityMapper->ScalarVisibilityOff();
        velocityActor->SetMapper(velocityMapper);
        velocityActor->GetProperty()->SetColor(0.0, 1.0, 0.0); //(R,G,B)
        renderer->AddActor(velocityActor);

    }

    void QuadratureActor::updateConfiguration(const DDauxIO<3>& auxIO)
    {
        
        std::cout<<"Updating QuadraturePoints..."<<std::flush;
        const auto t0= std::chrono::system_clock::now();
        
        vtkSmartPointer<vtkPoints> nodePoints = vtkSmartPointer<vtkPoints>::New();
        
        vtkSmartPointer<vtkDoubleArray> pkForces(vtkSmartPointer<vtkDoubleArray>::New());
        pkForces->SetNumberOfComponents(3);
        
        vtkSmartPointer<vtkDoubleArray> velocities(vtkSmartPointer<vtkDoubleArray>::New());
        velocities->SetNumberOfComponents(3);
        
        for(const auto& qp : auxIO.quadraturePoints())
        {
            nodePoints->InsertNextPoint(qp.r.data());
            pkForces->InsertNextTuple(qp.pkForce.data()); // arrow vector
            velocities->InsertNextTuple(qp.glideVelocity.data()); // arrow vector
        }
        
        pkPolydata->SetPoints(nodePoints);
        pkPolydata->GetPointData()->SetVectors(pkForces);
        pkPolydata->Modified();

        velocityPolydata->SetPoints(nodePoints);
        velocityPolydata->GetPointData()->SetVectors(velocities);
        velocityPolydata->Modified();

        std::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
    }

    void QuadratureActor::modify()
    {
        pkActor->SetVisibility(pkBox->isChecked());
        const double pkScaling(std::atof(pkScaleEdit->text() .toStdString().c_str()));
        pkGlyphs->SetScaleFactor(pkScaling);

        velocityActor->SetVisibility(velocityBox->isChecked());
        const double velocityScaling(std::atof(velocityScaleEdit->text() .toStdString().c_str()));
        velocityGlyphs->SetScaleFactor(velocityScaling);

        
        renderWindow->Render();
    }


} // namespace model
#endif
