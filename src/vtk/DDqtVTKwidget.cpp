/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DDqtVTKwidget_cpp_
#define model_DDqtVTKwidget_cpp_

#include <QFileDialog>

#include <DDqtVTKwidget.h>


namespace model
{
    
    DDqtVTKwidget::DDqtVTKwidget(QWidget*) :
    /* init */ mainLayout(new QGridLayout(this))
    /* init */,tabWidget(new QTabWidget(this))
    /* init */,openglWidget(new QVTKOpenGLStereoWidget(this))
    /* init */,renderWindow(vtkSmartPointer<vtkGenericOpenGLRenderWindow>::New())
    /* init */,renderer(vtkSmartPointer<vtkRenderer>::New())
    /* init */,ddBase(getWorkingDir())
    /* init */,workingDirLabel(new QLabel(QString::fromStdString(ddBase.simulationParameters.traitsIO.simulationFolder)))
    /* init */,meshActor(new SimplicialMeshActor(renderWindow,renderer,ddBase.mesh))
    /* init */,ddConfigVtk(new DDconfigVtk(renderWindow,renderer,openglWidget,ddBase))
    {
        renderer->SetBackground(1,1,1);
        renderWindow->AddRenderer(renderer);
        renderWindow->SetMultiSamples(4);
        
        openglWidget->setRenderWindow(renderWindow);

        
        tabWidget->addTab(ddConfigVtk, tr(std::string("Config").c_str()));
        tabWidget->addTab(meshActor, tr(std::string("Mesh").c_str()));

        mainLayout->addWidget(workingDirLabel,0,0,1,2);
        mainLayout->addWidget(tabWidget,1,0,1,1);
        mainLayout->addWidget(openglWidget,1,1,1,1);
        mainLayout->setColumnStretch(0, 3);
        mainLayout->setColumnStretch(1, 7);
        this->setLayout(mainLayout);
    }
    

    std::string DDqtVTKwidget::getWorkingDir() 
    {
        return QFileDialog::getExistingDirectory(this, tr("Open Directory"),
                                                             "/Users/giacomo/Documents/MoDELib2/tutorials/DislocationDynamics",
                                                             QFileDialog::ShowDirsOnly
                                                             | QFileDialog::DontUseNativeDialog ).toStdString();

    
    }


} // namespace model

#endif







