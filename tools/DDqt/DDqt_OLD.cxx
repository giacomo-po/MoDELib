/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2017 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


//#include "widget.h"

//#include <vtkAutoInit.h>
//VTK_MODULE_INIT(vtkRenderingOpenGL2);
//VTK_MODULE_INIT(vtkInteractionStyle);

//#include <QtOpenGL>
#include <QtWidgets/QApplication>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QWidget>
#include <QGridLayout>
#include <QPushButton>



#include <vtkSphereSource.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <QVTKOpenGLNativeWidget.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkNamedColors.h>
#include <vtkProperty.h>
#include <vtkNew.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkPropPicker.h>


#include <TextFileParser.h>
#include <SimplicialMesh.h>

#include <SimplicialMeshActor.h>

namespace model
{



class DDwidget : public QVTKOpenGLNativeWidget
{
//    Q_OBJECT
    
//    vtkSmartPointer<vtkGenericOpenGLRenderWindow> renderWindow;
    vtkSmartPointer<vtkRenderer> renderer;
    vtkSmartPointer<QVTKInteractor> renderWindowInteractor;
   vtkSmartPointer<vtkInteractorStyleTrackballCamera> style;

    const std::string polyFile;
    SimplicialMesh<3> mesh;
    SimplicialMeshActor meshActor;
    
public:
    DDwidget(QWidget *parent) :
    /* init */ QVTKOpenGLNativeWidget(parent)
//    /* init */,renderWindow(vtkSmartPointer<vtkGenericOpenGLRenderWindow>::New())
    /* init */,renderer(vtkSmartPointer<vtkRenderer>::New())
    /* init */,renderWindowInteractor(vtkSmartPointer<QVTKInteractor>::New())
    /* init */,style(vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New())
    /* init */,polyFile("/Users/giacomo/Documents/MoDELib2/tutorials/DislocationDynamics/periodicDomains/uniformLoadController/inputFiles/polycrystal.txt")
    /* init */,mesh(TextFileParser(polyFile).readString("meshFile",true),TextFileParser(polyFile).readMatrix<double>("A",3,3,true),TextFileParser(polyFile).readMatrix<double>("x0",1,3,true).transpose(),TextFileParser(polyFile).template readSet<int>("periodicFaceIDs",true))
    /* init */,meshActor(renderer,mesh)
    {
        renderer->SetBackground(1,1,1);
        this->renderWindow()->AddRenderer(renderer);

//        renderWindow->AddRenderer(ddRenderer);
        renderWindowInteractor->Initialize();

        renderWindowInteractor->SetRenderWindow(this->renderWindow());
        renderWindowInteractor->SetInteractorStyle(style);

//        this->setRenderWindow(renderWindow);
//        renderWindowInteractor->SetRenderWindow ( renderWindow );

//        this->show();
//        renderWindow->Render ();
//        renderWindowInteractor->Initialize();
        
//        this->setRenderWindow(renderWindow);
//        this->SetInteractor(renderWindowInteractor);
//        this->renderWindow()->AddRenderer(ddRenderer);
//        this->renderWindow()->SetInteractor(this->interactor());
//        this->SetInteractor(this->renderWindow()->interactor());
//
//        this->renderWindow()->GetInteractor()->SetInteractorStyle( style );
//        renderWindowInteractor->SetRenderWindow(this->renderWindow());
//        vtkNew<vtkInteractorStyleTrackballCamera> style;
//        vtkRenderWindow->SetInteractor(renderWindowInteractor);
//        renderWindowInteractor->SetInteractorStyle( style );
//        this->renderWindow()->GetInteractor()->Start();

//
//        renderWindowInteractor->SetInteractorStyle(style);
//        vtkNew<vtkNamedColors> colors;
//
//          // Create a sphere
//          vtkNew<vtkSphereSource> sphereSource;
//          sphereSource->SetCenter(0.0, 0.0, 0.0);
//          sphereSource->SetRadius(5.0);
//          // Make the surface smooth.
//          sphereSource->SetPhiResolution(100);
//          sphereSource->SetThetaResolution(100);
//
//          vtkNew<vtkPolyDataMapper> mapper;
//          mapper->SetInputConnection(sphereSource->GetOutputPort());
//
//          vtkNew<vtkActor> actor;
//          actor->SetMapper(mapper);
//          actor->GetProperty()->SetColor(colors->GetColor3d("Cornsilk").GetData());

//          vtkNew<vtkRenderer> renderer;
//            renderer->AddActor(actor);
//            renderer->SetBackground(colors->GetColor3d("DarkGreen").GetData());

//        vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
//        renderWindowInteractor->SetRenderWindow(this->renderWindow());
        
//        vtkSmartPointer<vtkInteractorStyleTrackballCamera> style = vtkSmartPointer<DDinteractionStyle>::New();
//    //      style->SetMyRenderWindow(renderWindow);
////        style->MyRenderWindow = this->renderWindow();
//        vtkRenderWindow->GetInteractor()->SetInteractorStyle( style );
////        style->SetRenderWindow(this->renderWindow());

//        vtkNew<vtkInteractorStyleTrackballCamera> style;
//
//        renderWindowInteractor->SetInteractorStyle(style);


        
 //       vtkNew<vtkGenericOpenGLRenderWindow> renderWindow;
//          vtkNew<vtkRenderWindow> renderWindow;
//          renderWindow->SetWindowName("Sphere");
//          renderWindow->AddRenderer(renderer);

//          renderWindowInteractor->SetRenderWindow(this->renderWindow());
//
////
//          renderWindow->Render();
//          renderWindowInteractor->Start();

    }
    
//    bool eventFilter(QObject *obj, QEvent *ev){
//      if (ev->type() == QEvent::KeyPress ) {
//        QKeyEvent* keyevent = dynamic_cast<QKeyEvent*>(ev);
//        if(keyevent->key() == Qt::Key_Q){
//          close();
//          // eventFilter intercepts the event.
//          return true;
//        }
//        else if(keyevent->key() == Qt::Key_R){
//          ResetCamera();
//          return true;
//        }
//        else{
//          vtkCamera* camera = renderer_->GetActiveCamera();
//          std::cout << tarckball_style->GetAutoAdjustCameraClippingRange()
//            << " : " << camera->GetClippingRange()[0]
//            << " : " << camera->GetClippingRange()[1] << std::endl;
//        }
//      }
//      // eventFilter doesn't intercept the event.
//      return false;
//    }

};

class MainWindow : public QMainWindow
{
    // https://kitware.github.io/vtk-examples/site/Cxx/Qt/RenderWindowNoUiFile/
// https://kitware.github.io/vtk-examples/site/Cxx/Qt/SideBySideRenderWindowsQt/
//https://vtk.org/doc/nightly/html/classQVTKOpenGLNativeWidget.html

//    vtkSmartPointer<vtkGenericOpenGLRenderWindow> renderWindow;
//    vtkSmartPointer<vtkRenderer> renderer;
//    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor;
//   vtkSmartPointer<vtkInteractorStyleTrackballCamera> style;

    
    QWidget *centralWidget;
    QGridLayout* mainLayout;
//    QVTKOpenGLNativeWidget* ddw1;
    DDwidget* ddw2;
    
    QPushButton* pushButton;

public:
    
    MainWindow() :
//    /* init */ renderWindow(vtkSmartPointer<vtkGenericOpenGLRenderWindow>::New())
//    /* init */,renderer(vtkSmartPointer<vtkRenderer>::New())
//    /* init */,renderWindowInteractor(vtkSmartPointer<vtkRenderWindowInteractor>::New())
//    /* init */,style(vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New())
    /* init */ centralWidget(new QWidget(this))
    /* init */,mainLayout(new QGridLayout(this))
//    /* init */,ddw1(new QVTKOpenGLNativeWidget(this))
    /* init */,ddw2(new DDwidget(this))
    /* init */,pushButton(new QPushButton(this))
    {
        resize(1920, 1080);

        
        pushButton->setText("Example");

//        vtkNew<vtkRenderer> renderer;
//        vtkNew<vtkGenericOpenGLRenderWindow> renderWindow;
//        renderWindow->AddRenderer(renderer);
////        ddw2->setRenderWindow(renderWindow);
//        ddw1->setRenderWindow(renderWindow);
////        vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
//        renderWindowInteractor->SetRenderWindow ( renderWindow );
//        ddw2->show();
//        renderWindow->Render ();
//        renderWindowInteractor->Initialize();
        
//        vtkSmartPointer<vtkInteractorStyleTrackballCamera> style;
//
//        ddw2->
//
        mainLayout->addWidget(ddw2,3,4,10,10);
        mainLayout->addWidget(pushButton,7,3,1,1);

        centralWidget->setLayout(mainLayout);
        
        setCentralWidget(centralWidget);

    }

};
}




int main(int argc, char *argv[])
{
    // needed to ensure appropriate OpenGL context is created for VTK rendering.
    //QSurfaceFormat::setDefaultFormat(QVTKOpenGLWidget::defaultFormat());
    QSurfaceFormat::setDefaultFormat(QVTKOpenGLNativeWidget::defaultFormat());
    QApplication a(argc, argv);
    model::MainWindow window;
    window.show();

    return a.exec();
}

